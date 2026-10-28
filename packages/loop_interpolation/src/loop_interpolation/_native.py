"""Adapter between a ``FiniteDifferenceInterpolator`` and the native core.

The native core is the Rust package ``loop_interpolation_core``. This module
reads the stored constraint data of the interpolator, builds a core
``Problem``, solves it, and writes the field back to ``interpolator.c``.

Two functions are public:

* ``check_supported(interpolator, setup_kwargs, solver_kwargs)`` raises
  ``NotImplementedError`` for a setting that the core cannot represent.
  The message contains the name of the setting.
* ``solve(interpolator, tol, solver_kwargs) -> bool`` runs the solve.

The module does not import the core at import time. The core is imported
inside ``solve``, so this module loads when the core is not installed.

The adapter reads the raw data with the ``get_*_constraints`` methods. It
never reads the sparse rows.

MAPPING from the Python data to the core meaning
------------------------------------------------
The core applies a weight verbatim to a row that it has scaled to unit L2
norm. These are the equivalent Python weights (the same mapping as
``packages/loop_interpolation_core/fixtures/make_fixture.py``).
Python column layouts are x y z | data | weight.

* Grid. ``counts`` is ``support.nsteps`` (number of NODES for each axis).
  The node order is x fastest, as in ``interpolator.c``. A rotated support
  is not supported.
* Regularisation. ``dxx dyy dzz`` are the same in Python and core. The
  cross-derivative weights are divided by 4 inside ``get_operators``, so
  the core weight is ``dxy / 4`` (same for ``dxz`` and ``dyz``).
* Border. ``dx dy dz`` are the same in Python and core.
* Value. Weight is ``cpw * column 4``.
* Gradient (direction only). Weight is ``gpw``. The weight column is NOT
  used by Python.
* Normal (vector with length). Weight is ``npw / 3 * column 6``. When the
  data has no weight column, the column weight is 1.
* Tangent. Weight is ``tpw``. The weight column is NOT used by Python.
* Interface. Weight is ``ipw`` (one value for all points). The pairing is
  "star": each point is paired with the first point of its id group.
* Value inequality. Columns x y z | lower | upper. The weight is 1.
* Inequality pairs. Columns x y z | rock id. A pair (a, b) gives the bound
  ``lower <= f(point of a) - f(point of b) <= upper`` for all point
  combinations of the two rock ids. The weight is 1. The first id of a
  pair is the upper side. A row with a nan or inf position is dropped.
  A pair can have ``n_a * n_b`` rows. The setup keyword
  ``inequality_pair_max_rows`` (default 100000) limits the rows of each
  pair. When ``n_a * n_b`` is more than the limit, the adapter uses that
  many random combinations and writes one warning. The random generator is
  ``numpy.random.default_rng`` with the seed ``inequality_pair_seed``
  (default 0), so two runs with the same input give the same rows. The
  adapter picks the indices and does not build the full product. It sends
  the rows to the core in chunks of 100000.
* Directional regularisation. The direction field is evaluated at the
  regularisation sample points (the nodes) and passed with its weight.
* A constraint family with weight 0 is not sent to the core. The core
  raises an error for a weight of 0.
* ``ridge_factor`` has no setting in the core. It is not mapped.
* Column scaling (``apply_scaling_matrix``) is ignored for any value. It
  changes the conditioning of the matrix, not the minimiser. One info line
  says that the adapter skips it.
"""

from __future__ import annotations

from time import perf_counter

import numpy as np
from loop_common.logging import get_logger as getLogger

logger = getLogger(__name__)

# Keys of solver_kwargs that the adapter accepts. Value: key of the core settings.
_SOLVER_KEYS = {
    "tol": "tol",
    "maxiter": "max_iter",
    "precision": "precision",
    "budget_bytes": "budget_bytes",
    "rho": "rho",
    "inner_max_iter": "inner_max_iter",
    "inner_tol": "inner_tol",
    "admm_abs_tol": "admm_abs_tol",
    "admm_rel_tol": "admm_rel_tol",
}

# Inequality pairs: default cap of rows for each pair, seed of the sampling,
# and rows for each call of add_pair_bound.
_PAIR_MAX_ROWS = 100_000
_PAIR_SEED = 0
_PAIR_CHUNK = 100_000

_REGULARISATION_KEYS = ("dxx", "dyy", "dzz", "dxy", "dxz", "dyz")
_BORDER_KEYS = ("dx", "dy", "dz")

_INSTALL_HINT = (
    "solver='native' needs the package loop_interpolation_core. "
    "Install it with: pip install loop-interpolation[native]"
)


def _unsupported(setting: str, reason: str = "") -> NotImplementedError:
    text = f"solver='native' does not support the setting '{setting}'"
    if reason:
        text += f": {reason}"
    return NotImplementedError(text)


def check_supported(interpolator, setup_kwargs, solver_kwargs) -> None:
    """Raise ``NotImplementedError`` for a setting that the core cannot represent.

    Parameters
    ----------
    interpolator : FiniteDifferenceInterpolator
    setup_kwargs : dict or None
        The keywords of ``setup_interpolator``.
    solver_kwargs : dict or None
        The keywords of the solver.
    """
    setup_kwargs = dict(setup_kwargs or {})
    solver_kwargs = dict(solver_kwargs or {})

    # Only the plain finite difference class is known to match the core.
    # A subclass must set the class attribute supports_native = True in its
    # own body. The flag is not inherited, so a subclass of a supported
    # subclass must set it again.
    # Fold, P1 and P2 interpolators add terms that the core cannot represent.
    from ._finite_difference_interpolator import FiniteDifferenceInterpolator

    cls = type(interpolator)
    if cls is not FiniteDifferenceInterpolator and not cls.__dict__.get("supports_native", False):
        raise _unsupported(
            "interpolator class",
            f"{cls.__name__} is not supported; the core works with "
            "FiniteDifferenceInterpolator on a uniform structured grid",
        )

    # Sparse rows that the user adds are not in the raw data. The core cannot
    # read them. setup_interpolator records the names of the rows that it
    # made itself. Each other name in "constraints" is an extra row. This
    # is true after a native setup and after a matrix setup. The rows of
    # "inequality_feature" are always extra.
    if "inequality_feature" in getattr(interpolator, "ineq_constraints", {}):
        raise _unsupported(
            "add_inequality_feature", "the core reads only the raw data, not sparse rows"
        )
    setup_names = getattr(interpolator, "_setup_constraint_names", None) or frozenset()
    if set(getattr(interpolator, "constraints", None) or {}) - setup_names:
        raise _unsupported(
            "add_constraints_to_least_squares",
            "the core reads only the raw data, not sparse rows",
        )

    support = interpolator.support
    if getattr(support, "dimension", None) != 3:
        raise _unsupported("2D support", "the core works on 3D structured grids")
    if type(support).__name__ == "RectilinearGrid" or not all(
        hasattr(support, name) for name in ("origin", "nsteps", "step_vector")
    ):
        raise _unsupported("support type", "the core works on uniform structured grids")
    rotation = getattr(support, "rotation_xy", None)
    if rotation is not None and not np.allclose(np.asarray(rotation, dtype=float), np.eye(3)):
        raise _unsupported("rotation_xy", "rotated supports are not supported")

    if setup_kwargs.get("region") is not None:
        raise _unsupported("region", "a region mask is not supported")
    if "operators" in setup_kwargs:
        raise _unsupported("operators", "custom operators are not supported")
    if setup_kwargs.get("use_regularisation_weight_scale", False) or getattr(
        interpolator, "use_regularisation_weight_scale", False
    ):
        raise _unsupported("use_regularisation_weight_scale")
    if solver_kwargs.get("constant_norm_iterations"):
        raise _unsupported("constant_norm_iterations")
    if setup_kwargs.get("regularisation_matrix_free", False):
        raise _unsupported("regularisation_matrix_free")

    mode = str(
        setup_kwargs.get(
            "interface_pair_mode", getattr(interpolator, "interface_pair_mode", "auto")
        )
    ).lower()
    if mode in ("all", "chain"):
        raise _unsupported("interface_pair_mode", f"mode '{mode}' is not supported; use 'star'")
    if mode == "auto":
        fallback = str(
            setup_kwargs.get(
                "interface_pair_fallback",
                getattr(interpolator, "interface_pair_fallback", "star"),
            )
        ).lower()
        if fallback == "chain":
            raise _unsupported("interface_pair_fallback", "the fallback 'chain' is not supported")

    for key in solver_kwargs:
        if key not in _SOLVER_KEYS and key != "constant_norm_iterations":
            raise NotImplementedError(
                f"solver='native' does not support the solver keyword '{key}'. "
                f"Supported keywords: {', '.join(_SOLVER_KEYS)}"
            )

    if getattr(interpolator, "apply_scaling_matrix", False):
        logger.info("solver='native': column scaling (apply_scaling_matrix) is skipped.")


def _core_settings(tol, solver_kwargs) -> dict:
    """Make the settings dict of the core from ``tol`` and ``solver_kwargs``."""
    settings = {}
    if tol is not None:
        settings["tol"] = float(tol)
    for key, value in solver_kwargs.items():
        if key == "constant_norm_iterations":
            continue
        if value is None:
            continue
        settings[_SOLVER_KEYS[key]] = value
    return settings


def _flat(a) -> np.ndarray:
    return np.ascontiguousarray(np.asarray(a, dtype=np.float64).ravel())


def _weights(w, n: int) -> np.ndarray:
    """Return a weight array with ``n`` values from a scalar or an array."""
    return np.ascontiguousarray(np.broadcast_to(np.asarray(w, dtype=np.float64), (n,)).copy())


def _rows(array) -> np.ndarray:
    """Return a 2D float array. An empty or missing array gives 0 rows."""
    if array is None:
        return np.zeros((0, 0))
    a = np.asarray(array, dtype=np.float64)
    if a.ndim != 2 or a.shape[0] == 0:
        return np.zeros((0, a.shape[1] if a.ndim == 2 else 0))
    return a


def _weight(interpolator, key: str, default: float = 1.0):
    return interpolator.interpolation_weights.get(key, default)


# How the matrix path treats bad data rows (FD interpolator, shared helpers):
# - A row with a nan or inf in the position or the data is dropped.
# - A negative weight is kept. The solver squares the weight, so only the size counts.
# - A weight of 0 gives a row with no effect. A nan weight is dropped.
# - A gradient or tangent vector of length 0 gives no constraint (dropped).
# - A normal vector of length 0 is kept. It constrains the gradient to 0.
# The native path does the same. A negative weight becomes its absolute value,
# because the core accepts only weights that are more than 0.
def _add_weighted(add, points, data, weights, drop_zero_norm=False, name="data") -> None:
    """Call ``add(points, data, weights)`` for the rows that the matrix path keeps."""
    data2 = data.reshape(data.shape[0], -1)
    keep = np.isfinite(weights) & (weights != 0)
    keep &= np.isfinite(points).all(axis=1) & np.isfinite(data2).all(axis=1)
    if drop_zero_norm:
        keep &= np.linalg.norm(data2, axis=1) > 0
    n_dropped = int(keep.size - np.count_nonzero(keep))
    if n_dropped:
        logger.debug("native: %d of %d %s rows dropped (bad data)", n_dropped, keep.size, name)
    if not np.any(keep):
        return
    add(
        _flat(points[keep]),
        _flat(data[keep]),
        _flat(np.abs(weights[keep])),
    )


def _add_data(problem, interpolator) -> None:
    """Send the point constraints to the core."""
    # value: weight = cpw * column 4
    value = _rows(interpolator.get_value_constraints())
    if value.shape[0] > 0:
        w = _weights(_weight(interpolator, "cpw"), value.shape[0]) * value[:, 4]
        _add_weighted(problem.add_value, value[:, :3], value[:, 3], w, name="value")

    # gradient: weight = gpw (the weight column is not used)
    gradient = _rows(interpolator.get_gradient_constraints())
    if gradient.shape[0] > 0:
        w = _weights(_weight(interpolator, "gpw"), gradient.shape[0])
        _add_weighted(
            problem.add_gradient,
            gradient[:, :3],
            gradient[:, 3:6],
            w,
            drop_zero_norm=True,
            name="gradient",
        )

    # normal: weight = npw / 3 * column 6 (1 when there is no column)
    normal = _rows(interpolator.get_norm_constraints())
    if normal.shape[0] > 0:
        w = _weights(_weight(interpolator, "npw") / 3.0, normal.shape[0])
        if normal.shape[1] > 6:
            w = w * normal[:, 6]
        _add_weighted(problem.add_normal, normal[:, :3], normal[:, 3:6], w, name="normal")

    # tangent: weight = tpw (the weight column is not used)
    tangent = _rows(interpolator.get_tangent_constraints())
    if tangent.shape[0] > 0:
        w = _weights(_weight(interpolator, "tpw"), tangent.shape[0])
        _add_weighted(
            problem.add_tangent,
            tangent[:, :3],
            tangent[:, 3:6],
            w,
            drop_zero_norm=True,
            name="tangent",
        )

    # interface: one weight (ipw)
    interface = _rows(interpolator.get_interface_constraints())
    if interface.shape[0] > 1:
        ipw = np.asarray(_weight(interpolator, "ipw"), dtype=np.float64)
        if ipw.size > 1 and not np.all(ipw == ipw.flat[0]):
            raise _unsupported(
                "interface per-point weights", "the core takes one weight for all interface points"
            )
        ipw = float(ipw.flat[0])
        # Same rules as the matrix path:
        # - A row with a nan id or a nan or inf position is dropped.
        # - A group with fewer than two points gives no pair and is dropped.
        # - A negative weight becomes its absolute value (see _add_weighted).
        # - A weight of 0 or a nan weight gives no constraint.
        valid = np.isfinite(interface[:, :4]).all(axis=1)
        ids, counts = np.unique(interface[valid, 3], return_counts=True)
        valid &= np.isin(interface[:, 3], ids[counts > 1])
        n_dropped = int(valid.size - np.count_nonzero(valid))
        if n_dropped:
            logger.debug(
                "native: %d of %d interface rows dropped (bad id, bad position or single point)",
                n_dropped,
                valid.size,
            )
        if np.isfinite(ipw) and ipw != 0 and np.any(valid):
            problem.add_interface(_flat(interface[valid, :3]), _flat(interface[valid, 3]), abs(ipw))


def _add_inequalities(problem, interpolator, setup_kwargs) -> None:
    """Send the value bounds and the pair bounds to the core."""
    bounds = _rows(interpolator.get_inequality_value_constraints())
    if bounds.shape[0] > 0:
        problem.add_value_bound(
            _flat(bounds[:, :3]),
            _flat(bounds[:, 3]),
            _flat(bounds[:, 4]),
            _flat(np.ones(bounds.shape[0])),
        )

    pairs_data = _rows(interpolator.get_inequality_pairs_constraints())
    if pairs_data.shape[0] == 0:
        return
    upper_bound = setup_kwargs.get("inequality_pair_upper_bound", np.finfo(float).eps)
    lower_bound = setup_kwargs.get("inequality_pair_lower_bound", -np.inf)
    pairs = setup_kwargs.get("inequality_pairs", None)
    max_rows = setup_kwargs.get("inequality_pair_max_rows", _PAIR_MAX_ROWS)
    if int(max_rows) < 1:
        raise ValueError(f"inequality_pair_max_rows must be 1 or more, got {max_rows}")
    max_rows = int(max_rows)
    rng = np.random.default_rng(setup_kwargs.get("inequality_pair_seed", _PAIR_SEED))
    # The core stops on a nan or inf position. Drop these rows.
    finite = np.isfinite(pairs_data[:, :3]).all(axis=1)
    if not finite.all():
        logger.debug("native: %d inequality pair rows dropped (bad position)", int((~finite).sum()))
        pairs_data = pairs_data[finite]
    ids = pairs_data[:, 3]
    if not pairs:
        # Same default as the matrix solver: each unordered pair of rock ids.
        unique = np.unique(ids)
        pairs = [(a, b) for i, a in enumerate(unique) for b in unique[i + 1 :]]
    for upper_id, lower_id in pairs:
        points_a = pairs_data[ids == upper_id, :3]
        points_b = pairs_data[ids == lower_id, :3]
        n_a, n_b = points_a.shape[0], points_b.shape[0]
        if n_a == 0 or n_b == 0:
            continue
        # One row for each combination of a point of each rock id. A pair with
        # more than max_rows combinations gets a random sample of max_rows rows.
        total = n_a * n_b
        if total > max_rows:
            logger.warning(
                "native: inequality pair (%s, %s) has %d point combinations. "
                "Only %d random combinations are used (inequality_pair_max_rows).",
                upper_id,
                lower_id,
                total,
                max_rows,
            )
            flat = _sample_indices(rng, total, max_rows)
        else:
            flat = np.arange(total)
        # Send the rows in chunks. Index arithmetic makes the chunk, so the
        # full set of combinations is never in memory.
        for start in range(0, flat.size, _PAIR_CHUNK):
            k = flat[start : start + _PAIR_CHUNK]
            n = k.size
            problem.add_pair_bound(
                _flat(points_a[k // n_b]),
                _flat(points_b[k % n_b]),
                _flat(np.full(n, lower_bound)),
                _flat(np.full(n, upper_bound)),
                _flat(np.ones(n)),
            )


def _sample_indices(rng, total: int, count: int) -> np.ndarray:
    """Return ``count`` different random integers in ``[0, total)``, sorted.

    ``rng.choice(total, replace=False)`` can allocate ``total`` values, and
    ``total`` can be very large. This function allocates about ``count``
    values. When ``total`` is more than 4 times ``count``, it draws with
    replacement and removes the repeats until ``count`` values remain.
    Then, a repeat is rare and a round or two is enough.
    """
    if total <= 4 * count:
        return np.sort(rng.choice(total, size=count, replace=False))
    picked = np.unique(rng.integers(0, total, size=count))
    while picked.size < count:
        extra = rng.integers(0, total, size=count - picked.size)
        picked = np.unique(np.concatenate([picked, extra]))
    return picked


def _add_regularisation(problem, interpolator, setup_kwargs) -> None:
    """Send the regularisation, border and directional weights to the core."""
    reg = {}
    for key in _REGULARISATION_KEYS:
        value = float(_weight(interpolator, key))
        reg[key] = value / 4.0 if key in ("dxy", "dxz", "dyz") else value
    problem.set_regularisation(reg)
    problem.set_border_regularisation(
        {key: float(_weight(interpolator, key)) for key in _BORDER_KEYS}
    )

    config = interpolator.resolve_regularisation_config(
        regularisation=setup_kwargs.get("regularisation", None),
        directional_regularisation=setup_kwargs.get("directional_regularisation", None),
    )
    if len(config.directional) == 0:
        return
    sample_points = np.asarray(interpolator.get_regularisation_sample_points(), dtype=float)
    for term in config.directional:
        if term.weight == 0:
            continue
        vectors = term.direction(sample_points) if callable(term.direction) else term.direction
        vectors = np.asarray(vectors, dtype=float)
        if vectors.shape != sample_points.shape:
            logger.warning(
                "%s: directional regularisation vectors must have shape %s, got %s. Skipping.",
                term.name,
                sample_points.shape,
                vectors.shape,
            )
            continue
        problem.add_directional_regularisation(_flat(vectors), float(term.weight))


def solve(interpolator, tol, solver_kwargs) -> bool:
    """Solve the interpolation problem with the native core.

    Parameters
    ----------
    interpolator : FiniteDifferenceInterpolator
    tol : float or None
        Relative residual tolerance. ``None`` keeps the default of the core.
    solver_kwargs : dict or None
        Keys: tol, maxiter, precision, budget_bytes, rho, inner_max_iter,
        inner_tol, admm_abs_tol, admm_rel_tol.

    Returns
    -------
    bool
        True when the solve converged.

    Notes
    -----
    ``interpolator.c`` changes only when the solve converged. When the solve
    does not converge, ``c`` keeps its old value. The partial field is in
    ``interpolator.latest_native_report["field"]``. A converged report has no
    ``field`` key. The matrix solvers (for example ``cg``) differ here. They
    write the last iterate to ``c`` also when they do not converge.

    The defaults of the core are ``tol=1e-8`` and ``maxiter=5000``. The
    ``cg`` matrix solver uses ``rtol=1e-8``. Its ``maxiter`` is the scipy
    default, 10 times the number of nodes. So the native solver can stop
    earlier on large grids. Set ``maxiter`` to get the same limit.
    """
    started = perf_counter()
    solver_kwargs = dict(solver_kwargs or {})
    setup_kwargs = dict(getattr(interpolator, "_setup_kwargs", None) or {})

    # Check the settings first. These errors must not need the core.
    check_supported(interpolator, setup_kwargs, solver_kwargs)

    try:
        import loop_interpolation_core as core
    except ImportError as err:
        raise ImportError(_INSTALL_HINT) from err

    support = interpolator.support
    problem = core.Problem(
        np.asarray(support.origin, dtype=np.float64),
        np.asarray(support.nsteps, dtype=np.int64),
        np.asarray(support.step_vector, dtype=np.float64),
    )
    _add_regularisation(problem, interpolator, setup_kwargs)
    _add_data(problem, interpolator)
    _add_inequalities(problem, interpolator, setup_kwargs)

    # tol: the solver_kwargs value wins over the tol argument.
    settings = _core_settings(solver_kwargs.get("tol", tol), solver_kwargs)
    field, report = problem.solve(settings)

    # The core can return float32. ``c`` is always float64.
    field = np.asarray(field, dtype=np.float64)
    if field.shape[0] != support.n_nodes:
        raise RuntimeError(
            f"The native solver returned {field.shape[0]} values. "
            f"The support has {support.n_nodes} nodes."
        )

    converged = bool(report.get("converged", False))
    if converged:
        interpolator.c = field
    else:
        # Keep the old ``c``. Give the partial field in the report.
        report["field"] = field
    interpolator.latest_native_report = report
    if not converged:
        logger.warning(
            "The native solver did not converge: %s iterations, residual %s.",
            report.get("iterations"),
            report.get("residual"),
        )
    for message in report.get("warnings", []):
        logger.warning("native solver: %s", message)

    total = perf_counter() - started
    interpolator.latest_solve_timing = {
        "solver": "native",
        "total": total,
        "total_seconds": total,
        "iterations": report.get("iterations"),
        "up_to_date": converged,
    }
    return converged
