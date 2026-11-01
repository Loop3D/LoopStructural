.. _native_solver:

Native Solver
=============
The native solver is a Rust solver. It solves the same discrete system as the Python solvers, but it uses less memory.
Use it when a model has many nodes and the Python solvers run out of memory.

The native solver uses the package ``loop_interpolation_core``. Install it with the optional extra:

.. code-block:: bash

    pip install loop-interpolation[native]

Select the native solver for an interpolator before the interpolator is set up:

.. code-block:: python

    interpolator = InterpolatorFactory.create_interpolator(..., solver="native")

The solver can also be set with ``InterpolatorBuilder.use_solver`` or with ``solve_system(solver=...)``.
Set the solver before ``setup``. When the solver is set before setup, the Python interpolator does not build its sparse rows.
This saves Python memory. When the solver is set after setup, the rows are already built. The result is still correct, but the memory saving is lost.
In this case, the solver logs one warning.

Solver options
--------------
The native solver accepts these keyword arguments in ``solver_kwargs``:

.. list-table:: Native solver options
    :widths: 25 75
    :header-rows: 1

    * - Option
      - Meaning
    * - tol
      - Relative residual tolerance. The default is 1e-8.
    * - maxiter
      - Maximum number of inner conjugate gradient iterations. The default is 5000.
    * - precision
      - Type of the solver vectors. ``"f64"`` (double, the default) or ``"f32"`` (single).
    * - budget_bytes
      - Memory budget in bytes. The solver stops with an error when the estimate is above the budget.
    * - rho
      - ADMM penalty. The default is 1.0.
    * - inner_max_iter
      - Maximum number of inner iterations of ADMM. The default is 3000.
    * - inner_tol
      - Tolerance of the inner iterations of ADMM. The default is 1e-6.
    * - admm_abs_tol
      - Absolute tolerance of the ADMM outer loop. The default is 1e-4.
    * - admm_rel_tol
      - Relative tolerance of the ADMM outer loop. The default is 1e-3.
    * - ridge_factor
      - Ridge factor. The solver adds ``ridge_factor**2`` to the diagonal of the normal equations. 0 removes the ridge term. The default is the ``ridge_factor`` of the interpolator (1e-8).
    * - add_ridge_regularisation
      - ``False`` removes the ridge term for this solve. The default is the setting of the interpolator.
    * - warm_start
      - ``True`` starts the solver from the last field of the interpolator, not from 0. If the last field is not finite, the solver starts from 0. The default is ``False``.
    * - x0
      - Start value, one value for each node, or a function that takes the support and returns one. It has priority over ``warm_start``.

The value of ``tol`` in ``solver_kwargs`` replaces the ``tol`` argument of the solve.
A start value needs 8 more bytes for each node until the solve starts. The memory estimate includes these bytes.
Use ``warm_start`` when you solve again after a small change, for example a new weight. The solve then needs fewer iterations.
Any other keyword is rejected with a ``NotImplementedError``. The message lists the accepted keywords.

Checking that the solve converged
---------------------------------
The solver stops after ``maxiter`` iterations even when it has not reached ``tol``. When this happens it logs a warning
(``The native solver did not converge``), sets ``interpolator.up_to_date`` to ``False``, and leaves the unfinished field in the interpolator.
Nothing is raised. On a test problem with about 55 000 unknowns the default limit of 5000 iterations was not enough: the field was far
from the Python solution, although the residual in the warning was close to ``tol``. Raising ``maxiter`` fixed it.
Check ``up_to_date`` after a native solve, and raise ``maxiter`` if it is ``False``.
The gallery example :doc:`/_auto_examples/5_interpolation/plot_6_native_solver` shows this, and the memory saving, on a worked problem.

Unsupported settings
--------------------
The native solver does not support these settings. When one of them is set, the solver raises ``NotImplementedError``. The message contains the name of the setting.

* A region mask (``region``).
* ``use_regularisation_weight_scale``.
* Constant-norm iterations (``constant_norm_iterations``) with a true value.
* ``interface_pair_mode`` with the value ``"all"`` or ``"chain"``. Use ``"star"``.
* ``interface_pair_fallback`` with the value ``"chain"``.
* Matrix-free regularisation (``regularisation_matrix_free``).
* Custom operators (``operators``).
* A 2D support.
* A rectilinear grid.
* A rotated support.

The native solver works on uniform 3D structured grids only.

Column scaling (``apply_scaling_matrix``) is ignored. The solver writes one information message to the log.
Column scaling changes the conditioning of the matrix. It does not change the solution.

Threads
-------
The Python wheel of the core uses several threads. The browser build uses one thread.
On a 100 by 100 by 100 grid, the wheel solved a test problem about 3 times faster than the Python matrix solver.
On one thread, the native solver was about 10 percent slower. The result does not depend on the number of threads (differences are in the last bits).

Memory estimate
---------------
The core estimates the memory of a solve before it starts. The estimate uses only the counts of the problem. It does not allocate memory.
The estimate has four parts, in bytes: ``vectors``, ``constraints``, ``working`` and ``total``.

The solver vectors dominate the estimate. For conjugate gradient in double precision, a solve needs about 40 bytes for each node (measured for a 100 by 100 by 100 grid).
This is true when each axis has 7 or more nodes and there is no directional regularisation. Then the solver applies the regularisation with one stencil for each class of node, and it needs no other buffer for each node.
For a smaller grid, or with directional regularisation, the solver uses the older operator, and a solve needs about 63 bytes for each node.
The ADMM solver, for problems with bounds, needs about the same memory. It keeps a few more values for each bound, not for each node. The constraint records need more memory in addition to these values.

With ``precision="f32"``, the vector part is half of the value for double precision.

The core module ``loop_interpolation_core`` returns this estimate from the method ``Problem.estimate(settings)``. The ``loop_interpolation`` package does not call this method.
Use the estimate to choose a grid size before a long solve, or to set ``budget_bytes``.
