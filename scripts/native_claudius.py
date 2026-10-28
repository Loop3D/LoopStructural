"""Build the claudius model with the native (Rust) solver.

The script builds the stratigraphy of the claudius demo data set with the
finite difference interpolator and ``solver="native"``. It can also build
the same model with the Python ``cg`` solver and compare the two fields.

Needs the native core::

    pip install loop-interpolation[native]

or, in a checkout::

    maturin develop --release -m packages/loop_interpolation_core/Cargo.toml

Examples::

    python scripts/native_claudius.py
    python scripts/native_claudius.py --nelements 1000000 --solver native
    python scripts/native_claudius.py --plot
    python scripts/native_claudius.py --solver native --screenshot claudius.png
"""

from __future__ import annotations

import argparse
import logging
import time

import numpy as np
import pandas as pd

# Import the adapter now, so that its logger exists when main sets the level.
from loop_interpolation import _native  # noqa: F401

import LoopStructural
from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius
from LoopStructural.modelling.features import GeologicalFeature

# Scalar field values of the horizons in the claudius data.
UNIT_VALUES = [0, 60, 250, 330, 600]


def build_model(
    data: pd.DataFrame, bb: np.ndarray, solver: str, nelements: int, tol: float, maxiter: int
) -> tuple[GeologicalModel, GeologicalFeature, float]:
    """Build the claudius stratigraphy with one solver.

    Return the model, the feature and the build time in seconds.
    """
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.data = data
    for i in range(len(UNIT_VALUES) - 1):
        model.stratigraphic_column.add_unit(
            f"unit_{i}", thickness=UNIT_VALUES[i + 1] - UNIT_VALUES[i], id=i
        )
    model.stratigraphic_column.group_mapping["Group_0"] = "strati"

    if solver == "native":
        # tol is the relative residual of the core. maxiter is the limit of
        # the CG iterations (the core default is 5000).
        solver_kwargs = {"tol": tol, "maxiter": maxiter}
    else:
        solver_kwargs = {"rtol": tol, "maxiter": maxiter}
    feature = model.create_and_add_foliation(
        "strati",
        interpolatortype="FDI",
        nelements=nelements,
        buffer=0.3,
        solver=solver,
        solver_kwargs=solver_kwargs,
    )

    started = time.perf_counter()
    model.update(progressbar=False)
    seconds = time.perf_counter() - started
    return model, feature, seconds


def report(name: str, feature: GeologicalFeature, seconds: float) -> None:
    """Print the grid, the time and the solver result of one build."""
    interpolator = feature.builder.interpolator
    nsteps = interpolator.support.nsteps
    print(f"\n{name}")
    print(f"  grid        {nsteps[0]} x {nsteps[1]} x {nsteps[2]} = {interpolator.dof} nodes")
    print(f"  build time  {seconds:.2f} s")
    native = getattr(interpolator, "latest_native_report", None)
    if name == "native" and native is not None:
        print(f"  converged   {native['converged']}")
        print(f"  iterations  {native['iterations']}")
        print(f"  residual    {native['residual']:.2e}")
        print(
            f"  core peak   {native['peak_bytes'] / 1e6:.1f} MB "
            f"({native['peak_bytes'] / interpolator.dof:.1f} bytes per node)"
        )
        for warning in native.get("warnings", []):
            print(f"  warning     {warning}")
    else:
        timing = interpolator.get_last_solve_timing()
        print(f"  converged   {timing.get('up_to_date')}")
        print(f"  solve time  {timing.get('solve_seconds', float('nan')):.2f} s")


def main() -> None:
    """Parse the arguments, build the models and print the results."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--nelements", type=int, default=200_000, help="grid cells (default 200000)")
    parser.add_argument(
        "--solver",
        choices=["native", "cg", "both"],
        default="both",
        help="solver to use; 'both' also compares the fields (default both)",
    )
    parser.add_argument("--tol", type=float, default=1e-8, help="relative residual (default 1e-8)")
    parser.add_argument("--maxiter", type=int, default=20_000, help="CG iteration limit (default 20000)")
    parser.add_argument("--plot", action="store_true", help="show the native model surfaces")
    parser.add_argument(
        "--screenshot", metavar="PATH", help="save the native model surfaces to a PNG file"
    )
    parser.add_argument("--verbose", action="store_true", help="show the INFO log messages")
    args = parser.parse_args()
    LoopStructural.setLogging("info" if args.verbose else "warning")
    # loop_interpolation sets a level on each of its loggers, so set them all.
    level = logging.INFO if args.verbose else logging.WARNING
    for name in list(logging.root.manager.loggerDict):
        if name.startswith("loop_interpolation"):
            logging.getLogger(name).setLevel(level)

    data, bb = load_claudius()
    solvers = ["native", "cg"] if args.solver == "both" else [args.solver]

    models = {}
    for solver in solvers:
        model, feature, seconds = build_model(
            data, bb, solver, args.nelements, args.tol, args.maxiter
        )
        report(solver, feature, seconds)
        models[solver] = model

    if len(models) == 2:
        # Compare the two scalar fields on a regular grid in the model box.
        xyz = models["native"].regular_grid(nsteps=(40, 40, 40), shuffle=False)
        native = models["native"].evaluate_feature_value("strati", xyz)
        cg = models["cg"].evaluate_feature_value("strati", xyz)
        field_range = float(np.nanmax(cg) - np.nanmin(cg))
        difference = float(np.nanmax(np.abs(native - cg)))
        print("\nnative vs cg")
        print(f"  largest difference  {difference:.3e}")
        print(f"  field range         {field_range:.3e}")
        print(f"  relative            {difference / field_range:.2e}")

    if args.plot or args.screenshot:
        from LoopStructural.visualisation import Loop3DView

        model = models.get("native", next(iter(models.values())))
        # Without --plot the view is drawn off screen, for the screenshot only.
        viewer = Loop3DView(model, off_screen=not args.plot)
        viewer.plot_model_surfaces(cmap="tab20")
        viewer.add_mesh(model.bounding_box.vtk().outline(), color="black", name="outline")
        if args.screenshot:
            viewer.screenshot(args.screenshot, window_size=(1400, 1000))
            print(f"\nsaved {args.screenshot}")
        if args.plot:
            viewer.display()


if __name__ == "__main__":
    main()
