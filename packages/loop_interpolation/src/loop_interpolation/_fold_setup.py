"""Shared setup helpers for fold-aware interpolators."""


def setup_with_fold_constraints(
    *,
    fold,
    kwargs: dict,
    interpolator_name: str,
    base_setup,
    add_fold_constraints,
    finalize_report,
):
    """Run common fold interpolator setup sequence.

    Parameters
    ----------
    fold
        Fold event instance or None.
    kwargs : dict
        Setup kwargs passed to setup_interpolator.
    interpolator_name : str
        Name used in error messages.
    base_setup : callable
        Parent setup_interpolator callable.
    add_fold_constraints : callable
        Method that applies fold-specific constraints.
    finalize_report : callable
        Callable that returns final diagnostics report.
    """
    if fold is None:
        raise RuntimeError(
            f"{interpolator_name}: no fold event set. Assign self.fold before calling setup_interpolator."
        )

    setup_kwargs = dict(kwargs)
    fold_weights = setup_kwargs.pop("fold_weights", {})
    base_setup(**setup_kwargs)
    add_fold_constraints(**fold_weights)
    return finalize_report()
