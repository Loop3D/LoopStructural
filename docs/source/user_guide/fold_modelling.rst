Fold modelling
==============

Fold interpolators constrain a scalar-field gradient using a fold frame, a fold
axis and a fitted limb-rotation profile. ``DFI`` uses linear tetrahedra;
``P2Fold`` uses quadratic tetrahedra. Select either when creating a folded
foliation or a folded fold frame. The latter uses the fold interpolator for its
first coordinate; its other two coordinates still use FDI.

Using P2Fold
------------

``P2Fold`` with its default ``gauss4`` sampling is useful when sparse orientation
data leave the fold constraints to control the field between observations. Its
four interior samples constrain the gradient throughout a quadratic tetrahedron.
In the Bermagui validation it reduced median fold-direction error with fewer
unknowns than DFI for both S2 and S0. This comparison measures agreement with a
shared fold event, rather than independent geological truth, and does not imply
that P2Fold is always the better interpolator.

Given an existing model and fold frame::

    folded = model.create_and_add_folded_foliation(
        "s0",
        fold_frame=frame,
        interpolatortype="P2Fold",
        nelements=1000,
        buffer=0.0,
        limb_wl=None,
        fold_weights={"sample_rule": "gauss4"},
        solver="cg",
        solver_kwargs={"maxiter": 20000, "rtol": 1e-8},
    )
    values = folded.evaluate_value(points)  # triggers the lazy build
    print(folded.builder.interpolator.dof)

The same ``interpolatortype`` and ``fold_weights`` options can be passed to
``create_and_add_folded_fold_frame``. The bundled synthetic example is
:ref:`sphx_glr__auto_examples_2_fold_plot_3_p2_fold.py`.

Sampling and mesh size
----------------------

``fold_weights["sample_rule"]`` accepts ``"gauss4"`` (default) or ``"centroid"``.
The library weights each sample residual by ``w * sqrt(q_i)``, with quadrature
weights summing to one in each tetrahedron. Adding more samples therefore does
not multiply the total normalised operator weight by the sample count.

A single centroid sample does not adequately control a varying P2 gradient when
regularisation is disabled. The unregularised Bermagui centroid runs had about
60 degrees median error and unstable field scales, even when CG converged.
Use gauss4 for this configuration; increasing mesh size did not rescue centroid.
This result does not establish the behaviour of every regularised centroid model.

Start with a small ``nelements`` request and inspect the realised unknown and
tetrahedron counts::

    interpolator = folded.builder.interpolator
    print(interpolator.support.n_elements, interpolator.dof)

The current support factory treats ``nelements`` as a structured-grid cell
budget, then splits each cell into five tetrahedra. P2 adds nodes on the edges.
It is not a literal tetrahedron count, and equal requests are not equal cost
between DFI and P2Fold. For the Bermagui box with ``buffer=0.0``:

.. list-table:: Realised P2 meshes
   :header-rows: 1
   :widths: 25 35 40

   * - Requested nelements
     - Tetrahedra
     - Unknowns
   * - 1,000
     - 5,850
     - 9,711
   * - 10,000
     - 56,420
     - 85,673
   * - 100,000
     - 525,690
     - 765,347

The largest gauss4 S0 case took about 93 seconds in ``solve_system`` and used
about 6.6 GiB peak resident memory on the validation machine. These are single
measurements, not runtime or memory guarantees. The
`Bermagui handoff <https://github.com/Loop3D/LoopStructural/blob/master/agent/HANDOFF_p2_fold.md>`_
records the corrected comparison summary and limitations. The scripts in
``scripts/validation`` reproduce the detailed results and figures locally.

Interpreting the result
-----------------------

At the smaller budget, S2 gauss4 had 4.69 degrees median error with 9,711 unknowns,
compared with DFI's 6.15 degrees at 12,992 unknowns. S0 had 11.14 degrees at 85,673
unknowns, compared with 20.63 degrees at 112,560 unknowns. S0 still had 19.4 percent
of inner probes above 30 degrees; inspect the error distribution and fitted
profiles rather than relying only on the median.

The scalar-field scale is arbitrary. The norm constraint need not be met, so the
field should not be interpreted as a distance function. A wavelength expressed
in field units is not portable between independently fitted models. Inspect each
rotation profile with its observations. The Bermagui S0 validation uses
``limb_wl=None`` (automatic fitting); starts between 2 and 5 were poor local
minima in the earlier experiments. Other datasets require their own checks.

Remaining issues include buffer application at both model-box and support
construction, lower-face NaNs in the linear tetrahedral support (F2), amplitude
damping (F3), and rotation-profile fit instability. Use the same explicit buffer
when comparing engines. ``buffer=0.0`` keeps this validation inside the model box;
boundary queries can return NaN. These issues are separate follow-up work.
