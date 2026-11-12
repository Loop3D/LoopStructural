6. Fault modelling scenarios
----------------------------
A set of synthetic fault scenarios, from one planar fault to fault networks,
faults and unconformities, and families of models for fault uncertainty.
The scenarios are taken from the fault meshing benchmark suite
(``fault_mesh_benchmarks``).

Every example uses the unit cube :math:`[0, 1]^3` with z up, and the same
workflow:

1. Define each fault by observations of its surface (points with normals)
   and by its kinematics (displacement, slip vector, axes and displacement
   profile).
2. Build the faults, youngest first when one fault displaces another.
3. Make synthetic stratigraphic observations that agree with the fault
   kinematics: restore random points through the faults, and give each
   point the value of a known stratigraphy in its restored position.
4. Interpolate the faulted stratigraphy, and look at it in 3D and on
   vertical sections.

Because the stratigraphy before faulting is known, the faults are the only
source of offset, and some examples compare the model with the exact
solution.

.. _fault-scenario-future-work:

Known limits and future work
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The examples show some things that LoopStructural cannot do yet, or does
only with a workaround. These are recorded for a future specification
(``specs/backlog/fault-scenario-gaps.md``):

* **Displacement profiles.** Only ``"BaseFault"`` and ``"BaseFault3D"`` have
  names. A constant (step) profile and a one-sided (hanging wall only)
  profile must be built from the classes in the private module
  ``_fault_function``, with ``scale=1.0`` (the default ``scale`` of
  ``FaultDisplacement`` is 0.5). The decay distance cannot be set for each
  wall separately without a custom ``CubicFunction``.
* **Tip line shape.** The profile is a product of a ``gy`` profile and a
  ``gz`` profile, so the tip line is a rectangle in the fault plane, not an
  ellipse.
* **Frame scaling.** ``gx`` is scaled by ``minor_axis / 2`` (a full width),
  but ``gy`` and ``gz`` by ``intermediate_axis`` and ``major_axis`` (radii).
* **Faults displaced by other faults.** A series point is restored twice
  through the younger fault (double restoration), and the slip direction of
  a displaced fault is slow to calculate. The examples use a workaround,
  ``fault_network_workaround``.
* **Fault age relations.** There is no way to give two faults the same age
  (conjugate faults that end at their intersection). A constant profile
  moves the whole hanging-wall half-space, including other faults.
* **Unconformities do not cut faults.** The surface of a fault that is older
  than an unconformity continues into the younger rocks.
* **Forward modelling.** ``FaultSegment.apply_to_points(reverse=True)`` does
  not move points forward. The examples invert the restoration by
  iteration.
* **Displacement queries.** ``GeologicalModel.evaluate_fault_displacements``
  returns the profile value, not the displacement.
* **Surfaces across faults.** Isosurfaces are extracted on a regular grid,
  so they are smeared across faults, and the horizons of a narrow fault
  block need a fine output grid over the whole model (no local
  refinement). The interpolated field itself is correct.
* **Fault uncertainty.** There is no sampler for fault parameters.
