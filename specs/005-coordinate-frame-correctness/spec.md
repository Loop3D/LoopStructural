# Feature spec: Coordinate and fault-frame correctness

**Feature**: `005-coordinate-frame-correctness` | **Created**: 2026-11-03 | **Status**: S1 merged (pull request #392); S2 in #396; S3 in #397; S4 planned; S5 and S6 need planning

**Input**: Fix the shared world-to-local coordinate defects and the fault-frame interpolation defects found in the two investigations.

## Context

Geologists supply observations and fault geometry in world coordinates. A mesh uses local coordinates to keep numerical calculations near its origin. Some paths mix these spaces. The resulting geology can depend on the map origin, the mesh rotation, or whether the field was saved and loaded.

The investigations found defects in reference constraints, fold constraints, fault construction, support bounds, persistence, data display, grid sampling, and diagnostics. Adding more observations can hide some failures. It does not repair the coordinate contract.

The smallest useful result is a correctly located and directed fault frame, built through the normal public API, with its shared constraint paths fixed. Other stories protect reuse, display, diagnostics, and the remaining audited utility defects. Each story must be independently testable and reversible.

This spec distinguishes three spaces:

| Space | Meaning |
|---|---|
| World | Physical point locations and directions in the model's map coordinates |
| Mesh local | Point locations and directions used by one interpolation support |
| Geological frame | Scalar field values such as the three fault coordinates |

Fault restoration changes a physical point. Its input and output remain in world coordinates. It is not a projection into mesh-local coordinates.

The current public contract already accepts world-coordinate observations and queries. This spec completes that contract. It does not restore the previous convention of shifting observations at model ingestion.

### Investigation evidence

The investigations used LoopStructural 1.7.0 at source commit `91326cce91c02e3d75a902f035f189214babfb07`. The checkout has since advanced. Planning must check each finding against the implementation selected for this work.

| Reproduction | Observed result |
|---|---|
| Translated equality region | 0 rows installed; 168 rows expected |
| Translated analytic fold | 60 degree orientation error |
| Direct interpolator save and reload | Valid world values became NaN for FDI, PLI, and P2 |
| Rotated support construction | 4 of 8 model-box corners outside the support |
| Display of a world observation | `(102,203,304)` displayed at `(202,403,604)` |
| Gradient diagnostic | Stored constraint entries changed by up to 5.322 |
| Assign an unchanged support maximum | Support shortened by 1 cell on each axis |
| Isolated sampling correction | Fold row residuals below `3e-14` for the three tested assemblers |

The fault investigation also isolated an incorrectly projected explicit centre, an overwritten buffer, and an ambiguous extent direction. Correcting their combined effect recovered the planar reference frame to numerical precision with an accurate solve.

These results motivate the requirements. They are not acceptance results for production fixes. Acceptance tests must exercise production paths without the investigation's replacement samplers.

## User scenarios and tests

### User story 1: Build a fault frame at its supplied location (priority: P1)

A geologist supplies fault observations, a centre, directions, and dimensions. The resulting frame describes that geometry at the supplied world location.

**Why this priority**: this is the original modelling failure. A displacement calculated from a misplaced or tilted frame describes the wrong geology.

**Independent test**: build the planar reference fault through the public API. Compare its values and gradients with the analytic frame. Repeat after translation.

**Acceptance scenarios**:

1. **Given** an explicit world centre, **when** the fault is built, **then** all generated observations use that centre without an extra projection.
2. **Given** observations whose inferred centre equals the explicit reference centre, **when** both faults are built, **then** their frame values and directions agree within SC-003.
3. **Given** a translated copy of the same geometry, **when** the translated queries are evaluated, **then** scalar values agree within SC-003. Gradients retain their world directions.
4. **Given** each R5 buffer request, **when** support construction completes, **then** its actual bounds match SC-005. Different requests produce different bounds.
5. **Given** caller-forced mesh geometry, **when** the factory completes, **then** a later default mesh operation does not overwrite it.
6. **Given** a centre, normal, and slip direction, **when** the frame is built, **then** its centre extent direction is consistent with those inputs. An explicit planar mode is not required.
7. **Given** curved fault observations with the same local centre geometry, **when** the frame is built, **then** its extent direction can vary away from the centre. Generated extent information does not impose a global plane.

### User story 2: Constrain features at the intended physical points (priority: P1)

A modeller relates one feature to another through orthogonality, equality, inequalities, or fold geometry. Each relationship is evaluated at the intended world locations.

**Why this priority**: these relationships share the defective boundary. A fault-only correction would leave other frames and folded features exposed.

**Independent test**: constrain analytic fields with spatially varying reference gradients and world-coordinate regions. Compare selected locations, directions, indices, targets, and row residuals.

**Acceptance scenarios**:

1. **Given** local support samples, **when** a reference feature is evaluated, **then** it receives their corresponding world locations.
2. **Given** a world reference gradient, **when** its constraint is assembled, **then** its direction is expressed in the target support's local basis.
3. **Given** a world-coordinate region, **when** equality constraints are selected, **then** every intended node is selected exactly once. Its local index remains associated with the sampled world value.
4. **Given** a feature that bounds another field, **when** inequalities are built, **then** each bound equals the reference value at that node's world location.
5. **Given** fold orientation, axis, normalisation, or directional regularisation, **when** its geometry is sampled, **then** both point locations and direction bases follow the same contract.
6. **Given** a rotated computational frame, **when** an exactly representable affine field is evaluated, **then** its world values and gradients agree with the unrotated control.
7. **Given** support expansion through fault restoration, **when** restored bounds are calculated, **then** restoration receives world points. The bounds contain their projections into the target support.

### User story 3: Save, copy, and display the same physical model (priority: P1)

A modeller saves an interpolator, reloads it, copies its configuration, or displays its data. These operations retain the physical meaning of its coordinate frame.

**Why this priority**: a correct solve is not useful if reuse or display changes its location.

**Independent test**: round-trip a translated and rotated solved field. Compare world evaluations, constraint locations, displayed observations, and saved transform metadata.

**Acceptance scenarios**:

1. **Given** a solved field and a nonidentity transform, **when** it is saved and loaded, **then** world values and gradients satisfy SC-006.
2. **Given** saved local constraint data, **when** the file is loaded, **then** loading does not apply the world-to-local transform to those arrays again.
3. **Given** an API that returns an empty copy, **when** the copy is created, **then** it retains the coordinate frame. It need not contain the source observations or coefficients.
4. **Given** an observation already stored in world coordinates, **when** it is displayed or exported as observation geometry, **then** its physical location is unchanged.
5. **Given** a world-coordinate fault centre or intrusion observation, **when** a downstream consumer uses it, **then** no extra inverse projection changes its location.
6. **Given** a legacy file with no recoverable transform or explicit caller declaration, **when** loading is requested on the 1.x line, **then** the file loads with the identity frame that older releases assumed, and one warning that end users see by default names the choices and says that the next major release raises an error. From the next major release, loading fails instead. The error identifies the missing transform and the required caller action. No partial field is published.
7. **Given** that legacy file and a caller-supplied transform or explicit identity declaration, **when** loading is retried, **then** its world evaluations satisfy SC-006. The loader does not guess identity without a warning.
8. **Given** a newly saved interpolator, **when** its payload is inspected, **then** it declares its format version, transform, and constraint coordinate space. These declarations survive 2 save/load cycles.
9. **Given** a populated interpolator, **when** its transform is changed, **then** the documented FR-020 policy applies. Rejection preserves prior state; so does the rejection of each later use of the stale state. Acceptance preserves world evaluation after rebuilding. Restoration caches that exist (spec 003 adds them) are populated before a separate check for stale use.
10. **Given** a completed stage and its documented rollback target, **when** the stage is reversed, **then** old-format R1 fields remain usable. An unsupported new-format field fails without changing valid state.
11. **Given** a stable method that this work deprecates, **when** a caller uses it, **then** it returns the same result as before and emits one `DeprecationWarning` that names its replacement. No call inside the library emits that warning.

### User story 4: Trust support geometry and diagnostics (priority: P2)

A modeller checks residuals and samples a field on a grid. These checks refer to the same physical locations as the data and exported geometry.

**Why this priority**: wrong diagnostics can conceal a modelling defect and can currently change the constraints being checked.

**Independent test**: check an analytic field, resize a support to unchanged bounds, and compare cell-centre samples with their associated grid cells.

**Acceptance scenarios**:

1. **Given** a rotated model box, **when** a support is made or resized, **then** it encloses all projected corners.
2. **Given** unchanged requested bounds, **when** the bounds are assigned again, **then** the support does not lose a cell or change its geometry.
3. **Given** a world sampling grid, **when** cell centres are generated, **then** each point lies at the centre of its associated cell.
4. **Given** a gradient or value diagnostic, **when** it runs, **then** it leaves observations, vectors, weights, coefficients, and model state unchanged.
5. **Given** a row that cannot be installed, **when** assembly completes, **then** its exclusion is reported with a reason and a count.
6. **Given** an invalid transformation matrix, **when** it is supplied as a rotation, **then** it is rejected before valid model state changes.

### User story 5: Verify displacement and auxiliary transformations independently (priority: P2)

A developer compares frame interpolation, displacement, and coordinate utilities with independent analytic references. Each test identifies which operation failed.

**Why this priority**: deform-then-restore agreement alone can conceal a wrong frame or fail because of an independent integration error.

**Independent test**: compare an interpolated planar frame with an analytic frame using the same displacement evaluation. Separately verify point-transform round trips.

**Acceptance scenarios**:

1. **Given** an interpolated frame and its analytic reference, **when** displacement is evaluated away from discontinuities, **then** their difference satisfies SC-004.
2. **Given** points on a fault discontinuity, **when** error results are reported, **then** they are counted separately from the smooth-region result.
3. **Given** an invalid or zero gradient, **when** displacement is evaluated, **then** the existing protection against division by zero remains effective.
4. **Given** an ordinary point cloud and a supported auxiliary transformation, **when** forward and inverse transformations run, **then** their round trip satisfies SC-010.
5. **Given** analytic evidence that a golden result contains a coordinate defect, **when** its correction is reviewed, **then** only justified affected results change. Unaffected baselines and acceptance tolerances remain unchanged.
6. **Given** the backend coverage matrix, **when** acceptance checks run, **then** every required cell has a result. Missing optional dependencies have a reason. No mandatory Python cell is skipped.

### Edge cases

- Nonzero, negative, and map-scale origins; independent translations on the three axes.
- Identity, 30 degree, and 90 degree computational rotations; rotations that exchange or reverse projected bounds.
- Two-dimensional and three-dimensional point-transform operations; a single point and a point array.
- A target and reference feature with different supports and different local transforms.
- Curved reference gradients; a restricted reference support; finite and invalid samples mixed in one batch.
- Zero-weight relationships, zero-length sampled directions, NaN values, and empty selected regions.
- A support with nonuniform axis spacing; repeated resizing; different meanings of node and cell counts at an export boundary.
- Explicit and inferred fault centres; nonorthogonal input fault directions; curved observations; conflicting extent observations.
- Buffers of 0, 0.2, and 0.5; forced mesh geometry; constraints on and outside support boundaries.
- A transform changed after data, coefficients, or restored-point caches already exist.
- Legacy files without transform metadata load with the identity frame and one warning on the 1.x line, unless the caller supplies a transform or explicitly declares identity. From the next major release they fail until the caller does so. Recoverable metadata is used without guessing.
- Malformed metadata and unsupported future format versions fail without publishing a partial field or changing prior valid state.
- Conflicting caller transform information and saved metadata fail with an error naming the conflict. An explicit identity declaration is not an implicit default.
- A diagnostic with no valid vectors; a constant scalar field; residual normalisation with zero field range.
- Invalid rotation matrices, including scale, shear, singular matrices, nonfinite entries, and wrong dimensions.
- A backend not installed or a feature it does not support. This spec does not add support by silently ignoring the feature.

## Requirements

### Functional requirements

- **FR-001**: Public model and geological-feature operations MUST accept world-coordinate points and return gradient and displacement directions in world coordinates.
- **FR-002**: Model observations, generated fault geometry, and builder observations MUST retain their world-coordinate locations until the interpolation boundary.
- **FR-003**: Each boundary MUST convert point locations exactly once. Translation MUST NOT be applied to a direction, gradient, scalar value, or geological-frame coordinate.
- **FR-004**: Feature-facing reference callbacks MUST receive world-coordinate points. Internal local-coordinate callbacks MUST have an explicit and documented contract.
- **FR-005**: Directions sampled from a reference feature MUST be converted into the target support's basis before local assembly. The reference support's transform MUST NOT substitute for the target transform.
- **FR-006**: Equality region selection and scalar targets MUST use the same world locations, with the corresponding target-local node indices retained.
- **FR-007**: Feature inequality bounds MUST use reference values sampled at the target nodes' world locations.
- **FR-008**: Fold geometry, directional regularisation, quadrature samples, masks, and normal-sign alignment MUST use consistent location and direction spaces.
- **FR-009**: Fault restoration used for support expansion MUST operate on world points. The resulting local support MUST contain the restored world samples after projection.
- **FR-010**: Explicit and inferred fault centres MUST have the same world-coordinate meaning throughout construction and downstream use. Inference from R2's symmetric observations MUST preserve all 3 translated coordinates, including elevation.
- **FR-011**: Supplied centre, normal, and slip directions MUST define a consistent local extent direction at the centre without requiring an explicit planar mode. Use the existing frame handedness: for R2 the extent direction is `+y`. Curved geometry MUST remain free to vary away from the centre.
- **FR-012**: Named fault-buffer and forced-geometry settings MUST retain their requested effect through every factory and builder operation. The support MUST match the requested bounds within the declared mesh-rounding allowance. Containment alone is insufficient.
- **FR-013**: Support construction and resizing MUST enclose the projection of every corner of the requested world box.
- **FR-014**: Node and cell counts MUST remain consistent during support resizing. Reapplying unchanged bounds MUST preserve support geometry.
- **FR-015**: Cell-centre sampling MUST agree with the world geometry and cell association of the returned grid. No centre may lie beyond its associated cell because of mixed bases or endpoint spacing.
- **FR-016**: Interpolator persistence MUST retain its transform and the declared coordinate space of saved constraints. A format version MUST identify the new metadata.
- **FR-017**: Loading saved local arrays MUST NOT apply the public world-input conversion to them a second time.
- **FR-018**: APIs that copy an interpolator configuration MUST retain its coordinate frame, including APIs whose documented result is empty. An API that has no such copy (the piecewise linear interpolator has no `copy`) stays unsupported and keeps its documented rejection.
- **FR-019**: A legacy file with no recoverable transform MUST require a caller-supplied transform (the local origin and the local rotation, both) or an explicit identity declaration. A partial transform is not enough. Missing, malformed, unsupported, or conflicting transform information MUST fail with an actionable error before a partial field is published. The loader MUST NOT silently assume identity. On the 1.x line, a legacy file with no choice at all is the one exception to the error: it loads with the identity frame that older releases assumed, and it emits one `FutureWarning` that names the choices and says that the next major release raises an error. The code that allows this MUST carry a note to remove it at the next major release (C-018).
- **FR-020**: Changes to a transform MUST NOT leave dependent constraints, coefficients, or caches reusable in their previous coordinate frame. Rebuilding the dependent state, rejecting the change before mutation, or rejecting every later use of the stale state are acceptable documented policies. In each policy the stored constraints, coefficients, and caches stay unchanged until the caller rebuilds them.
- **FR-021**: A matrix accepted as a rotation MUST be finite, have the correct dimensions, and represent a proper rigid rotation. Scale and shear MUST NOT silently use rigid gradient-transformation rules.
- **FR-022**: Displayed observations, fault-centre geometry, and intrusion classification MUST use their stored world locations without an extra inverse projection.
- **FR-023**: Gradient and value diagnostics MUST compare quantities in one coordinate basis and MUST NOT mutate the state they inspect.
- **FR-024**: Verification MUST report physical data residuals and assembled constraint residuals separately. Excluded rows MUST be counted and explained.
- **FR-025**: The audited auxiliary point-transformation utilities MUST support forward/inverse round trips for their documented point-cloud shapes. Duplicate public entry points MUST retain their compatible import paths.
- **FR-026**: Fault-frame validation MUST check scalar values, gradient magnitudes, gradient directions, constraint residuals, and displacement against analytic references.
- **FR-027**: Integration error, fault discontinuities, and outside-support behaviour MUST be assessed separately from smooth frame-interpolation error. Existing zero-gradient protection MUST remain effective.
- **FR-028**: The public model API, return shapes, supported import paths, and inputs with known coordinate meaning MUST be preserved. On the 1.x line, a legacy file keeps its old load path (FR-019, C-018), so no stage breaks the public API there. The change of that path to an error happens at the next major release, which the maintainer plans. `ROADMAP.md` and `COMPAT.md` record the plan and the removal of the fallback. Any required public API breaking change MUST follow the constitution and roadmap decision process.
- **FR-029**: Analytic evidence MAY justify an explicitly reviewed correction to an affected golden result. Each change MUST identify the defect and its expected physical effect. Unaffected baselines and acceptance tolerances MUST remain unchanged. Tests MUST NOT preserve a demonstrated defect merely because an old result contains it.
- **FR-030**: Each required backend and capability in the coverage matrix MUST satisfy the coordinate contract through its integration boundary. Optional exclusions MUST identify the unavailable dependency or unsupported capability. A new native, browser, or RBF solver MUST NOT become a prerequisite for existing Python fixes.
- **FR-031**: Production-path regression tests MUST reproduce each defect before its fix and verify it after the fix. Replacement investigation samplers MUST NOT be the acceptance path.
- **FR-032**: The work MUST ship in independently testable stages. Each stage MUST have a documented rollback and MUST preserve compatibility entry points.
- **FR-033**: A stable public method that this work stops using inside the library MAY be deprecated if its behaviour and signature stay unchanged. A deprecation MUST emit a runtime `DeprecationWarning` that names the replacement. It MUST be logged in `COMPAT.md` with its removal rule. No call inside the library MAY trigger it. The first deprecations are `GeologicalModel.rescale` (local to world; the replacement is the bounding-box reprojection) and `GeologicalModel.scale` (world to local; the replacement is the bounding-box projection). Both go in the same release, with one lifecycle.

### Ownership and boundaries

The existing architecture reference is ROADMAP.md. It has no machine-readable `arch` block. This spec therefore has no `arch-impact` block or generated architecture check.

| Existing module or package | Responsibility in this spec |
|---|---|
| Common geometry | Owns the world/local rigid-transform meaning, transformed bounds, and geometry conventions |
| Interpolation | Owns the world-input boundary, local stored constraints and assembly, reference sampling, and interpolator persistence |
| Geological modelling | Owns world observations, generated fault geometry, geological-frame semantics, restoration, and feature consumers |
| Grid and display consumers | Preserve the physical geometry and associations produced by those owners |

This extends existing components. It does not propose a new service, process, or standalone package.

The main scenario crosses a public API boundary and a saved-file boundary:

1. The caller sends world observations and fault geometry to geological modelling.
2. Geological modelling sends world observations and the frame transform to interpolation.
3. Interpolation uses common geometry to obtain local points and directions for assembly.
4. Interpolation requests reference feature geometry at corresponding world locations.
5. The local solver returns field coefficients to interpolation.
6. The caller requests values, gradients, diagnostics, or displacement at world locations.
7. Evaluation returns world-directed results; persistence records the transform and data-space meaning for reuse.

### Related specs and compatibility

No unfinished solver implementation is a prerequisite for the existing Python corrections. A documentation cascade to specs 003 and 004 is a prerequisite for starting the plan of spec 005.

| Spec | Relationship |
|---|---|
| 001-low-memory-interpolator | The low-level core accepts an unrotated computational grid. Its caller continues to own external coordinate conversion. This spec must make that ownership explicit at the integration boundary. |
| 002-rbf-interpolator | Its public value/gradient integration must preserve world coordinates. Internal RBF conditioning remains the RBF solver's responsibility and must not be applied twice. The new solver is not a dependency. |
| 003-fault-evaluation-performance | Before planning 005, amend FR-006 and its dependent preservation criteria with C-004's reviewed correctness exception. Keep performance targets and unaffected baselines unchanged. Include transform changes in restored-point cache invalidation. |
| 004-matrix-free-directional-constraints | Before planning 005, amend AC-009 and AC-010 with C-004's reviewed correctness exception for affected results and fingerprints. Keep bounds and unaffected fingerprints unchanged. Clarify conversion of world directions into the computational basis without changing the low-level core contract. |

This review changes only spec 005. It does not perform the required cascade. The cascade must apply the accepted exception and resolve the conflicts before planning starts. C-004 is settled and does not need another policy decision.

## Success criteria

The values below are proposed acceptance targets for review. They are not claims that the unfixed implementation meets them.

Define `L` as the longest model-box side, `W` as the maximum absolute world-coordinate component, and `eps` as double-precision machine epsilon. Define a position tolerance `Tp = 1e-10*L + 32*eps*max(W,L)`. Define `F` as the larger of 1 and the maximum absolute reference scalar value. Unless stated otherwise, numeric controls use double precision.

### Reference fixtures

Use world anchor `O` and world offsets `(x',y',z') = point - O`. Offsets are physical world-axis components, not mesh-local coordinates. Test anchors `(0,0,0)`, `(100,200,300)`, `(-100,-200,-300)`, and `(500000,6500000,100)`.

For each computational rotation, use the world model minimum as the local-frame origin. Rotate 0, 30, or 90 degrees about z. Physical-rotation controls rotate geometry, query points, and directions together. They enclose the rotated world corners in the model box.

#### R1: Affine field and lifecycle

- World box: `O` to `O + (8,10,12)`.
- Field: `f = 2 + 0.2*x' - 0.3*y' + 0.4*z'`. Gradient: `(0.2,-0.3,0.4)`.
- Evaluation and persistence support: 5 by 6 by 7 regular cells, or the corresponding supported tetrahedral or quadratic mesh over the same bounds.
- Exact evaluation control: set the nodal field from the analytic world field. This isolates coordinate conversion from solver accuracy.
- Fit control: 27 unit-weight value observations on the product of `x'={1,4,7}`, `y'={1,5,9}`, and `z'={1,6,11}`.
- Fit support: 216 requested grid cells, or the corresponding tetrahedral mesh. Add no extra observations to hide a misplaced datum.
- World query offsets: `(2.2,3.7,4.3)` and `(6.1,8.2,9.1)`.
- Use the SC-003 solve settings for fit controls. For a constant-norm backend, set the norm target to `sqrt(0.29)`, the analytic gradient norm.

#### R2: Planar fault and displacement

- World model box: `O - (10,10,10)` to `O + (10,10,10)`.
- Centre: `O`. Repeat with no explicit centre; the symmetric observations must infer `O`, including its z coordinate.
- Normal: `(1,0,0)`. Slip direction: `(0,0,-1)`. Extent direction: `(0,1,0)`.
- Major, minor, and intermediate axes: `20`, `10`, and `20` world-length units.
- Input: 25 unit-weight observations on `x'=0`, with both `y'` and `z'` in `{-8,-4,0,4,8}`. Each has surface-coordinate value 0.
- Supply no extra slip-coordinate or extent-coordinate observations. The production geometry builder must generate the required local direction information.
- Use the normal-vector geometry mode, buffer `0.2`, and no forced mesh geometry. Disable extra trace and dip anisotropy terms.
- Request 20000 grid cells for FDI and 1000 grid cells for PLI. These requests are not actual tetrahedron counts.
- Use curvature weight 1, each frame orthogonality weight 1, and orthogonality sampling step 10.
- Expected fields: `gx=x'/20`, `gy=-z'/20`, and `gz=y'/20`.
- Value and gradient probes: the 6859 points on `{-9,-8,...,9}^3`.
- Displacement section: `y'=0`, with `x'` and `z'` on `{-9,-8.5,...,9}`. There are 1369 points.
- Smooth-region mask: `|x'| >= 0.5`, leaving 1332 points. Report the 37 excluded fault-plane points separately.
- Prescribed slip: 3 world-length units. The displacement comparison uses the library displacement paths, with the default fault profile and 20 integration steps (C-019).
- The scalar comparison uses the library scalar path. It uses no quadrature rule (C-019).
- The vector comparison uses the library displacement path (`FaultSegment.apply_to_points`) for both frames, with the same algorithm for both frames. It must not add a new library displacement law.
- Also check the library scalar displacement path (`FaultSegment.evaluate_displacement`) with the same function on both frames. Its comparison measures frame effects separately from vector integration.
- The analytic frame uses the 3 expected fields above and their exact gradients. If a displacement function needs coordinate length scales, use `(20,20,20)` for the analytic frame.

The R2 acceptance path builds the frame through production APIs. The library displacement paths consume that frame. They must not replace its builder or sampler. Include the reference helpers with the implementation's test fixtures; tests must not import an external investigation workspace.

#### R3: Reference relationships

- World box and initial regular cells: as in R1.
- Reference field: `h=x' + 0.05*y'^2 + z'`. Gradient: `(1,0.1*y',1)`.
- Equality region: `x' > 4`. For the unrotated FDI and PLI grids, select exactly 168 of 336 nodes.
- For a rotated grid or a quadratic mesh, compute the expected selection independently from the world-node locations. Do not reuse the unrotated count.
- Feature inequality bounds: the analytic value of `h` at each target node's world location.
- Full reference support: the entire test box. Restricted support: the closed world box from `O+(2.1,2.1,2.1)` to `O+(5.9,7.9,9.9)`.
- Samples outside restricted support return an invalid reference value or direction. Count exactly the corresponding excluded rows.
- Test barycentres and the supported centroid and four-point tetrahedral quadrature rules. Check region and mask callbacks at the same declared world samples.
- Fold controls cover orientation, axis, nonzero normalisation targets, normal-sign alignment, and each supported directional regularisation family.
- For an affine orthogonality control, use field `f=y'` and world direction `(1,0,0)`. Its dot-product target is 0 under every rigid transform.
- For a nonzero-target control, use that same field with world direction `(0,1,0)` and target 1. Rotation must preserve that target.

#### R4: Curved geometric compatibility

Use the box `O-(10,10,10)` to `O+(10,10,10)` and radius `R=40`. Define `X=R+x'` and `r=sqrt(X^2+y'^2)`.

| Coordinate | Analytic value | Analytic world gradient |
|---|---|---|
| Surface | `(r-R)/20` | `(X/(20*r), y'/(20*r), 0)` |
| Slip | `-z'/20` | `(0,0,-1/20)` |
| Extent | `R*atan2(y',X)/20` | `(-R*y'/(20*r^2), R*X/(20*r^2), 0)` |

Use these fields as analytic reference geometry through the production relationship APIs. Sample the 125 points on `{-6,-3,0,3,6}^3`. Include centre `O` and probes `O+(0,-6,0)` and `O+(0,6,0)`.

At the centre, supply normal `(1,0,0)` and slip direction `(0,0,-1)`. Their extent direction is `(0,1,0)`. The analytic extent directions at the two off-centre probes differ by more than 15 degrees.

Check automatic direction constraints introduced by the centre-direction correction against the analytic gradients at each constraint's own world location. Repeating the centre direction over the support must fail this check.

R4 checks geometric compatibility and sampling. It does not require an exact nonlinear solution on a linear mesh. Analytic reference fields are inputs to production relationships; replacement production samplers are prohibited.

#### R5: Buffer and forced geometry

- World model box: `O-(10,10,10)` to `O+(10,10,10)`.
- Centre, normal, and slip direction: as in R2. Repeat with the centre inferred from symmetric observations.
- Input: 25 unit-weight surface observations on `x'=0`, with `y'` and `z'` in `{-4,-2,0,2,4}`.
- Major, minor, and intermediate axes: `12`, `4`, and `12` world-length units. Generated observations lie inside the unbuffered model box.
- Force the base mesh geometry to the world model box. Request buffers `b=0`, `0.2`, and `0.5` through the public fault factory.
- Use the same requested mesh counts as R2 for FDI and PLI.
- Before projection and mesh rounding, the requested bounds are `O +/- (10+20*b)` on each axis. Thus half-widths are 10, 14, and 20.
- Apply all 3 computational rotations. The expected local box is the enclosure of all 8 projected requested world corners.
- Allow only outward mesh rounding by at most 1 cell per local face. Use the actual local cell spacing for this allowance.
- Require distinct bounds for all 3 requests. A fixed default buffer and a fixed maximum-sized support both fail.

#### Required backend coverage

| Backend or capability | Required fixtures and operations | Gate |
|---|---|---|
| Python FDI and PLI/P1 | R1 evaluation and fitting; R2 fault frame and displacement; R3 relationships; R5 bounds; supported persistence and copies; transform lifecycle | Mandatory for both engines |
| Python P2 | R1 evaluation with the exact nodal field (C-017); R3 reference and row assembly; bounds, persistence, and documented copy operations; transform lifecycle | Mandatory |
| Python FDI-fold, PLI discrete-fold, and P2-fold | R3 point and direction sampling, masks, sign alignment, nonzero targets, and each supported regularisation family; R4 geometry | Mandatory for every existing family; preserve documented unsupported cases |
| Python FDI and PLI constant-norm variants | R1 with its analytic norm target; world queries and gradients; transform lifecycle; inherited affected relationships | Mandatory |
| Common 2D and 3D geometry | Point and direction round trips, corner enclosures, resize stability, cell centres, and auxiliary transformations | Mandatory |
| Existing native FDI or fold integration | R1 and the R3 operations it already supports; world/local conversion at its Python boundary | Run when the existing extension is available; record absence explicitly |
| Optional Surfe integration | R1 world inputs, values, and gradients through its existing adapter | Run when its dependency is available; record absence explicitly |

An unsupported capability must retain its documented rejection. Missing optional dependencies must not skip mandatory Python checks. A future native or RBF engine is not a prerequisite.

- **SC-001**: Rigid point round trips differ by at most `Tp` for 2 and 3 dimensions. Direction round trips have maximum component error at most `1e-12` for unit inputs.
- **SC-002**: R1 value error is at most `1e-8*F`. Its gradient component error is at most `1e-8*max(1,|gradient|)`. The same limits apply after translation and computational or physical rotation. For P2, these limits apply to the exact nodal field of R1 (C-017).
- **SC-003**: R2 has maximum scalar-coordinate error at most `1e-6`, relative gradient-magnitude error at most `1e-4`, and maximum direction error at most 0.01 degree. Sample its 6859 specified probes with both FDI and PLI. Use CG relative tolerance `1e-10`, absolute tolerance `1e-12`, and at most 20000 iterations. Explicit and inferred centres satisfy the same limits.
- **SC-004**: On R2's 1332 off-fault section points, each displacement comparison differs from its analytic-frame result by at most 0.03 world-length units. This is 1 percent of slip 3. Use the same displacement algorithm for both frames. Count fault-plane and invalid-support samples separately.
- **SC-005**: Run R5 for 3 buffers, 3 rotations, and explicit and inferred centres with both FDI and PLI. All requested corners and generated observations are enclosed within `Tp`. Each face differs from its expected bound by at most 1 cell plus `Tp`, outward only. All 3 buffer requests produce distinct bounds. Zero valid observations are lost because of a space mismatch.
- **SC-006**: A save/load cycle preserves R1 values within `1e-12*F` and gradient components within `1e-12*max(1,|gradient|)`. Stored local constraint arrays are unchanged after 2 cycles. Both payloads declare the same supported format version, transform, and constraint coordinate space.
- **SC-007**: All valid reference samples in R3 receive the intended world location within `Tp`. Their local unit directions have component error at most `1e-10`. Equality selections have 0 missing or duplicate valid nodes. Invalid samples have an exact reported exclusion count.
- **SC-008**: R3 affine direction rows have maximum residual at most `1e-10` under each assembler's documented normalisation. Check zero and nonzero targets for FDI, PLI, and P2, including supported quadrature paths. The physical nonzero target remains 1 after rotation.
- **SC-009**: Diagnostics cause 0 changes to stored positions, vectors, weights, coefficients, and cache-validity state. A matching analytic unit direction gives cosine agreement within `1e-10`.
- **SC-010**: Test each documented auxiliary forward/inverse transformation with 1 point and point clouds of 4 and 50 points. Maximum round-trip position error is at most `Tp`. Untransformed columns remain exactly unchanged.
- **SC-011**: Reapplying unchanged support bounds 10 times causes 0 changes to origin, maximum, spacing, and counts. Every grid cell has exactly 1 centre; all centres match their exported physical cells within `Tp`.
- **SC-012**: Invalid rotations and invalid saved metadata produce 1 actionable failure per attempted operation and 0 changes to valid prior state. From the next major release, an ambiguous legacy file fails without an explicit caller choice; on the 1.x line it loads with the identity frame and one warning. With a supplied transform or identity declaration, its values and gradients satisfy SC-006.
- **SC-013**: R4 reference samples and the new automatic direction targets agree with their analytic unit gradients within `1e-10` per component. The centre extent is `(0,1,0)`. The two off-centre sampled extent directions differ by at least 15 degrees. A constant centre direction applied at those probes fails acceptance.
- **SC-014**: 100 percent of affected production paths have regression coverage. Before merge, `pytest` and `ruff check LoopStructural packages` must both pass. Affected package tests and API/import compatibility checks must also pass. Pre-change comparisons diagnose existing failures; they do not waive a failing required gate.
- **SC-015**: 100 percent of changed golden results have analytic evidence and an explicit review record. Exactly 0 unaffected baseline files change. Acceptance tolerances are relaxed in exactly 0 cases to make a correction pass.
- **SC-016**: For populated R1, change the local origin to `O+(3,-2,1)` and rotate 30 degrees about z. Rejection causes 0 state changes. If the policy rejects each later use of the stale state, each use of the populated interpolator raises 1 error. It causes 0 state changes. A rebuild in the new frame preserves world evaluation within SC-002. Acceptance causes 0 stale constraint or coefficient uses and preserves world evaluation within SC-002 after rebuilding. If restoration caches exist (spec 003 adds them), repeat on R2 after populating them. After an accepted change, require 0 stale cache hits and world restored points within `Tp` of an uncached control.
- **SC-017**: Each stage records 1 rollback target and passes 1 rollback compatibility check. Old-format R1 inputs retain their documented load path and SC-006 accuracy. An unsupported new-format input causes 1 actionable failure and 0 mutations of valid state. Test an unknown future format version as a separate failure case.
- **SC-018**: Every mandatory backend-matrix cell has a passing test result and 0 skips. Optional exclusions identify their exact dependency or unsupported capability. A missing optional backend is never reported as having passed.
- **SC-019**: Before planning starts, both affected specs 003 and 004 contain the accepted scoped amendments. There are 0 remaining conflicts about corrected baselines, fingerprints, cache invalidation, or the coordinate-conversion boundary. This gate does not require either solver feature to be implemented.

### Requirement coverage

`U2.3` means user story 2, acceptance scenario 3. The backend matrix applies to each relevant test below.

| Requirement | Acceptance scenario or criterion |
|---|---|
| FR-001 | U1.3; U2.6; SC-002 to SC-004 |
| FR-002 | U1.1; SC-003 |
| FR-003 | U2.1 and U2.2; U3.2; SC-001 |
| FR-004 | U2.1; SC-007 |
| FR-005 | U2.2; SC-007 and SC-008 |
| FR-006 | U2.3; SC-007 |
| FR-007 | U2.4; R3 analytic bounds |
| FR-008 | U2.5; SC-007, SC-008, and SC-018 |
| FR-009 | U2.7; restored world-point enclosure |
| FR-010 | U1.1 to U1.3; SC-003 |
| FR-011 | U1.6 and U1.7; SC-003 and SC-013 |
| FR-012 | U1.4 and U1.5; R5 and SC-005 |
| FR-013 | U4.1; SC-005 |
| FR-014 | U4.2; SC-011 |
| FR-015 | U4.3; SC-011 |
| FR-016 | U3.8; SC-006 and SC-017 |
| FR-017 | U3.2; SC-006 |
| FR-018 | U3.3; backend-matrix copy checks |
| FR-019 | U3.6 and U3.7; SC-012 and SC-017 |
| FR-020 | U3.9; SC-016 |
| FR-021 | U4.6; SC-001 and SC-012 |
| FR-022 | U3.4 and U3.5; world-location comparison within `Tp` |
| FR-023 | U4.4; SC-009 |
| FR-024 | U4.5; SC-003, SC-007, and SC-008 |
| FR-025 | U5.4; SC-010 |
| FR-026 | U5.1; SC-003, SC-004, and SC-008 |
| FR-027 | U5.1 to U5.3; SC-004 |
| FR-028 | SC-012, SC-014, and SC-017 |
| FR-029 | U5.5; SC-015 and SC-019 |
| FR-030 | U5.6; SC-018 |
| FR-031 | SC-014; pre-fix and post-fix regression evidence |
| FR-032 | U3.10; SC-017 |
| FR-033 | U3.11; SC-014 |

## Out of scope

- A new fault-displacement law or a trishear algorithm added to the library.
- A new numerical integrator, fault-discontinuity scheme, or promise of zero deform/restore round-trip error.
- Automatic reconstruction of missing information from an ambiguous saved file.
- General affine scaling or shear support. This spec validates rigid transforms rather than silently applying incorrect gradient mathematics.
- Imposing a planar field on every curved fault or frame.
- A new RBF, native, browser, or GPU solver; performance work beyond checking that correctness changes do not remove existing supported paths.
- A rewrite of the intrusion workflow, QGIS plugin, model graph, or all export formats. The audited coordinate consumers and shared geometry are in scope.
- Changes to other spec documents, branches, worktrees, commits, or release publication during this interview.
- The trace and dip anisotropy settings of a fault after a property change. A rebuild after a property change uses anisotropy weights of 1.0, also when the caller set 0. FR-012 covers only the fault buffer and forced mesh geometry. Issue #394 records this defect (C-005).
- Solver settings that `build` ignores. Top-level `rtol`, `atol` and `maxiter` have no effect, and no warning occurs; only `solver_kwargs` applies them. This is not a coordinate-frame defect. Issue #393 records it (C-006). The tests of this spec give the SC-003 solve settings through `solver_kwargs`.
- The accuracy of the P2 value-only fit of R1. In every frame it misses the SC-002 limits, also in the unrotated control. It depends on the solver and the regularisation, not on a coordinate space. Issue #399 records it (C-017).

## Risks and rollback

Coordinate fixes intentionally change outputs that encoded the defect. This can affect stored golden results and user workarounds that supplied extra observations. The user approved reviewed changes supported by analytic evidence; unaffected baselines remain fixed. Migration notes must explain the corrected world-coordinate behaviour and the effect on such workarounds.

A saved solved field may lack the information needed to recover its physical placement. From the next major release, loading requires a caller-supplied transform or an explicit identity declaration, and the error explains that requirement. On the 1.x line, a legacy file with no choice loads with the identity frame and one warning (C-018).

Stages must retain compatible entry points. A stage can be reverted independently if its acceptance checks fail. Files written in a new format must retain their declared version. A rollback must reject unsupported new files instead of silently treating them as legacy identity files.

Each stage must demonstrate its rollback against the documented reader version and input format. This includes the approved explicit-transform path for ambiguous legacy inputs. It must not claim that an old reader can interpret a new format it does not support.

## Clarifications

### Established from the request and project rules

- **2026-11-03, scope**: The user requested a spec for both the shared coordinate problems and the fault-frame interpolation problems discussed in this chat.
- **2026-11-03, compatibility**: The constitution protects the public model API. ROADMAP.md records that the modelling API uses world coordinates and the interpolation layer owns projection.
- **2026-11-03, architecture**: Existing common geometry, interpolation, and geological modelling packages own the work. No new package is proposed. ROADMAP.md has no `arch` block.
- **2026-11-03, workflow**: The accepted constitution was read. Strict project preflight passed using the source CLI workaround recorded in the project's setup notes. No setup policy was changed.

### Interview decisions answered on 2026-11-03

- **C-001, name**. Question: which name should this spec use? Answer: keep `005-coordinate-frame-correctness`.
- **C-002, extent geometry**. Question: how should supplied centre, normal, and slip directions constrain generated extent geometry? Answer: use those directions at the centre and allow curvature away from it. An explicit planar mode is not required.
- **C-003, legacy persistence**. Question: how should an old interpolator load when its transform cannot be recovered? Answer: require the caller to supply the transform or explicitly declare identity. Do not assume identity with a warning. (C-018 allows the warning-only fallback on the 1.x line until the next major release.)
- **C-004, known-defect baselines**. Question: may analytic evidence justify a change to a golden result containing a coordinate defect? Answer: allow reviewed corrections to affected baselines and preserve unaffected results.

All four interview decisions are resolved. The requirements, acceptance scenarios, edge cases, and success criteria above incorporate the answers.

### Review corrections applied on 2026-11-03

- Replaced the weaker zero-new-failures gate with the constitution's passing test and lint gates.
- Added self-contained R1 to R5 fixtures, a curved-direction check, and two-sided buffer checks with bounded mesh rounding.
- Added saved-format, populated-transform, cache, rollback, and backend-matrix acceptance coverage.
- Made the accepted cascade to specs 003 and 004 a pre-plan prerequisite. This review has not changed those documents.
- Review verdict remains **Needs changes** only because that external prerequisite is still outstanding. No new user-policy clarification is required.

### Decisions after stage S1, answered on 2026-11-04

The S1 build (pull request #392) found two defects that the spec did not decide, and the maintainer decided the plan point OP2.

- **C-005, anisotropy after a property change**. Question: must FR-012 also keep the trace and dip anisotropy settings through a rebuild? Answer: no. This is out of scope for spec 005. Issue #394 records it.
- **C-006, ignored solver settings**. Question: must this spec make `build` apply or reject top-level solver settings? Answer: no. It is not a coordinate-frame defect. Issue #393 records it.
- **C-007, R2 and R5 run time (plan OP2)**. Question: where does the mandatory R2 and R5 matrix run? Answer: in the default `pytest` run, with no marker. In S1, the four new test files take about 5 minutes. The analytic file takes 264 s. The full root suite takes about 14.5 minutes.

### Decisions after stages S2 to S4, answered on 2026-11-05

The maintainer decided these points during the S2 to S4 builds and plans. The plan holds their detail.

- **C-008, constant-norm FDI defect (S2)**. `ConstantNormFDIInterpolator` raised `IndexError` for every frame, because `StructuredGrid.element_size` is one scalar. The defect does not depend on the frame. A separate task fixed it. No requirement changed.
- **C-009, policy for a transform change (S3, plan OP4)**. The interpolator rejects each later use of its stored data in another frame, and changes no state. S3 neither rejects the change itself nor rebuilds automatically. FR-020, U3.9 and SC-016 name this policy.
- **C-010, old reader (S3, plan OP5)**. A reader from before S3 fails on a new-format file with a `TypeError` that names `format_version`, and makes no object. The maintainer accepts this as the failure that SC-017 asks for at the rollback target.
- **C-011, partial legacy transform (S3 end review)**. A legacy file needs the local origin and the local rotation, both, or an identity declaration. A partial input raises an error that names the missing part. FR-019 states this.
- **C-012, P2 fit accuracy (S3, plan OP7)**. The P2 value-only fit of R1 misses SC-002 in every frame. The S3 rebuild test uses the exact nodal control for P2. C-017 decided the scope of SC-002 for P2.
- **C-013, deprecation of `rescale` (S4 plan)**. `GeologicalModel.rescale` keeps its behaviour and warns. FR-033 and U3.11 state this.
- **C-014, version of the legacy-load break**. The first draft of FR-028 asked for a major version with a `BREAKING CHANGE` commit footer. C-015 replaced it.
- **C-015, release of the legacy-load break (2026-11-05, maintainer)**. Question: how does the break of FR-019 ship? Options were a major release with a `BREAKING CHANGE` footer on the squash commit of pull request #397, a hold of stage S3 for a 2.0 line, or a minor release with a note. Answer: a minor release with a note. (C-018 removes the break from 1.x, so no exception is needed there; the rest of this entry applies to the major release.) The rule of `ROADMAP.md` ("1.x stays truly backward compatible") gets an explicit exception. The constitution asks for a "major-version note" for a breaking change. The note is the statement of this decision in `ROADMAP.md`, `COMPAT.md` and the release notes. If the maintainer reads the constitution as a demand for a major release, the maintainer amends it by pull request.
- **C-016, deprecation of `scale` (2026-11-05, maintainer)**. Question: deprecate `GeologicalModel.scale` with `rescale`? Options were: both together, `rescale` only, or `scale` one minor release later. Answer: both together. Nothing in the library calls `scale`, so no library code changes for it. The warning names the bounding-box projection as the replacement. `COMPAT.md` gets one row for each method.
- **C-017, SC-002 for the P2 fit (2026-11-05, maintainer)**. Question: does SC-002 apply to the P2 value-only fit of R1? The fit misses SC-002 in every frame: the gradient error is 3.9e-7 against a limit of 1e-8, also in the unrotated control; the value error is 7.8e-8 against a limit of 7.5e-8 unrotated, and 2e-7 to 7.6e-7 at 30 degrees (the earlier figure of 1.8e-8 was a ratio, not the error). Options were: the exact nodal field only, a frame-invariance bound for the fit, or a fix of the P2 fit in this spec. Answer: the exact nodal field only. SC-002 for P2 checks the coordinate conversion. The fit accuracy is a solver property and goes to a separate issue (#399), as C-005 and C-006 did. The backend table and the out-of-scope list say so.
- **C-018, backward-compatible legacy load (2026-11-05, maintainer)**. The maintainer will release a major version, but not yet. Question: can the change of FR-019 stay backward compatible until then? Answer: yes. On the 1.x line, a legacy file with no choice loads with the identity frame that older releases assumed, and it emits one `FutureWarning` (visible to end users by default, because a wrong placement can be silent). The new keyword arguments (a supplied transform, an identity declaration) work now, and a partial transform is refused. The code carries a note to remove the fallback at the next major release, when the loader raises the error of FR-019. This replaces the sentence "do not assume identity with a warning" of C-003 for the 1.x line only. From the major release, C-003 applies in full. No break ships in 1.x, so C-015 needs no exception to the 1.x rule.

### Decision at the S6 plan, answered on 2026-11-06

- **C-019, displacement evaluator of R2 (2026-11-06, maintainer; plan OP18)**. Question: R2 named an existing test-side displacement evaluator with trishear parameters (apex angle 50 degrees, propagation/slip ratio 1.5, concentration 1, 20 integration steps, tip coordinate 0, propagation sign -1, and 16 quadrature points). The repository has no such evaluator, and the out-of-scope list bans a trishear algorithm in the library. How does the R2 displacement comparison run? Answer: it uses the library displacement paths, with the same algorithm for both frames. The vector comparison uses `FaultSegment.apply_to_points` with 20 steps. The scalar comparison uses `FaultSegment.evaluate_displacement`. The slip is 3 and the profile is the default fault profile. The trishear parameters and the 16 quadrature points are dropped. The R2 fixture text changed to match. SC-004 does not change (1332 smooth points, bound 0.03). No library displacement law is added.

## Next step

S1 is merged (#392). S2 is in pull request #396, and S3 in #397 (stacked on #396). The S4 plan is in `plan.md`. No clarification mark is open. Run the spec-tasks skill for stage S4. Stages S5 and S6 need the spec-planner skill (plan OP1).
