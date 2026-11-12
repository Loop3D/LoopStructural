# Labelled tetrahedral meshes with MMG3D (experimental)

`make_mmg_volume.py` makes, for a benchmark case, one conforming tetrahedral
mesh of the unit cube in which every tetrahedron has a stratigraphic unit and a
fault block, plus one OBJ per fault and one OBJ per horizon. The OBJ surfaces
are faces of the volume mesh, so they share vertices with the volume and with
each other.

The purpose is to learn what a fit-for-purpose geology mesher (planned in
Rust, following the logic of MMG) must do differently. The
[issue list](#mmg-issues-and-what-the-rust-mesher-should-do) is the main
result.

```sh
# from this folder, with the LoopStructural .venv; MMG 5.8 is found in
# ../fault_mmg_experiment/.mmg-runtime, or set MMG3D=/path/to/mmg3d_O3
python make_mmg_volume.py --case unconformity_faulted --reference <path>/benchmark
```

## Method

MMG3D `-ls` inserts the zero set of **one** scalar field per call. The mesher
therefore makes a sequence of calls. Each call uses an `LSReferences` table:
the materials (tetrahedron references) where the surface exists are split into
two **new** references; all other materials are `nosplit`.

A Python material table (`materials.json`) keeps, for every MMG reference:
the side of each inserted fault, the raw `gx` side, the series region, and the
range of units that the material can still contain. After the last cut, each
reference has one fault block and one unit.

Cut order, youngest structure first, so that each structure is inserted only
where it exists:

1. For each series (youngest first): its faults, in `fault_build_order`
   (younger and host faults first); then its unconformity surface. A series
   without an unconformity gets all material that is still unassigned.
2. The horizons of each series, one call per isovalue and per fault block.

Each field is the **smooth extension of one fault block** (`fields.forced_sides`),
evaluated with LoopStructural at the current MMG vertices. The field is never
sampled across a fault, so there is no smearing. An abutting fault is cut only
in materials on its host side. A fault that an unconformity truncates is cut
only in the material below that unconformity.

A fault with a tip is cut as its full extension `gx = 0`, then a tip surface
`t = 0` in the footwall only. The geological footwall block is `gx < 0, t > 0`.
The extension and the cap stay in the mesh as interfaces between materials
whose geology is the same; they are not exported.

Surfaces come from tetrahedron faces between materials: a fault face separates
two different sides of a fault *and* two raw `gx` sides; a horizon face
separates two materials of the same fault block with units on each side of the
isovalue.

After each cut MMG remeshes the whole mesh (`-hsiz 0.05 -hausd 0.002
-hgrad 1.4`); a last `-optim` pass improves quality. `--no-adapt` gives pure
cuts (`-noinsert -nomove -noswap`): exact, but with slivers.

## Outputs per case

| File | Content |
| --- | --- |
| `volume.vtu` | Tetrahedra. Cell data: `unit`, `series`, `block:<series>`, `side:<fault>`, `material` (MMG reference), `quality`, `min_dihedral_deg` |
| `volume.mesh` | Same mesh, Medit format, reference = material |
| `materials.json` | Unit table, material table (reference → sides, region, unit), and the lineage of every reference |
| `surfaces/fault_<f>.obj` | One surface per fault (normal toward the hanging wall) |
| `surfaces/horizon_<series>_<k>.obj` | One surface per horizon (normal toward the younger unit); same names as the benchmark |
| `surfaces/*.vtp` | Same surfaces with cell data `block` |
| `report.json` | Volume and surface QA, and every MMG call (command, time, size) |

`volume.vtu` and `volume.mesh` (5-10 MB each) are not committed; run the
command again to make them.

## Test set and results

The six cases below test all the hard meshing features (see the
"tests" column). Starting grid 20³ Kuhn cube, default MMG options. The runs
for the last five cases were made at the same time (five processes), so the
times are only approximate.

| Case | Tests | MMG calls | Tets | Surfaces | Volume checks | Unexplained interfaces | Min quality | Min dihedral | Max area error vs reference | Max vertex distance (fault / horizon) | Mesh time |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `unconformity_faulted` | unconformity, truncated fault, displaced fault, 2 series | 19 | 121,364 | 9/9 | pass | 0 | 0.024 | 1.3° | 3.8e-5 | 1.6e-4 / 4.0e-5 | 217 s |
| `squeezed_wedge` | abutment at 20°, wedge thinning to zero | 16 | 112,159 | 7/7 | pass | 0 | 0.175 | 8.0° | 1.1e-7 | 6.0e-8 / 3.7e-5 | 206 s |
| `narrow_horst` | block 0.022-0.058 wide (thinner than `hsiz`) | 15 | 112,866 | 6/6 | pass | 0 | 0.295 | 14.2° | 3.2e-8 | 1.1e-7 / 2.0e-5 | 195 s |
| `listric_fault` | curved fault, rollover | 10 | 113,818 | 5/5 | pass | 0 | 0.180 | 10.6° | 2.0e-3 | 1.7e-3 / **2.2e-2** | 139 s |
| `folded_faulted` | folded horizons, reverse fault, displaced fault, variable throw | 18 | 164,672 | 7/7 | pass | 0 | 0.015 | 1.5° | **9.0e-3** | 2.0e-4 / 4.9e-3 | 705 s |
| `fault_tip` | tip loop inside the volume | 10 | 136,014 | 5/5 | pass | **3,079** (see below) | **1.6e-5** | **0.0007°** | 3.5e-3 | 2.9e-8 / 1.9e-3 | 227 s |

"Volume checks" are: total volume 1, no unlabelled or inverted tets, no face
with more than two tets, no boundary face away from the cube. For
`unconformity_faulted`, every inner boundary edge of each surface is also an
edge of another surface (0 orphan edges); `f_pre` ends on `f_post` and on the
unconformity.

Notes on the cases:

* **`fault_tip`** (first run, before the cap fix). The 3,079 "unexplained"
  faces are the tip cap (3,070 faces, area 1.36; at that time the QA did not
  know about caps) plus **9 faces of almost zero area where the unit changes
  across the cap**. These 9 faces are a real defect: the horizon is cut first
  in the `+` block, MMG remeshing then moves those vertices along the cap within
  `hausd`, and the cut in the `-` block later finds a slightly different zero
  crossing. The two pieces of the horizon do not meet at the same vertex, and
  slivers form. All 28 tets with a dihedral angle below 1° touch these faces or
  the fault. The fault area is 0.35% smaller than the reference: the tip loop
  is a polygon of mesh edges.
  *Fix in the code, not tested yet* (the rerun was stopped because the machine
  was low on memory): when a surface is cut in a second block, vertices that
  are already on that surface and on a tip cap get the value 0 exactly, so MMG
  reuses them. The QA now counts cap faces separately (`tip_cap_faces`).
* **`listric_fault`**: the tip cut was skipped (the tip loop is outside the
  cube), so it does not test tips. `horizon_strat_1` has vertices up to 2.2e-2
  from the true surface (10× `hausd`); the cause was not found (the check was
  stopped because of low memory). Possible causes: vertices moved by
  remeshing near the curved fault, or a weak QA estimate where the gradient is
  small.
* **`folded_faulted`**: both tip cuts were skipped (tips outside the cube).
  Folded horizons are 0.4-0.9% smaller than the reference and up to 4.9e-3 from
  the true surface (2.5× `hausd`): MMG's surface model flattens the folds. Part
  of the area difference can come from the marching-cubes reference. One
  horizon call took 147 s, others 5-60 s.
* **`narrow_horst`**: the narrow block was kept with no special size control:
  the level-set cut inserts the interface whatever the target size is.

## MMG issues and what the Rust mesher should do

From these runs, the earlier `../fault_mmg_experiment`, and the MMG 5.8.0
source (`src/mmg3d/mmg3d2.c`).

| # | MMG behaviour | Effect seen | Rust mesher |
| --- | --- | --- | --- |
| 1 | `-ls` takes one scalar per vertex and one isovalue per call. | 10-19 calls per case; about one per horizon per fault block. | Cut many surfaces in one pass: a list of (surface id, field, region) jobs, applied in one remesh loop. |
| 2 | One value per vertex: a faulted (discontinuous) field cannot be given. | Each fault block needs its own call, with `nosplit` on all other blocks. | Store field values per (vertex, block), or evaluate the field of the tet's block on demand. A vertex on a fault has one value per side. |
| 3 | The field is sampled once on the current mesh, then MMG remeshes with its own surface model (no callback to the true field). | Vertices drift from the true surface: up to 2.5× `hausd` on folds, 10× on one listric horizon; folded areas 0.4-0.9% low. | Keep the implicit function in the loop: project new and moved interface vertices onto the true zero set (Newton step on the block field). |
| 4 | Each `-ls` call remeshes the **whole** mesh; there is no local mode. | 5-150 s per call; most of the time is spent on material that the cut does not touch. | Remesh only a band around the new interface; freeze the rest. |
| 5 | An interface is always a closed boundary between two materials; it cannot end inside the volume. | Fault tips need a non-geological extension and cap (area 0.65 + 1.36 in `fault_tip`), and these make slivers. | Support open internal surfaces (boundary-of-surface edges inside the volume), as tags on faces, not only as material changes. MMG's `-opnbdy` (keep triangles between two domains of the same reference) is the closest option; not tested. |
| 6 | Pieces of one surface from separate calls do not share their intersection points; there is no "this vertex is already on surface S". | `fault_tip`: horizon pieces from two blocks meet at the cap with almost-zero faces and 0.0007° slivers. | Keep a vertex → surfaces incidence table. A vertex already on S has value 0 for S by definition. |
| 7 | Material references carry no meaning; each split needs fresh references; the CLI gives every new interface triangle reference 10. | All meaning is kept in a Python table (`materials.json`); old interface tags must be renamed after each call. | Make the region label a structured value (fault-side bitset + unit id), not an integer that the caller must map. |
| 8 | `snpval_ls`: values with \|v\| < 1e-6 (absolute, not scaled to the field) snap to 0; the cut point on an edge is clamped to [1e-6, 1-1e-6] of the edge. | Behaviour depends on the units of the field; a vertex can never be exactly reused by a cut unless its value is exactly 0. | Use tolerances relative to the edge length and the field gradient; reuse an existing vertex when the cut point is within a tolerance. |
| 9 | `snpval_ls`: around a tet with quality < 1e-6, if any vertex value is < 1e-3, all four values are set to -1e-3. | Silent sign changes: in `fault_mmg_experiment` a horizon moved and MMG returned success. | Never change the field values; repair bad tets before or after the cut, and report every change. |
| 10 | `cuttet_ls`: required edges and triangles are split with only a warning ("Required entity ignored"). | You cannot freeze an interface and cut a neighbour block separately, then join the meshes. | Treat protected edges as hard constraints: cut the other way or refuse, and report. |
| 11 | `setref_ls`: a tet with all four values 0 is only an `assert` (not checked in release builds). | A surface that lies on existing faces (for example two equal horizons) gets an undefined label. | Handle "surface already present" as a valid case: relabel, do not cut. |
| 12 | Quality near junctions is not repaired: vertices on a line where several interfaces meet cannot move. | Min dihedral 1.3° (`unconformity_faulted`) and 1.5° (`folded_faulted`) at junctions; `-optim` does not fix them. | Allow junction vertices to slide along the junction curve (the intersection of two implicit surfaces), with projection. |
| 13 | CLI only (the Windows conda build has no Python bindings): ASCII Medit files for every call. | IO adds time and loses all non-mesh data between calls. | Library API with the mesh in memory. |

## Next steps

- When memory is free: rerun `fault_tip` to test the cap fix, and find the
  cause of the 2.2e-2 vertex distance on `listric_fault` `horizon_strat_1`.
- Run one case at a time (each run uses LoopStructural plus MMG; five at once
  used up the memory).
