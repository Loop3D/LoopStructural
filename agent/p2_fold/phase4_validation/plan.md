# Phase 4: validation, docs, release

## Goal
Confirm that the benefit seen in the prototype is real after the weight corrections. Then document and release.

## Why re-run
The handoff table used two wrong rules: a useless volume factor, and weight `w / n_samples` (energy `w^2 / n`). Gauss4 had 1/4 of the centroid strength. The balance between fold constraints and regularisation was different. The numbers (4 to 8 deg for gauss4, 6 to 12 deg for linear DFI) may change.

## Validation
Data: Bermagui S2 and S0 (`Data/Processed/TrimDataS4.csv`, `S3`, `S2`, `S0_LA2`). Box: origin (238275, 5965000, -50), maximum (238500, 5965250, 50). The scripts are in the author's Honours project (`debug_s0/scripts/`), not in this repo.

Runs:
1. Linear DFI.
2. `P2Fold` with `centroid`.
3. `P2Fold` with `gauss4`.

Rules for a fair comparison:
- Same data, same fixed seeds, same explicit `buffer` for all engines. The double buffer issue changes mesh size.
- Compare only inside one run. The model chain varies by about 3 deg.
- Use S0 with `limb_wl=None` (auto). Values 2 to 5 give a poor local minimum.
- Metric: angle between the unit gradient and `cross(fold_direction, fold_axis)`. Report median and share above 30 deg.
- Record unknowns, tets, solve time and peak memory.

## Decision gate
- If `gauss4` is better than DFI at equal or fewer unknowns: continue.
- If it is not: stop. Write the new numbers in the handoff. Ask the user before you document `P2Fold` as the better choice. Keep it as an experimental option.

## Docs
Short note, with an example:
- When to use `P2Fold`.
- Mesh size: `nelements` is tets. 1e4 gives 86k unknowns. 1e5 gives 765k unknowns and about 2.5 min.
- `sample_rule` option.
- The field scale is arbitrary (the norm constraint is almost never met). Direction is reliable. A wavelength in field units is not portable between runs.
- Open issues: buffer, F2, F3, rotation-profile fit.

## Release
- Commits use `feat:` for the class and `fix:` for Phase 0 items. Release-please sets the version.
- No Co-Authored-By line.
- Update the handoff file with the corrected table, or mark it as superseded.
- Update `agent/README.md`.

## Done when
- The new table is in the handoff or in the docs.
- The decision gate passed, or the user accepted the result.
- Docs are merged.
