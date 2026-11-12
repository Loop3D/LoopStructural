"""Labelled tetrahedral meshes of the benchmark cases with MMG3D (experimental).

Each case is cut by a sequence of MMG3D level-set calls, one smooth
block-local field per call, with an LSReferences material table so that only
the materials where a surface exists are split. Writes, per case:

    <out>/<case>/volume.vtu        tetrahedra, cell data unit, series, block:<series>,
                                   side:<fault>, material, quality, min_dihedral_deg
    <out>/<case>/volume.mesh       the same mesh in Medit format (reference = material)
    <out>/<case>/materials.json    material table: MMG reference -> fault sides, unit
    <out>/<case>/surfaces/fault_<f>.obj            one file per fault
    <out>/<case>/surfaces/horizon_<series>_<k>.obj one file per horizon
    <out>/<case>/report.json       QA and every MMG call

    python make_mmg_volume.py --case unconformity_faulted
    python make_mmg_volume.py --case all --cells 24 --hsiz 0.04

MMG 5.8 (mmgsuite) is needed: set MMG3D or pass --mmg, or keep the
runtime in ../fault_mmg_experiment/.mmg-runtime. See MMG_VOLUME.md.
"""

import argparse
import json
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from fault_benchmarks.builder import build_case, quiet_loopstructural  # noqa: E402
from fault_benchmarks.cases import CASES, get_case  # noqa: E402
from fault_benchmarks.mmg_volume import (  # noqa: E402
    MmgOptions,
    VolumeMesher,
    VolumeOptions,
    write_case,
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", action="append", required=True, help="case name or 'all'")
    ap.add_argument("--out", default="outputs/mmg_volume")
    ap.add_argument("--reference", default="benchmark",
                    help="benchmark folder with <case>/surfaces/*.obj for area checks")
    ap.add_argument("--cells", type=int, default=20, help="starting grid cells per axis")
    ap.add_argument("--hsiz", type=float, default=0.05, help="MMG target edge size")
    ap.add_argument("--hausd", type=float, default=0.002, help="MMG Hausdorff tolerance")
    ap.add_argument("--hgrad", type=float, default=1.4)
    ap.add_argument("--hmin", type=float, default=None)
    ap.add_argument("--no-adapt", action="store_true",
                    help="pure cuts (-noinsert -nomove -noswap): exact, but slivers")
    ap.add_argument("--no-optim", action="store_true", help="skip the final -optim pass")
    ap.add_argument("--keep-files", action="store_true", help="keep every MMG input/output")
    ap.add_argument("--mmg", default=None, help="path to mmg3d_O3")
    args = ap.parse_args()
    names = list(CASES) if "all" in args.case else args.case
    quiet_loopstructural()
    opts = VolumeOptions(
        cells=args.cells, final_optim=not args.no_optim,
        mmg=MmgOptions(executable=args.mmg, adapt=not args.no_adapt, hsiz=args.hsiz,
                       hausd=args.hausd, hgrad=args.hgrad, hmin=args.hmin,
                       keep_files=args.keep_files))
    out_root = Path(args.out)
    index = {}
    ok = True
    for name in names:
        t0 = time.perf_counter()
        print(f"[{name}] building LoopStructural model ...", flush=True)
        try:
            built = build_case(get_case(name))
            t1 = time.perf_counter()
            mesher = VolumeMesher(built, out_root / name / "mmg_runs", opts)
            mesher.run()
            t2 = time.perf_counter()
            ref = Path(args.reference) / name / "surfaces"
            rep = write_case(mesher, out_root / name, ref if ref.exists() else None)
            rep["timings_seconds"] = {"build": t1 - t0, "mesh": t2 - t1,
                                      "write_qa": time.perf_counter() - t2}
            (out_root / name / "report.json").write_text(json.dumps(rep, indent=1, default=str))
            v = rep["volume"]
            good = (v["unlabelled_tets"] == 0 and v["nonpositive_tets"] == 0
                    and v["faces_with_more_than_two_tets"] == 0
                    and v["boundary_faces_not_on_cube"] == 0
                    and v["interface_faces"]["unexplained"] == 0
                    and abs(v["volume"] - 1.0) < 1e-9 and not rep["missing_surfaces"])
            index[name] = {"status": "ok" if good else "check", "mmg_calls": rep["mmg_calls"],
                           "tets": v["n_tets"], "min_quality": v["quality"]["min"],
                           "surfaces": len(rep["surfaces"]),
                           "missing": rep["missing_surfaces"],
                           "seconds": round(time.perf_counter() - t0, 1)}
            print(f"[{name}] {index[name]}", flush=True)
            ok &= good
        except Exception as exc:  # noqa: BLE001 - record and continue with other cases
            traceback.print_exc()
            index[name] = {"status": "failed", "error": f"{type(exc).__name__}: {exc}",
                           "seconds": round(time.perf_counter() - t0, 1)}
            ok = False
    out_root.mkdir(parents=True, exist_ok=True)
    prev = out_root / "index.json"
    merged = json.loads(prev.read_text()) if prev.exists() else {}
    merged.update(index)
    prev.write_text(json.dumps(merged, indent=1))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
