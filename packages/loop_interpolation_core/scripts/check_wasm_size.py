"""Check the gzip size of the browser wasm (SC-010).

Run with: python -I packages/loop_interpolation_core/scripts/check_wasm_size.py
Exit code 1 if the gzip size is above the limit.
"""
import gzip
import pathlib
import sys

LIMIT = 2_000_000  # bytes, gzip level 9 (SC-010)


def main() -> int:
    # Find pkg/ from the script file, not from the current directory.
    pkg = pathlib.Path(__file__).resolve().parent.parent / "pkg"
    files = sorted(pkg.glob("*_bg.wasm"))
    if not files:
        print(f"error: no *_bg.wasm in {pkg}: run build_wasm.sh first", file=sys.stderr)
        return 2
    status = 0
    for path in files:
        raw = path.read_bytes()
        packed = gzip.compress(raw, compresslevel=9)
        print(f"{path.name}: raw {len(raw)} bytes, gzip {len(packed)} bytes")
        if len(packed) > LIMIT:
            print("Over the limit: see OP2 in plan.md")
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
