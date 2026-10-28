#!/usr/bin/env bash
# Build the browser package (wasm + JS glue) into pkg/.
# Steps: cargo build, wasm-bindgen, wasm-opt (optional).
set -euo pipefail

# Work from the package root, whatever the current directory is.
cd "$(dirname "$0")/.."

NAME=loop_interpolation_core

# Step 1: build the wasm with the size profile.
cargo build --profile wasm-size --target wasm32-unknown-unknown --features wasm

# Step 2: find the cargo target dir. The user config can set a shared dir.
# Use the first of python3, python that works.
PY=python
command -v python3 >/dev/null 2>&1 && python3 -c "" >/dev/null 2>&1 && PY=python3
TARGET_DIR="$(cargo metadata --format-version 1 --no-deps \
  | "$PY" -c 'import json,sys; print(json.load(sys.stdin)["target_directory"])')"
# Git Bash on Windows: change a Windows path to a POSIX path if cygpath exists.
if command -v cygpath >/dev/null 2>&1; then
  TARGET_DIR="$(cygpath -u "$TARGET_DIR")"
fi
WASM="$TARGET_DIR/wasm32-unknown-unknown/wasm-size/$NAME.wasm"
if [ ! -f "$WASM" ]; then
  echo "error: wasm file not found: $WASM" >&2
  exit 1
fi

# Step 3: the wasm-bindgen CLI version must equal the crate version in Cargo.lock.
CRATE_VER="$(awk '/^name = "wasm-bindgen"$/ {getline; gsub(/[^0-9.]/, "", $0); print; exit}' Cargo.lock)"
CLI_VER="$(wasm-bindgen --version | awk '{print $2}')"
if [ "$CRATE_VER" != "$CLI_VER" ]; then
  echo "error: wasm-bindgen CLI is $CLI_VER but Cargo.lock has $CRATE_VER." >&2
  echo "fix: cargo install wasm-bindgen-cli --version $CRATE_VER --locked" >&2
  exit 1
fi

# Step 4: make the JS glue and the final wasm.
rm -rf pkg
wasm-bindgen --target web --out-dir pkg --out-name "$NAME" "$WASM"

# Step 5: shrink the wasm in place, if wasm-opt is present.
if command -v wasm-opt >/dev/null 2>&1; then
  wasm-opt -Oz -o "pkg/${NAME}_bg.wasm" "pkg/${NAME}_bg.wasm"
else
  echo "wasm-opt not found: skipped"
fi

echo "Built pkg/ in $(pwd)/pkg"
