# loop-interpolation-core

## Purpose

This package is the native core of the low-memory interpolator. It is a
Python extension module, built with maturin. It is named
`loop_interpolation_core` on import. The package solves the finite
difference interpolation problem without building the sparse constraint
rows in Python.

The package is optional. The `loop-interpolation` package works without it.

## Specification

The requirements and the design are in
[specs/001-low-memory-interpolator/spec.md](../../specs/001-low-memory-interpolator/spec.md).

## Source

Source: loop-rs crates/interpolation-fd, commit bff51342da8ebff7ca3c87e07c1d4fbd31279c8b, branch lachlan/drillhole
