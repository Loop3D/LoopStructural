# loop-model-core

Native (Rust) model evaluation for LoopStructural.

## Fault kernel

The fault kernel restores points through the faults of a LoopStructural
feature. It also evaluates values and gradients. LoopStructural uses it
automatically. Install it with `pip install loopstructural[native]`.

The kernel version is `FAULT_KERNEL_VERSION = 1`. If the version does not
match, LoopStructural logs one warning and uses the Python path.

The classes are `FaultChain` (`restore`, `frame_value`, `feature_gradient`)
and `GridField` (`values`, `gradients`).

The native path is on by default. The switch is
`LoopStructuralConfig.native_faults`. It wins over the environment variable
`LOOPSTRUCTURAL_NATIVE_FAULTS=0`.

The kernel releases the Python GIL and uses a thread pool. Set
`LoopStructuralConfig.native_threads` to the thread count. `None` or `0`
uses all cores.

For the spec dict, see [specs/007-native-fault-kernel/plan.md](../../specs/007-native-fault-kernel/plan.md) (D2b).
