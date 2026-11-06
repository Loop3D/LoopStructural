# Changelog

## [0.2.0](https://github.com/Loop3D/LoopStructural/compare/loop-common-v0.1.0...loop-common-v0.2.0) (2026-11-06)


### Features

* add fold frame coordinate validation, P2 support fixes and dev-tools setup ([bb2962b](https://github.com/Loop3D/LoopStructural/commit/bb2962b4d8c02b18d6c0efa78ed9f672a13a2e33))
* trilinear interpolation on PropertyGrid, via StructuredGrid3DGeometry ([#369](https://github.com/Loop3D/LoopStructural/issues/369)) ([8eb8b66](https://github.com/Loop3D/LoopStructural/commit/8eb8b662cf9802d756ddc7c15ea9d49d8fd793c9))


### Bug Fixes

* coordinate and fault-frame correctness, stage S1 (spec 005) ([#392](https://github.com/Loop3D/LoopStructural/issues/392)) ([03a6e58](https://github.com/Loop3D/LoopStructural/commit/03a6e58bdd2f528fc193368df0c3e2617f2d329b))
* make scalar_field / block-model grids line up with their values ([#364](https://github.com/Loop3D/LoopStructural/issues/364)) ([916246e](https://github.com/Loop3D/LoopStructural/commit/916246e11888b921758e0f67e68db58b57b71021))
* make StructuredGrid3DGeometry agree with the structured support it duplicates ([#366](https://github.com/Loop3D/LoopStructural/issues/366)) ([722f50a](https://github.com/Loop3D/LoopStructural/commit/722f50a74be85d7c0fd2d576c1027a7fe4c7d4b5))
* point transformations, fault-frame validation and backend matrix, stage S6 (spec 005) ([#403](https://github.com/Loop3D/LoopStructural/issues/403)) ([6ec0263](https://github.com/Loop3D/LoopStructural/commit/6ec0263876ee9f1ca983eeb206d6d27cb972a4d9))
* resize paths, non-mutating diagnostics and residual report, stage S5 (spec 005) ([#402](https://github.com/Loop3D/LoopStructural/issues/402)) ([e1c10ef](https://github.com/Loop3D/LoopStructural/commit/e1c10ef68bbc320d21e8a4befd9b648be20b11c0))
* **supports:** keep points on tetrahedron faces inside despite rounding ([#404](https://github.com/Loop3D/LoopStructural/issues/404)) ([a5cb5ab](https://github.com/Loop3D/LoopStructural/commit/a5cb5ab25960990d1c36dd6af9556be6d4008455))

## [0.1.0](https://github.com/Loop3D/LoopStructural/compare/loop-common-v0.0.2...loop-common-v0.1.0) (2026-10-20)


### Features

* **interpolation:** P2 phase 1 - P2FoldInterpolator ([#346](https://github.com/Loop3D/LoopStructural/issues/346)) ([a1cae84](https://github.com/Loop3D/LoopStructural/commit/a1cae847a0f9ef31cdf499a93a89d9e2742654d7))


### Bug Fixes

* **interpolation:** P2 phase 0 prerequisites (bugs 3, 4, 5, 6, 11, F1) ([#345](https://github.com/Loop3D/LoopStructural/issues/345)) ([4194616](https://github.com/Loop3D/LoopStructural/commit/41946168b2f60d48508037743d0fc9f98cf7d979))
* **supports:** add to_dict to P2UnstructuredTetMesh so P2 interpolators serialise ([a37895d](https://github.com/Loop3D/LoopStructural/commit/a37895dbf669ce52ba565e8423a32a267bcd0813))

## [0.0.2](https://github.com/Loop3D/LoopStructural/compare/loop-common-v0.0.1...loop-common-v0.0.2) (2026-08-04)


### Bug Fixes

* apply chain rule in StructuredGrid2D element gradient ([08e2e19](https://github.com/Loop3D/LoopStructural/commit/08e2e193fa0de78e1cdaed9ef01c7ce90c4b7ce9))
* apply chain rule in StructuredGrid2D element gradient ([08e2e19](https://github.com/Loop3D/LoopStructural/commit/08e2e193fa0de78e1cdaed9ef01c7ce90c4b7ce9))
* apply chain rule in StructuredGrid2D element gradient ([3d43f5c](https://github.com/Loop3D/LoopStructural/commit/3d43f5ce1b0de3830ec93b814175414820867f4e)), closes [#289](https://github.com/Loop3D/LoopStructural/issues/289)
* remove release-as pin now that loop-common/loop-interpolation have released ([b7f2109](https://github.com/Loop3D/LoopStructural/commit/b7f21093f9e6bdd9001a347c14d1c20ef26f1d57))

## [0.0.1](https://github.com/Loop3D/LoopStructural/compare/loop-common-v0.0.1...loop-common-v0.0.1) (2026-08-03)


### Features

* restructuring loopstructural  ([eaedc35](https://github.com/Loop3D/LoopStructural/commit/eaedc35a6234efa3ee695a91f53c8716031eedaf))
* upgrade bounding box to use origin/maximum and affine transformation ([12f34f1](https://github.com/Loop3D/LoopStructural/commit/12f34f101e83d7159dc489fca227cc6cdb49b799))


### Bug Fixes

* cleaning up PR ([af7385e](https://github.com/Loop3D/LoopStructural/commit/af7385e62211f440966b9275743428e437eeb25d))
* ruff error fixes ([b9ff94d](https://github.com/Loop3D/LoopStructural/commit/b9ff94d5acaa5855cf1037838a20ba6bab93c565))
* updating imports for testing ([27c0836](https://github.com/Loop3D/LoopStructural/commit/27c08361a70cafdc6cb81ebe582a475cdfd8a57e))


### Miscellaneous Chores

* release 1.6.4 ([f06616f](https://github.com/Loop3D/LoopStructural/commit/f06616f8fac0ca3cfc58377524245952f56e686b))
* release 1.6.5 ([246e48d](https://github.com/Loop3D/LoopStructural/commit/246e48d86a99e9d1e96ab9a2d9567374ffcf8622))
