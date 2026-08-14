# Changelog

## [1.7.1](https://github.com/Loop3D/LoopStructural/compare/v1.7.0...v1.7.1) (2026-08-14)


### Bug Fixes

* forcing release ([10e89dc](https://github.com/Loop3D/LoopStructural/commit/10e89dc8a1137e065da97e01e4cbb2e59edd4642))

## [1.7.0](https://github.com/Loop3D/LoopStructural/compare/v1.6.28...v1.7.0) (2026-08-03)


### Features

* add api registry ([847634f](https://github.com/Loop3D/LoopStructural/commit/847634ff42692d348fdfcda2d7885506e7b28f6e))
* add recipe serialization methods to GeologicalModel and corresponding tests ([6a56793](https://github.com/Loop3D/LoopStructural/commit/6a5679320549f2ffbaaf2ed46c350ba296f544b6))
* enhance API stability by registering external symbols and adding stability checks ([c13d11b](https://github.com/Loop3D/LoopStructural/commit/c13d11b3b6287a1675103e5206a7ad32a40c21cc))
* implement fault cycle detection in geological features ([228a0aa](https://github.com/Loop3D/LoopStructural/commit/228a0aa53a1ffeba11f5dfe77594f0d209b466e4))
* implement JSON serialization for GeologicalModel with recipe methods and tests ([ca9cab6](https://github.com/Loop3D/LoopStructural/commit/ca9cab6f546eeb1596c0164efae91ac8ed964260))
* restructuring loopstructural  ([eaedc35](https://github.com/Loop3D/LoopStructural/commit/eaedc35a6234efa3ee695a91f53c8716031eedaf))
* upgrade bounding box to use origin/maximum and affine transformation ([12f34f1](https://github.com/Loop3D/LoopStructural/commit/12f34f101e83d7159dc489fca227cc6cdb49b799))


### Bug Fixes

* add datatypes deprecation ([719763a](https://github.com/Loop3D/LoopStructural/commit/719763a0d334685da2b32a4a427216dbf4c9a32f))
* build P1Unstructured2d from a bounding box, not just explicit arrays ([758cc8e](https://github.com/Loop3D/LoopStructural/commit/758cc8ef3474483bb5307a0a233ab78a9aae5556))
* build P2Unstructured2d from a bounding box and add proper value/gradient evaluation ([2197255](https://github.com/Loop3D/LoopStructural/commit/2197255abb534480ce0d6a4f4c1b5ae3229a85f2))
* build P2UnstructuredTetMesh from a bounding box, not just explicit arrays ([0c8133f](https://github.com/Loop3D/LoopStructural/commit/0c8133f46f7daa0a4e05bffa093e7c3dd960b75e))
* change TypeError to ValueError for fold_frame validation ([d3aac52](https://github.com/Loop3D/LoopStructural/commit/d3aac52de15ff8879d49439c9574b7cfef5869c1))
* cleaning up PR ([af7385e](https://github.com/Loop3D/LoopStructural/commit/af7385e62211f440966b9275743428e437eeb25d))
* **codestyle:** enforcing/checking coding style ([01bcef5](https://github.com/Loop3D/LoopStructural/commit/01bcef59e0df1b19ee1da911d6394949536ef777))
* correct barycentric coordinate columns and chunking bug in 2D unstructured meshes ([52e0740](https://github.com/Loop3D/LoopStructural/commit/52e0740a30946d0e09bad48012f07de872364e09))
* correct FaultTopology adjacency bookkeeping and dict round-trip ([fc8bb8b](https://github.com/Loop3D/LoopStructural/commit/fc8bb8b91659856d0cae1111417db90a02facd50))
* correct get_quadrature_points array shape in P2Unstructured2d ([8ca6b7b](https://github.com/Loop3D/LoopStructural/commit/8ca6b7b9b9b14a3370915eae09bb34c397770a72))
* correct P2 interpolator gradient constraint indexing and NaN check ([f7047ac](https://github.com/Loop3D/LoopStructural/commit/f7047ac7f50fedb6db7b4a34edff102be4c15f78))
* don't drop single-point value constraints in P1/P2 interpolators ([38c4f8e](https://github.com/Loop3D/LoopStructural/commit/38c4f8e6b82e1c8b478a66648f5b05206f22d683))
* ensure isovalues are correctly set for zero input in LoopIsosurfacer.fit() ([2b703ac](https://github.com/Loop3D/LoopStructural/commit/2b703acebbf5123cca06e87e574ca01864af2e57))
* guard against divide-by-zero in constant norm constraint ([8314164](https://github.com/Loop3D/LoopStructural/commit/8314164b745df97102722614fc9df4f450c9b8e1))
* harden exception handling and unsafe model deserialization ([e5bece7](https://github.com/Loop3D/LoopStructural/commit/e5bece7f1a8068541fb0fcef01d8f6a96403dcea))
* implement P2Unstructured2d.evaluate_shape_d2 using self.hessian ([c9c8a8a](https://github.com/Loop3D/LoopStructural/commit/c9c8a8a652ea19d9d26b72c153281f34eb2da129))
* process disjoint chunks in tetra get_element_for_location ([16c803c](https://github.com/Loop3D/LoopStructural/commit/16c803cf2eb6bedd4d45137e5124876a9ff76b8f))
* raise clear error for unimplemented GOCAD VSet point export ([2d17fef](https://github.com/Loop3D/LoopStructural/commit/2d17fef62e7a8c03ddda6f678db981e00f4ac2a4))
* raise instead of silently no-op on unsupported OMF structured grid export ([3d9c1b7](https://github.com/Loop3D/LoopStructural/commit/3d9c1b7572d5dbf5675b36af609930241c6654b0))
* raise ValueError instead of bare except + BaseException in data setter ([15e1f15](https://github.com/Loop3D/LoopStructural/commit/15e1f15220ebb711bbe148e6f7e5af5ef0fabd2f))
* reject nan/inf constraints and stop leaking solver_kwargs state ([d16cc38](https://github.com/Loop3D/LoopStructural/commit/d16cc38d3244b1ac12f2488188bd5b1cb2ab6d9b))
* resolve basement handling in StratigraphicColumn.clear() and improve observer method tracking ([50b207d](https://github.com/Loop3D/LoopStructural/commit/50b207d94af51996107bae043f82ca1acdbfaf97))
* ruff error fixes ([b9ff94d](https://github.com/Loop3D/LoopStructural/commit/b9ff94d5acaa5855cf1037838a20ba6bab93c565))
* stop sharing BaseFeature.regions list across instances ([4f1608f](https://github.com/Loop3D/LoopStructural/commit/4f1608fd2ce7ccfb66a6c19c007cd8c2cac0f175))
* sync __all__ with LoopStructural's actual public API ([b7ce574](https://github.com/Loop3D/LoopStructural/commit/b7ce574dc34adef8e96c27f6b95212f767e68f9e))
* update parameter name from 'group' to 'groupname' in add_surface_to_geoh5 function ([a0daa26](https://github.com/Loop3D/LoopStructural/commit/a0daa2678da9f42e4202e56e2fe688206b7a8a65))
* update unconformity handling in GeologicalModel and add performance example for evaluation ([0226b5c](https://github.com/Loop3D/LoopStructural/commit/0226b5cf408df11d002ab51de27e1764924354d5))
* updating imports for testing ([27c0836](https://github.com/Loop3D/LoopStructural/commit/27c08361a70cafdc6cb81ebe582a475cdfd8a57e))
* use self.dimensions instead of hardcoded 3D column slicing in P1/P2 constraints ([c9448bc](https://github.com/Loop3D/LoopStructural/commit/c9448bc8da01f3952c21f00d417535aea67d1eda))


### Performance Improvements

* vectorize identified per-row Python loops in intrusion and svariogram helpers ([33099f2](https://github.com/Loop3D/LoopStructural/commit/33099f29d1b50323e71761a165643f20b1e97752))


### Documentation

* fill placeholder docstrings, fix stale param names, wire utils into API ([2c47646](https://github.com/Loop3D/LoopStructural/commit/2c47646a721917dbe870d5c6ca3fb99336cf320d))
* update compatibility shims and roadmap for extracted package integration ([3939bb6](https://github.com/Loop3D/LoopStructural/commit/3939bb6be0b7d6036f433b60cd04c22b2e6fb2e4))

## [1.6.28](https://github.com/Loop3D/LoopStructural/compare/v1.6.27...v1.6.28) (2026-06-03)


### Bug Fixes

* add default for z ([12bdcd2](https://github.com/Loop3D/LoopStructural/commit/12bdcd2b1f9dcf8a5c3e19aa62dd68a11bd03f7f))
* adding add_group and add points from dataframe for geoh5 + some tests ([9d29d31](https://github.com/Loop3D/LoopStructural/commit/9d29d31c28f619438fe41af806d9332181c0eb87))
* adding gocad export for structured grid support ([228ca78](https://github.com/Loop3D/LoopStructural/commit/228ca78889c9a263c83da54837342adc69a5b083))
* adding post_init to ValuePoints/VectorPoints to validate as numpy arrays ([df938d1](https://github.com/Loop3D/LoopStructural/commit/df938d18c9ed5524628e5bd7925536e7b9687f12))
* vtk export had maximum &gt; grid maximum ([5f37322](https://github.com/Loop3D/LoopStructural/commit/5f37322e0b9c27cb68b0628410e586dd8d40a200))

## [1.6.27](https://github.com/Loop3D/LoopStructural/compare/v1.6.26...v1.6.27) (2026-04-14)


### Bug Fixes

* add caching of reused values ([54fba19](https://github.com/Loop3D/LoopStructural/commit/54fba19df600210e6c14fb895d20b7f00451f663))
* numpy speed improvements ([5e29fa1](https://github.com/Loop3D/LoopStructural/commit/5e29fa15b13a0864baa83ef9c95feef99d4cae01))
* remove redundant call to mask nan values in discrete interpolator ([ce082f6](https://github.com/Loop3D/LoopStructural/commit/ce082f6e49bda22f8e02346e4c9cea9a81ba2bc0))

## [1.6.26](https://github.com/Loop3D/LoopStructural/compare/v1.6.25...v1.6.26) (2026-02-11)


### Bug Fixes

* add threshold variable and remove duplicate code ([4d60155](https://github.com/Loop3D/LoopStructural/commit/4d60155ebe079bfee4f23b1cb6a607778728de0f))
* calculate fault normal from plane of points as well as trace ([71c6f17](https://github.com/Loop3D/LoopStructural/commit/71c6f17169326daf847dcabfcf4a995a99d8de9c))
* ensure all data are copies of original datastructure ([f409ad8](https://github.com/Loop3D/LoopStructural/commit/f409ad88b130cc770bd98ec9f459909edf037da8))
* incorrect boolean ([5d9e04b](https://github.com/Loop3D/LoopStructural/commit/5d9e04b3532cd6a05dd3d1bf18ddafd4410626a0))

## [1.6.25](https://github.com/Loop3D/LoopStructural/compare/v1.6.24...v1.6.25) (2026-01-31)


### Bug Fixes

* prevent int overflow if bb min/max is defined as int ([65d2819](https://github.com/Loop3D/LoopStructural/commit/65d2819fcf6fe0b0e0440f487d3cf75e629341b7))
* speed up model eval by collecting all conformable units ([9b4a719](https://github.com/Loop3D/LoopStructural/commit/9b4a7192fb75e5f2e002e8079891d0a7dd5487bc))

## [1.6.24](https://github.com/Loop3D/LoopStructural/compare/v1.6.23...v1.6.24) (2026-01-15)


### Bug Fixes

* change default regularisation to 1.0 as 0.1 was causing overfitting ([fd3fd20](https://github.com/Loop3D/LoopStructural/commit/fd3fd201e743d421708d0999b367a77cb290696c))

## [1.6.23](https://github.com/Loop3D/LoopStructural/compare/v1.6.22...v1.6.23) (2025-11-05)


### Bug Fixes

* add dipdirection2vector helper ([717d6c9](https://github.com/Loop3D/LoopStructural/commit/717d6c940665791c44aeac0e8e1737b3b8cd48e0))
* add function to clean nan vertices and post_init to validate data ([f69f6d8](https://github.com/Loop3D/LoopStructural/commit/f69f6d8132ec5f14d62785183f1f5c447414c937))
* add horizontal model with change in thickness example ([6a5ee85](https://github.com/Loop3D/LoopStructural/commit/6a5ee857f90f02dd178e83d7b3d685538bc7b58c))
* add inequalities to interpolator builder ([2ee38ac](https://github.com/Loop3D/LoopStructural/commit/2ee38acde97bab5a966c1995d7f70fc98a40651c))
* add methods to fold rotation angle feature ([efb46a2](https://github.com/Loop3D/LoopStructural/commit/efb46a23e867eecaedb930a768f8c616ab2313d1))
* add scaling/transformation matrix to bounding box ([984ab1c](https://github.com/Loop3D/LoopStructural/commit/984ab1c4551e92243e50c0ab7a3ce14eb580bfa0))
* adding analytical fold builder ([6d238aa](https://github.com/Loop3D/LoopStructural/commit/6d238aa01bb3ed08eaf571dad4e40ce0ac05b86f))
* adding groupname to save option to geoh5 ([d8a8c63](https://github.com/Loop3D/LoopStructural/commit/d8a8c63cd9ce0269c8aad55df9dc492b17895e8d))
* adding min/max angle to trig profile ([c097372](https://github.com/Loop3D/LoopStructural/commit/c09737245c31a0b18ee6b35d2271edbc672bf939))
* check dimensions of bounding box constructor ([39b694b](https://github.com/Loop3D/LoopStructural/commit/39b694b281502de8a7997ccf010dcb5d9ee4a17f))
* evaluate gradient for structural frame calls interpolator, now working ([596e59f](https://github.com/Loop3D/LoopStructural/commit/596e59f1d28d68a54be321253c13cd8f02c94175))
* remove fold rotation inversion when axis and cross product are not in the same direction. ([3d4d90f](https://github.com/Loop3D/LoopStructural/commit/3d4d90f90ed0dd54db74cf354cf4710fd9f99321))
* remove requirement for featurename in dataframe used to directly construct feature ([09940ba](https://github.com/Loop3D/LoopStructural/commit/09940ba6e91c9293aa98261f8d9ba195bd2d50de))
* run interpolator before trying to create vtk object ([b230152](https://github.com/Loop3D/LoopStructural/commit/b230152e9ec9137c1f1279ffddc6d30848b8f3b2))


### Documentation

* copilot review/update of docstrings to numpy format ([#273](https://github.com/Loop3D/LoopStructural/issues/273)) ([89a3c5f](https://github.com/Loop3D/LoopStructural/commit/89a3c5f8ae01b7b80c2e9bffd33e7f6339cc81e1))

## [1.6.22](https://github.com/Loop3D/LoopStructural/compare/v1.6.21...v1.6.22) (2025-09-08)


### Bug Fixes

* add export of inequalities ([05a160e](https://github.com/Loop3D/LoopStructural/commit/05a160ee142cb8f49b715104c55a82023f569927))
* dictionary serialisation was reversing order ([8bb63bb](https://github.com/Loop3D/LoopStructural/commit/8bb63bbeb49ed83745f1307210db056a53492931))
* disable inequality data ([c9d16e9](https://github.com/Loop3D/LoopStructural/commit/c9d16e9717b470c29a11c93dd43664de0c99e5c4))

## [1.6.21](https://github.com/Loop3D/LoopStructural/compare/v1.6.20...v1.6.21) (2025-09-02)


### Bug Fixes

* add model reference when converting from feature to frame ([ef07373](https://github.com/Loop3D/LoopStructural/commit/ef07373a816e60f4e931f4f32975ded13eb528f6))
* change fault axis to info not warning ([68abea7](https://github.com/Loop3D/LoopStructural/commit/68abea77ceaad03db2255da0a7f6bcd65ee9df4c))
* clear column actually removes elements. Also load from dict adds in reversed order to maintain correct order ([35ff6a6](https://github.com/Loop3D/LoopStructural/commit/35ff6a6427b32cb0ee45c3cedae998477b49a491))

## [1.6.20](https://github.com/Loop3D/LoopStructural/compare/v1.6.19...v1.6.20) (2025-08-25)


### Bug Fixes

* add convert from feature to structural frame ([5b34346](https://github.com/Loop3D/LoopStructural/commit/5b34346e0b930f07044244916ded99d55cb4b3e3))
* add type to P1 and P2 ([71ce492](https://github.com/Loop3D/LoopStructural/commit/71ce4922f917e59913fe948b63f632beda8cbaaf))
* add wrapper to convert between a feature and structural frame ([676c0ef](https://github.com/Loop3D/LoopStructural/commit/676c0ef4a01b777d5c934e55a3202a4f29ea266e))
* allow model to set any type of feature using basefeature subclass check ([97122a5](https://github.com/Loop3D/LoopStructural/commit/97122a5bdf0b2e49b3e56fd7a2b22c40bd94314e))
* allow no isovalue for surfaces. Take middle as value ([178d78b](https://github.com/Loop3D/LoopStructural/commit/178d78bbe775e77e0294875e3f56c7d7ce41d3fa))
* change kwargs to build args ([cb0af22](https://github.com/Loop3D/LoopStructural/commit/cb0af22ef96ddf21e06ca6ae0250a45de52de220))
* linting ([a32ee4f](https://github.com/Loop3D/LoopStructural/commit/a32ee4f4973f701638eede4caef1a66e1f4baf2d))
* rename optional data argument to data for consistency ([9b23bd3](https://github.com/Loop3D/LoopStructural/commit/9b23bd30bb1fcf4af4f7ca656b41e88855afcdff))
* return the first feature of a structural frame ([6972062](https://github.com/Loop3D/LoopStructural/commit/6972062ac60f3e97618f26139649db2167bdcfd9))

## [1.6.19](https://github.com/Loop3D/LoopStructural/compare/v1.6.18...v1.6.19) (2025-08-14)


### Bug Fixes

* add default parameters datastructure ([b71c2f8](https://github.com/Loop3D/LoopStructural/commit/b71c2f8d90c240bf2ecc99c70b021c5045f4b5e4))
* add get stratigraphic column cmap to stratigraphic column class ([a03e716](https://github.com/Loop3D/LoopStructural/commit/a03e71646430db6958149b2e3c23be2b5a59e33e))
* allow adding object into model at different indexes ([e1f5cde](https://github.com/Loop3D/LoopStructural/commit/e1f5cde19e033e301db5206ddeb658270bbc5886))
* check id type when creating and use add id to/from dict ([eef12b4](https://github.com/Loop3D/LoopStructural/commit/eef12b4c2569b2957688d385da434270513ccb95))
* convert feature builder to folded feature builder ([fcf4f76](https://github.com/Loop3D/LoopStructural/commit/fcf4f7637245f31487413beb9c1b593efb587a5d))
* don't add basement to the stratigraphic column made from dictionary ([7c98c95](https://github.com/Loop3D/LoopStructural/commit/7c98c95195daedbf4eb33703cd30918d9c7e07a4))
* don't try to access contacts if they are none ([37792a3](https://github.com/Loop3D/LoopStructural/commit/37792a3098ad1574698982d5bd8cee519097e665))

## [1.6.18](https://github.com/Loop3D/LoopStructural/compare/v1.6.17...v1.6.18) (2025-08-04)


### Bug Fixes

* add logging for wavelength guess ([67bebb1](https://github.com/Loop3D/LoopStructural/commit/67bebb1052a2948fb84e95364c06c1de99167667))
* removing typealias ([d1626f5](https://github.com/Loop3D/LoopStructural/commit/d1626f5e75e7edfe582feb2a464450dfd23cc4b2))
* store min/max unit value in unit and keep this up to date ([5af3815](https://github.com/Loop3D/LoopStructural/commit/5af3815c566c8db88dbd0c63b5fc2e00d0d4303d))
* update strartigraphic column/stratigraphic id for new column ([55c303b](https://github.com/Loop3D/LoopStructural/commit/55c303bd3aabdac2a9d51beb143b44f13c67beaa))

## [1.6.17](https://github.com/Loop3D/LoopStructural/compare/v1.6.16...v1.6.17) (2025-07-30)


### Bug Fixes

* add an observer/notify datastructure ([ea4b9f0](https://github.com/Loop3D/LoopStructural/commit/ea4b9f000e4ef43cf349bbb152a4452dc40a4351))
* add get methods for getting specific relationships ([30f01c8](https://github.com/Loop3D/LoopStructural/commit/30f01c837ff776d919d728877a3e4aa4aecd80e3))
* add observer pattern imports for enhanced notification capabilities ([f237d27](https://github.com/Loop3D/LoopStructural/commit/f237d2729073a5c56f2a73d45436f758259c757a))
* add remove fault ([bc4ca17](https://github.com/Loop3D/LoopStructural/commit/bc4ca179ccaeeb59f85c44b5339e50ad7f21cd8a))
* add setter/getter for stratigraphic column ([f771cdd](https://github.com/Loop3D/LoopStructural/commit/f771cdd12d101622d1e1130c94f7726d84fd96bd))
* adding a fault topology datastructure and link with stratigraphic column ([57f362d](https://github.com/Loop3D/LoopStructural/commit/57f362d33bf0fb07a04bae655b4fdbd26ad66a50))
* change strat/fault relationship datastructure to a dictionary with tuple keys ([2a4503b](https://github.com/Loop3D/LoopStructural/commit/2a4503b4973d7e28bf38253741fcae8b75300dbf))
* clarify naming for individual isosurfaces based on input name, don't add isovalue when not needed ([447fb17](https://github.com/Loop3D/LoopStructural/commit/447fb17e64bc42e5154aef36d94c62ad91ac6f78))
* enhance FaultTopology class with notification support for relationship changes ([330c662](https://github.com/Loop3D/LoopStructural/commit/330c6624e8c2f1cfabf4c5391b4c375bf9d23d2d))
* implement Observer pattern with Observable and Disposable classes ([3942920](https://github.com/Loop3D/LoopStructural/commit/394292036bd70aaa9294da5c285f232c1b8bae3e))
* integrate Observable pattern into StratigraphicColumn for enhanced notification support ([6619b80](https://github.com/Loop3D/LoopStructural/commit/6619b80191f3055d656c63f3256ee1e56e034442))
* remove observers before pickle ([c50afdb](https://github.com/Loop3D/LoopStructural/commit/c50afdb0a5ed4b50b4d015e651a3db7b0cbb5162))
* remove raise warning when no weights provided to set_normal_constraints ([33146d5](https://github.com/Loop3D/LoopStructural/commit/33146d570a2733acf51b3841579e72a2b79a1ddd))
* store unitname fault topology instead of group fault ([5a8c0c5](https://github.com/Loop3D/LoopStructural/commit/5a8c0c546f82458cae57b7c30e64f2857efa621a))
* use groupname not group for stratigraphy/fault relationship ([94a563c](https://github.com/Loop3D/LoopStructural/commit/94a563cc7bc88722dc4131f89e2874d109ae586c))

## [1.6.16](https://github.com/Loop3D/LoopStructural/compare/v1.6.15...v1.6.16) (2025-07-21)


### Bug Fixes

* add ability to pass handler to loopstructural logger ([dfebd48](https://github.com/Loop3D/LoopStructural/commit/dfebd487c72830b5e5899660663db93b33e05e0b))
* add new stratigraphic column implementation ([a495bb1](https://github.com/Loop3D/LoopStructural/commit/a495bb16102916b7d53b15b1ab60b939ee2e3440))
* remove initialisation ([66a6f1d](https://github.com/Loop3D/LoopStructural/commit/66a6f1d965c59bc0058e272ae671f4ba1aa90756))
* update stratigraphic column from dictionary ([2f6e0ed](https://github.com/Loop3D/LoopStructural/commit/2f6e0edd2950f497595bfbcfc264d38922a6be0b))

## [1.6.15](https://github.com/Loop3D/LoopStructural/compare/v1.6.14...v1.6.15) (2025-07-08)


### Bug Fixes

* add option to pass a dataframe directly to the create and add methods ([f756054](https://github.com/Loop3D/LoopStructural/commit/f7560545048331137f35172ebd1324895f12faf6))
* add parameter separator to clean up api ([7d8b0b1](https://github.com/Loop3D/LoopStructural/commit/7d8b0b1c2f0c30b4d9119d9438de3ded44b04bd8))
* add visualisation to plot gradient norm ([e774f95](https://github.com/Loop3D/LoopStructural/commit/e774f95934b8fbec84ff37e73a098d36c6620f23))
* adding constant norm interpolators ([29570f1](https://github.com/Loop3D/LoopStructural/commit/29570f1543f2a4aa04efda2fb4d894a7dc1ed9bf))
* allow data to be specified in create and add function. ([f6db4a5](https://github.com/Loop3D/LoopStructural/commit/f6db4a5b6aabc75b7a061014b93288ec91791c57))
* update vector scaling ([0a1e18e](https://github.com/Loop3D/LoopStructural/commit/0a1e18e9740ab8fc72fe2a27231bf6dc8fe2fece))
* updating bounding box project/reproject to just translate to origin ([7606864](https://github.com/Loop3D/LoopStructural/commit/760686432710c92ec2653cafe226012c64c24551))

## [1.6.14](https://github.com/Loop3D/LoopStructural/compare/v1.6.13...v1.6.14) (2025-05-28)


### Bug Fixes

* adjust fault function so that gradient is 0 at the edges ([6fdec1a](https://github.com/Loop3D/LoopStructural/commit/6fdec1a5c1bd0f0ef8a47bfcd3b2ee7623595ab1))
* change surfe import behaviour warning ([5d04a92](https://github.com/Loop3D/LoopStructural/commit/5d04a9298ab1db82b13a71cbb333e68141c29a52))
* default regularisation should be 0.1 for both FDI and PLI. ([7c20b88](https://github.com/Loop3D/LoopStructural/commit/7c20b885550569217c6882dc0bd206f38c1d4f00))
* hide processor import error until its used ([a949ade](https://github.com/Loop3D/LoopStructural/commit/a949ade9fa64366ec6ababb74a81894463cb44e1))
* hide surfe warning ([7aad85b](https://github.com/Loop3D/LoopStructural/commit/7aad85b41ddff799271879a925261eec92689ee8))
* include nsteps in bounding box initialization when creating buffer ([9f859c9](https://github.com/Loop3D/LoopStructural/commit/9f859c984c32cbf266a65a3ae9c6fc56514734a9))
* norm constraint magnitude will be honoured by interpolator ([171f78a](https://github.com/Loop3D/LoopStructural/commit/171f78a98381be5e801ca6d4994e5e57cbc75a60))
* rename azimuthplunge to plungeazimuth ([4605a15](https://github.com/Loop3D/LoopStructural/commit/4605a151f27c3d8071e9e8e206ec34a26d9f373f))
* update tests ([7ca8810](https://github.com/Loop3D/LoopStructural/commit/7ca88101a70155a815ab843e51d1497cd5105467))

## [1.6.13](https://github.com/Loop3D/LoopStructural/compare/v1.6.12...v1.6.13) (2025-04-29)


### Bug Fixes

* add faulted vector calc to lambda feature ([300b575](https://github.com/Loop3D/LoopStructural/commit/300b575181f727b419c6d4e5bf755d4c909a101f))
* remove default initialisation with mutable. ([bda68a7](https://github.com/Loop3D/LoopStructural/commit/bda68a71d22db2fa921f1534e7491dd6995a59a5))
* use fault normals/slip vectors from data if available. ([7264222](https://github.com/Loop3D/LoopStructural/commit/726422220c2e6c0aba49bf0ad7c57b8713e7f585))

## [1.6.12](https://github.com/Loop3D/LoopStructural/compare/v1.6.11...v1.6.12) (2025-04-16)


### Bug Fixes

* add fault function with 0 gradient ([6bd359f](https://github.com/Loop3D/LoopStructural/commit/6bd359f00e921ad5340242af4c02725d46118b2d))
* refactor fault ellipsoid plotter to not use fault builder ([025e286](https://github.com/Loop3D/LoopStructural/commit/025e28635847ae8a5796160ac6784d1b72ddd968))
* vector point visualisation bug where nan values exist ([9c70825](https://github.com/Loop3D/LoopStructural/commit/9c7082535fb5561f79258038e369d53903fa66d2))

## [1.6.11](https://github.com/Loop3D/LoopStructural/compare/v1.6.10...v1.6.11) (2025-04-06)


### Bug Fixes

* colours correct for surfaces ([ea3709a](https://github.com/Loop3D/LoopStructural/commit/ea3709a25c87f4b4f1fbddee77d7a66526298a66))

## [1.6.10](https://github.com/Loop3D/LoopStructural/compare/v1.6.9...v1.6.10) (2025-04-04)


### Bug Fixes

* add fault pitch ([a05d277](https://github.com/Loop3D/LoopStructural/commit/a05d2773663cd6c826fde90ff2465e79de17aa6f))
* Adding colours to surfaces ([2d40563](https://github.com/Loop3D/LoopStructural/commit/2d40563364dfb229fcd3a22b1c5e6e6bc841de6d))

## [1.6.9](https://github.com/Loop3D/LoopStructural/compare/v1.6.8...v1.6.9) (2025-03-30)


### Bug Fixes

* adding copy method for lambda ([6a8d940](https://github.com/Loop3D/LoopStructural/commit/6a8d940a0989a046663dcd3dd4ced6f443d895a6))
* adding from dict method for bb ([71c2dcc](https://github.com/Loop3D/LoopStructural/commit/71c2dccc60855478e04e3fc2f8aa45c969f402a2))
* allow scalar field of feature without a model ([2e36743](https://github.com/Loop3D/LoopStructural/commit/2e36743f1c090bb63efaa6468f2b5388e59fda4a))

## [1.6.8](https://github.com/Loop3D/LoopStructural/compare/v1.6.7...v1.6.8) (2025-02-20)


### Bug Fixes

* updating scaling for plotting ([#219](https://github.com/Loop3D/LoopStructural/issues/219)) ([78ccbd3](https://github.com/Loop3D/LoopStructural/commit/78ccbd3edbb67d49b4c21222bc066fbdd82c4dac))

## [1.6.7](https://github.com/Loop3D/LoopStructural/compare/v1.6.6...v1.6.7) (2025-02-03)


### Bug Fixes

* fault orientation init as empty df rather than nan ([c004d9f](https://github.com/Loop3D/LoopStructural/commit/c004d9f84e65a636faa0566c26797749a42da577))
* update matplotlib cmap for deprecation ([#215](https://github.com/Loop3D/LoopStructural/issues/215)) ([8d7e9f9](https://github.com/Loop3D/LoopStructural/commit/8d7e9f9e6f873befd705473dcacbec0492f85187))

## [1.6.6](https://github.com/Loop3D/LoopStructural/compare/v1.6.5...v1.6.6) (2025-01-23)


### Bug Fixes

* add parent directory to export dir ([3b51f8f](https://github.com/Loop3D/LoopStructural/commit/3b51f8fc398c8d61c182811d4a1478306fd825a3))
* adding interpolator builder. ([169545b](https://github.com/Loop3D/LoopStructural/commit/169545b620046a983a6e2744b80273cc14060f13))
* adding local bb option to isosurfacer ([acc5b95](https://github.com/Loop3D/LoopStructural/commit/acc5b95869accf563ce0d151603f62bc37e9800b))
* adding random hex colour utility function ([62359d4](https://github.com/Loop3D/LoopStructural/commit/62359d46e860b4c944f64081302e2802ee8e3472))
* inactive faults no longer get cropped by unconformities ([4211b9e](https://github.com/Loop3D/LoopStructural/commit/4211b9e118a1f2a0d902974028c553449b0bc10c))
* remove get interpolator (replaced with factory) ([fc5d22a](https://github.com/Loop3D/LoopStructural/commit/fc5d22ade1d2e292c0aef04ccc13f6e69f98c8be))

## [1.6.5](https://github.com/Loop3D/LoopStructural/compare/v1.6.4...v1.6.5) (2024-12-17)


### Miscellaneous Chores

* release 1.6.4 ([f06616f](https://github.com/Loop3D/LoopStructural/commit/f06616f8fac0ca3cfc58377524245952f56e686b))
* release 1.6.5 ([246e48d](https://github.com/Loop3D/LoopStructural/commit/246e48d86a99e9d1e96ab9a2d9567374ffcf8622))

## [1.6.4](https://github.com/Loop3D/LoopStructural/compare/v1.6.4...v1.6.4) (2024-12-17)


### Miscellaneous Chores

* release 1.6.4 ([f06616f](https://github.com/Loop3D/LoopStructural/commit/f06616f8fac0ca3cfc58377524245952f56e686b))
