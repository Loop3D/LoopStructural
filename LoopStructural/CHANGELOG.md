# Changelog

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
