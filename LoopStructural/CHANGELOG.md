# Changelog

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
