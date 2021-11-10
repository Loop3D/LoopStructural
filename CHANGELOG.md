# Changelog

### [1.3.9](https://www.github.com/Loop3D/LoopStructural/compare/v1.3.8...v1.3.9) (2021-11-10)


### Bug Fixes

* added points argument to faultframe builder ([4588d3c](https://www.github.com/Loop3D/LoopStructural/commit/4588d3c6436666d951dfc84755701b0b7cb95140))
* bugfix thickness=False not working for processor ([45e3f8f](https://www.github.com/Loop3D/LoopStructural/commit/45e3f8f12dbae1a7343eda27819b55d48f321cef))
* fault frame build using scaled vector ([87bc4e2](https://www.github.com/Loop3D/LoopStructural/commit/87bc4e29dff30d234732526c97210c7b3d01ad2b))
* fdi bug was weigthing grad by volume ([7cd6296](https://www.github.com/Loop3D/LoopStructural/commit/7cd6296fd46b23643fbca66f653e18828632f069))
* passing verb argument to pyamg solve ([52c277d](https://www.github.com/Loop3D/LoopStructural/commit/52c277d819ab202bb03e5e7116719dbcc3ea2778))
* reducing default fault buffer from .4 to .2 ([328013e](https://www.github.com/Loop3D/LoopStructural/commit/328013e09210098e38f02ee6fad65f10c54fa067))
* removing _ from pli constraint names ([08473c3](https://www.github.com/Loop3D/LoopStructural/commit/08473c3d90844e74ba7fb68a3b7bc7290ce56ac8))
* typos ([c67fc54](https://www.github.com/Loop3D/LoopStructural/commit/c67fc54b3aedd84a48ff709d32d69ba3c5489f99))
* updating docker image ([ff2becd](https://www.github.com/Loop3D/LoopStructural/commit/ff2becd59c248a461838a5385adcf885edc5f3b5))
* weight can be float or int ([d2c469f](https://www.github.com/Loop3D/LoopStructural/commit/d2c469fc8e5279de116f066ba61d2833dbc16e7f))


### Documentation

* updating docs to use releaseplease changelog ([cbc8690](https://www.github.com/Loop3D/LoopStructural/commit/cbc86909efe5480c6c1d4bcef874c9b2d4dcea1d))

### [1.3.8](https://www.github.com/Loop3D/LoopStructural/compare/v1.3.7...v1.3.8) (2021-11-08)


### Bug Fixes

* adding cython code for getting neighbour pairs ([65a0a44](https://www.github.com/Loop3D/LoopStructural/commit/65a0a44e182adc16b62d6e38bb9cf941b777e3ef))
* adding origin to scalar field ([8bfeb54](https://www.github.com/Loop3D/LoopStructural/commit/8bfeb54e5475e9d9b15fabe9b291373e9747777f))
* bugfix for pli weighting ([1ce7380](https://www.github.com/Loop3D/LoopStructural/commit/1ce7380d4ba000bad004bbcbb308b8889b18fe31))
* faults are more generic ([4fc4627](https://www.github.com/Loop3D/LoopStructural/commit/4fc462736a0d7a8ff5e1ebfa012b116123087933))
* fold builder was adding data twice ([770e666](https://www.github.com/Loop3D/LoopStructural/commit/770e66656f8aa5ea381130ac3822222bddd2a77b))
* fold rotation angle calculated incorrectly ([aea51bf](https://www.github.com/Loop3D/LoopStructural/commit/aea51bf2610a4ca284db258d8f6ece2f4b597ac6))
* generalising support functions for grid and tetra ([761fd03](https://www.github.com/Loop3D/LoopStructural/commit/761fd03671126b3d8200f19126ea12fb5b5e5166))
* gradient constraints in pli weighted by vol ([9b487fd](https://www.github.com/Loop3D/LoopStructural/commit/9b487fd588830ec008facc922241db46a1e832d2))
* kwargs weren't being passed to fold builder ([d5135c5](https://www.github.com/Loop3D/LoopStructural/commit/d5135c5642bb25d34406bbc65349725f1d88b60f))
* matrix assembly based on constraints dict ([20407e4](https://www.github.com/Loop3D/LoopStructural/commit/20407e4a051611cbcdd7a01486a31e2aa9ee45cb))
* need to put exception type to get message ([16c01ea](https://www.github.com/Loop3D/LoopStructural/commit/16c01eaf971769af6b7c6abe664a002faed42551))
* nonlinear will replace discrete class ([86816d1](https://www.github.com/Loop3D/LoopStructural/commit/86816d12d1d8b883f2529d9a114c382a005a2bb3))
* normalise rows of interpolation matrix ([ae315d2](https://www.github.com/Loop3D/LoopStructural/commit/ae315d2122808594863e7903e6be221a5e0a48c1))
* normalise vector for cg ([f8cf221](https://www.github.com/Loop3D/LoopStructural/commit/f8cf221ce3b707e676f16aa04a634fb3ff45d15b))
* plot vector field crashing ([e35f6b2](https://www.github.com/Loop3D/LoopStructural/commit/e35f6b214ccf448022f53e7b0d892b9e022c766b))
* pytest failed ([9136c99](https://www.github.com/Loop3D/LoopStructural/commit/9136c991e7feee6a1e7ffce76745d7c74c9a0357))
* removing sign change for rotation angle ([7ab33fd](https://www.github.com/Loop3D/LoopStructural/commit/7ab33fd63742f4350bf729537bd5b423d6f84274))
* speed up ([e5b53d8](https://www.github.com/Loop3D/LoopStructural/commit/e5b53d87d2cd9bf3afa6d1390bcc6d593f885467))
* speeding up interface constraints ([3f4a845](https://www.github.com/Loop3D/LoopStructural/commit/3f4a845be37fbb85b07260ae8d62dc7b46aaf26c))
* weighting fdi using volume of element ([665e4fe](https://www.github.com/Loop3D/LoopStructural/commit/665e4fe70bbcd310ca3dfb2b957cfe1f25af2b3c))


### Documentation

* adding debugging guide ([540b5ac](https://www.github.com/Loop3D/LoopStructural/commit/540b5ac381bfe3a3b92170f0c1d3b5a2ba5e07de))
* updating to issue form, using scipy for a template ([70b3a67](https://www.github.com/Loop3D/LoopStructural/commit/70b3a67581d2c3802dbed225c6e7057633d9eff3))

### [1.3.7](https://www.github.com/Loop3D/LoopStructural/compare/v1.3.6...v1.3.7) (2021-10-13)


### Bug Fixes

* abutting fault polarity was being calculated inconsistently. ([c77681f](https://www.github.com/Loop3D/LoopStructural/commit/c77681f466070587e235ccaa1ff6f6fd7b87db92))
* adding folded fold frame creates a fold frame not structural frame ([36aa4b3](https://www.github.com/Loop3D/LoopStructural/commit/36aa4b34adfd625766e5ad61e1212f20fabfcdbf))
* call to update feature when initialising rotation angle plotter ([b97f017](https://www.github.com/Loop3D/LoopStructural/commit/b97f017e583f0190cce030b940a8a6888d436a35))
* setting default for viewer to model = none ([8aec0e4](https://www.github.com/Loop3D/LoopStructural/commit/8aec0e493dd16e6725e2eed1ad2b859549580cc1))
* support box is now rescaled for plot ([6723790](https://www.github.com/Loop3D/LoopStructural/commit/672379059ad97a63edf6ec77740d4d37f8a332b9))

### [1.3.6](https://www.github.com/Loop3D/LoopStructural/compare/v1.3.5...v1.3.6) (2021-10-12)


### Bug Fixes

* removing invalid classifiers ([5d8de87](https://www.github.com/Loop3D/LoopStructural/commit/5d8de8782cbfb41ebf229d5e87ac53f6c06f65c3))

### [1.3.5](https://www.github.com/Loop3D/LoopStructural/compare/v1.3.4...v1.3.5) (2021-10-11)


### Bug Fixes

* adding aabb ([512f0a3](https://www.github.com/Loop3D/LoopStructural/commit/512f0a3ec2dddfe0032a03a1d5f24c1a039ddc9e))
* adding gradient option for m2l wrapper ([cdce6e8](https://www.github.com/Loop3D/LoopStructural/commit/cdce6e8887f771d274ef68ad6ad94b25f745a3d5))
* adding option to define custom mesh builder for pli ([b21fa8c](https://www.github.com/Loop3D/LoopStructural/commit/b21fa8c309590c6dd85ea28f0c456a7acb854df2))
* changed to boolean logic for aabb ([f5a5f9b](https://www.github.com/Loop3D/LoopStructural/commit/f5a5f9b811220f534c2b667a652d063ad530294f))
* cleaning up ([262a89d](https://www.github.com/Loop3D/LoopStructural/commit/262a89da9be5aabb3462b5a35dfea830eda6c728))
* incorrect indexing for FDI grad constraints ([d6b8280](https://www.github.com/Loop3D/LoopStructural/commit/d6b82808e8f2716de313a9695293902b499bc27f))
* move cg call to interpolator ([53492a2](https://www.github.com/Loop3D/LoopStructural/commit/53492a2436c6b90cf1ebe7a10283a5d8e5e49749))
* moved generic methods to base class ([23ec788](https://www.github.com/Loop3D/LoopStructural/commit/23ec7889500863bf6f3a43e8b8b5c9a30433a190))
* names kwarg wasn't used for multiple slices ([daebcf0](https://www.github.com/Loop3D/LoopStructural/commit/daebcf04ca0558c2f19ac60646d2e7adca8ef08a))
* pli grad constraint weights were divided by 3 ([cedcffb](https://www.github.com/Loop3D/LoopStructural/commit/cedcffb2577257baa3bd30dc5e38ff7c5451b37a))
* removing old lavavu wrapper, name wasn't ([9226f40](https://www.github.com/Loop3D/LoopStructural/commit/9226f40cdebd2d65ef4d723dfbc5c736a003051c))
* renaming mesh to support for PLI ([2d07317](https://www.github.com/Loop3D/LoopStructural/commit/2d0731740958bd8ef284aa7f7f8d73f6d45734de))
* setup.py codec and filename error ([479ae2b](https://www.github.com/Loop3D/LoopStructural/commit/479ae2b9f9805fb2061da9ae9aa7e5dac9c2f49f))
* unstructured mesh seems to work ([3bd5e24](https://www.github.com/Loop3D/LoopStructural/commit/3bd5e248a2be44812c0d689c60ac0fa33e1c2e31))
* updating setup.py to include metadata for pypi ([fea4317](https://www.github.com/Loop3D/LoopStructural/commit/fea4317b11a7bdbeffadc8c69f9ffc08ada4b3ca))

### [1.3.4](https://www.github.com/Loop3D/LoopStructural/compare/v1.3.3...v1.3.4) (2021-10-05)


### Bug Fixes

* adding loop specific exceptions and project file ([a7664d2](https://www.github.com/Loop3D/LoopStructural/commit/a7664d261f3c6524c3e29f2c0e082ef6b5512813))
* adding non-linear constraint template ([03b9e73](https://www.github.com/Loop3D/LoopStructural/commit/03b9e734dd058dc774dc4c755efc0a7e247e5d92))
* adding non-linear constraint template ([1be8a45](https://www.github.com/Loop3D/LoopStructural/commit/1be8a45473828c90dc56eb23f698e8f20970a609))
* boolean operator in surfe wrapper ([f11816e](https://www.github.com/Loop3D/LoopStructural/commit/f11816e9d1ed92076e5d9719aead0a4e3ad206a0))
* bugfix for gradient constraints ([5fbbb08](https://www.github.com/Loop3D/LoopStructural/commit/5fbbb08733c2fbdf2975a143d996ae18a4711959))
* constant fold axis was referencing fold_frame ([2050b68](https://www.github.com/Loop3D/LoopStructural/commit/2050b68b6ac837cb5113fdce481a78540fc1bdf5))
* intersection visualisation was using the ([1cf531b](https://www.github.com/Loop3D/LoopStructural/commit/1cf531b51a82228133d47d2be6485d742e180e6d))
* method names for FDI/PLI are consistent ([aebba23](https://www.github.com/Loop3D/LoopStructural/commit/aebba23b84c894e1bee75200547514793b32b4ae))
* structural frame weights were being overwritten by ([15e7f23](https://www.github.com/Loop3D/LoopStructural/commit/15e7f23fa4207ce4ded440a99637e0d950010558))

### [1.3.3](https://www.github.com/Loop3D/LoopStructural/compare/v1.3.2...v1.3.3) (2021-09-29)


### Bug Fixes

* kernel parameter for surfe wasn't being applied ([4e99bbf](https://www.github.com/Loop3D/LoopStructural/commit/4e99bbf942467d39d82f1ee8acfe52174b45bb7b))

### [1.3.2](https://www.github.com/Loop3D/LoopStructural/compare/v1.3.1...v1.3.2) (2021-09-28)


### Bug Fixes

* adding missing functions to lavavu wrapper ([8a7b92c](https://www.github.com/Loop3D/LoopStructural/commit/8a7b92ca149c4798a977181c77bbc7338c2d3457))
* adding tangents to surfepy wrapper ([1c0917e](https://www.github.com/Loop3D/LoopStructural/commit/1c0917e4b821383281fb8b507c44bf55fa9f3103))

### [1.3.1](https://www.github.com/Loop3D/LoopStructural/compare/v1.3.0...v1.3.1) (2021-09-28)


### Bug Fixes

* updating workflows to run on edited  ([5b051f4](https://www.github.com/Loop3D/LoopStructural/commit/5b051f4b58e56c08a6a0f5f9366733879be7761c))

## [1.3.0](https://www.github.com/Loop3D/LoopStructural/compare/v1.2.4...v1.3.0) (2021-09-27)


### Features

* updating viewer to use base class ([7c314a0](https://www.github.com/Loop3D/LoopStructural/commit/7c314a057c7217d6e120eb2398a7955e32dfa62f))


### Bug Fixes

* adding base plotting class for easy inheritance ([1a9614d](https://www.github.com/Loop3D/LoopStructural/commit/1a9614d6f3ebf63e073f261de5346a49c326b275))
* adding builders for fold and folded frame ([a6b61fb](https://www.github.com/Loop3D/LoopStructural/commit/a6b61fb0660f37d890377a7f30a4b00b4f13a879))
* adding callback function to builder update ([1eca6f5](https://www.github.com/Loop3D/LoopStructural/commit/1eca6f59e18f3f4b116e6f3c7396860ee26dae47))
* adding callback to model.update ([1fe5fad](https://www.github.com/Loop3D/LoopStructural/commit/1fe5fad5a8c4ec70ef691899e942299481ede189))
* adding callbackfunction to isosurface ([038df51](https://www.github.com/Loop3D/LoopStructural/commit/038df51d1f9fedc4499012350b02e7debbda50ef))
* adding check to see if a feature is valid. ([c9bd3b0](https://www.github.com/Loop3D/LoopStructural/commit/c9bd3b0e1871e4b6d8644ea26dbba4f6ad3cc80c))
* adding exception if fold frame hasn't ([59a9d66](https://www.github.com/Loop3D/LoopStructural/commit/59a9d66ec1a86a22d86fe0556d1401dde8db5afb))
* adding exception when modelling faults with surfe ([15806d7](https://www.github.com/Loop3D/LoopStructural/commit/15806d74a8841cb6532f95f678021c26145279be))
* adding function to generate cmap and range for ([6354712](https://www.github.com/Loop3D/LoopStructural/commit/6354712a5907ac22c4a59d7269760b00779afa74))
* adding lavavu subclass ([b3c042f](https://www.github.com/Loop3D/LoopStructural/commit/b3c042fb87a376d107a3b39e551731f19fc2e6d3))
* adding loop exception for interpolator ([f43a583](https://www.github.com/Loop3D/LoopStructural/commit/f43a583f8e793c2690f45246fb6d2b9c1da08e95))
* adding set item for fold frame to change ([7d975bd](https://www.github.com/Loop3D/LoopStructural/commit/7d975bd3a7bc0dcbf2ecb7336e193f50a05b96d9))
* adding vtk export ([7e7d63a](https://www.github.com/Loop3D/LoopStructural/commit/7e7d63a2044713eba77e8176915d591a3b6a621d))
* changing structural frames to have a setup ([d6fcdea](https://www.github.com/Loop3D/LoopStructural/commit/d6fcdea0efee978f5aba8845c2edd6aee5d22139))
* checking type for vector plotter + adding name ([19c7d5e](https://www.github.com/Loop3D/LoopStructural/commit/19c7d5e7f0bfbccd78aec6fabe52c39e9ecba5b0))
* cmap not being passed for add_isosurface ([491fb34](https://www.github.com/Loop3D/LoopStructural/commit/491fb34ceb16290f08016dac605ac43234baab83))
* creating specific builder for folded foliations ([b40ba9a](https://www.github.com/Loop3D/LoopStructural/commit/b40ba9a90294cc79ddbddc19bc1ded6591f43d07))
* end of day fold frame ([6c49981](https://www.github.com/Loop3D/LoopStructural/commit/6c4998126fb12356cb31e173bd3638c4854efb8b))
* error to catch when not enough data ([bdd93ea](https://www.github.com/Loop3D/LoopStructural/commit/bdd93ea928b428cb4ca2da3313eb7603b78d3afd))
* import vtk exception ([01005b3](https://www.github.com/Loop3D/LoopStructural/commit/01005b3fb7a0483d0b2947f84d482810db1e3617))
* process data doesn't allow bad formatted ([b9f8364](https://www.github.com/Loop3D/LoopStructural/commit/b9f83649cfebe825acc8f5f2e8839867e6a06cb4))
* removing function to calculate intersection, ([9080877](https://www.github.com/Loop3D/LoopStructural/commit/9080877707456a86fe36649218ef27e29a660bf1))
* removing kwargs that have value of None ([ce05b5f](https://www.github.com/Loop3D/LoopStructural/commit/ce05b5f7e142e49dfc7bb5d67bf635bc47349d0d))
* removing normalisation from ([3216d07](https://www.github.com/Loop3D/LoopStructural/commit/3216d07efdca58f445afdfcc06de631c10c8f654))
* removing specific builder for folded fold frame ([9eed793](https://www.github.com/Loop3D/LoopStructural/commit/9eed793cb9be6b62422b7e1237f7e77ea841ff83))
* sections can be painted with a feature ([905eeb1](https://www.github.com/Loop3D/LoopStructural/commit/905eeb1406d3cf0c9110de3bef8288577a37807a))
* transitioning fold frame building into ([7d12b3b](https://www.github.com/Loop3D/LoopStructural/commit/7d12b3b92d4f7f83f61553257716e26375024b06))
* typo in parameter ([917711a](https://www.github.com/Loop3D/LoopStructural/commit/917711a4beeee4b90a203aae0de6a52772766c88))

### [1.2.4](https://www.github.com/Loop3D/LoopStructural/compare/v1.2.3...v1.2.4) (2021-09-07)


### Bug Fixes

* changing ci to run on published ([7c043f0](https://www.github.com/Loop3D/LoopStructural/commit/7c043f0fa70ec51b4b5728db5f0da9432e9cae37))
* pip publish run on release published ([bbc0e38](https://www.github.com/Loop3D/LoopStructural/commit/bbc0e38fcb16e695b0f35bf910cf3f203ef854fc))

### [1.2.3](https://www.github.com/Loop3D/LoopStructural/compare/v1.2.2...v1.2.3) (2021-09-07)


### Bug Fixes

* :green_heart: axis label missing from PLI mesh ([d52ad3f](https://www.github.com/Loop3D/LoopStructural/commit/d52ad3fd50318717fe555c69136322e21fd26442))
* adding release-please action ([be808dd](https://www.github.com/Loop3D/LoopStructural/commit/be808dd59b5503796a44b4835295bdc00c7829f7))
* changing branch for release-please to master ([f0ce0d4](https://www.github.com/Loop3D/LoopStructural/commit/f0ce0d46df050fb903063ede91412db3cb882c00))
* changing version inclusion ([6efca29](https://www.github.com/Loop3D/LoopStructural/commit/6efca29d0c08a23f363314c7fcc5f29c025447c0))
