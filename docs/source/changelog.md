# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [1.1.0](https://github.com/Loop3D/LoopStructural/releases/tag/1.1.0) - 2021-06-08

<small>[Compare with 1.0.91](https://github.com/Loop3D/LoopStructural/compare/1.0.91...1.1.0)</small>

### Added
- Addig logger to files that didn't have it ([1f5efeb](https://github.com/Loop3D/LoopStructural/commit/1f5efeb976ea5547f439fcc06a65e1284fbdbb0f) by Lachlan Grose).
- Adding background page for documentation ([4d71586](https://github.com/Loop3D/LoopStructural/commit/4d715863e553abcb894102569573bc43120bac06) by Lachlan Grose).
- Adding option to perturb geometries from m2l ([9c5f3fb](https://github.com/Loop3D/LoopStructural/commit/9c5f3fbd2d733b64c914dd1cad590ec5e975bc1c) by Lachlan Grose).
- Adding scalar field class for visualising regular grid ([98888e6](https://github.com/Loop3D/LoopStructural/commit/98888e6c6dd5d6fbc55483ac92aafeff1d768296) by Lachlan Grose).
- Adding thickness probabilities to m2l function ([fb8a9a2](https://github.com/Loop3D/LoopStructural/commit/fb8a9a2ed8f006fe4c3be186e0d2d640c44d0172) by Lachlan Grose).
- Adding autocalculated min and max to lambda feature ([cf271c4](https://github.com/Loop3D/LoopStructural/commit/cf271c47a6ab40628ee95e24d3b3a34422ed3491) by Lachlan Grose).
- Adding analytical geological feature ([c3c006f](https://github.com/Loop3D/LoopStructural/commit/c3c006f7bfd5b58dac1f0c89c8c417ec44ce4a48) by Lachlan Grose).
- Adding debugging documentation page ([5f0726b](https://github.com/Loop3D/LoopStructural/commit/5f0726b59a2b2998c57df808674b5ee3d99c1953) by Lachlan Grose).
- Adding argument to skip features for debugging m2l ([295cbf3](https://github.com/Loop3D/LoopStructural/commit/295cbf37f700186f9ed816678078f3b47f502eb6) by Lachlan Grose).
- Adding ability to colour faults by displacement magnitude ([3188ba5](https://github.com/Loop3D/LoopStructural/commit/3188ba55f8e910456627a94d5bee9b19347f758e) by Lachlan Grose).
- Adding displacement magnitude to fault print string ([d046abb](https://github.com/Loop3D/LoopStructural/commit/d046abb8c48c412efa59c44c218d295b8da56d4c) by Lachlan Grose).
- Adding function to install equality constraints to builder ([6538614](https://github.com/Loop3D/LoopStructural/commit/6538614b9e198b48197d55f86d367047598c0d21) by Lachlan Grose).
- Added splay/abut as constraint that is evaluated ([7889671](https://github.com/Loop3D/LoopStructural/commit/78896717db98e1e0e172818d1a11baf5753b8e16) by Lachlan Grose).
- Adding series/fault attribute for geological model ([c8d7094](https://github.com/Loop3D/LoopStructural/commit/c8d7094bd331401de904db3684290ba0a1f6e3b8) by Lachlan Grose).
- Adding nsteps attribute to model vis ([9e51b98](https://github.com/Loop3D/LoopStructural/commit/9e51b98458f820b96f6a925f1be4b91658634e61) by Lachlan Grose).
- Adding automatic estimate of the solver tolerance using the bounding box geometry and 1e-10 as a reference ([21bb923](https://github.com/Loop3D/LoopStructural/commit/21bb9236502e663be3826e4bf4dc7a9232d0b0cd) by Lachlan Grose).
- Adding attribute to model which contains list of faults ([680100b](https://github.com/Loop3D/LoopStructural/commit/680100b04334cb75b98b8d84053603a9ca06a96e) by Lachlan Grose).
- Adding fault topology calculation per location in model ([90ed747](https://github.com/Loop3D/LoopStructural/commit/90ed747cf969cea6b307f754c55dea961c872af1) by Lachlan Grose).
- Adding graphviz to development dockerfile ([d119497](https://github.com/Loop3D/LoopStructural/commit/d1194975de4865388951d1d79848c6d2ea6ecc1c) by Lachlan Grose).
- Adding placeholder to rotate interpolation support ([3d65e37](https://github.com/Loop3D/LoopStructural/commit/3d65e370949688964944dd6e269b477295547026) by Lachlan Grose).
- Adding callback_function to add_isosurface this is called with verts, tri, name for each surface and can be used for exporting the surface in custom format or extracting vertices of a surface ([0c0277f](https://github.com/Loop3D/LoopStructural/commit/0c0277ffa918747288e690e8283f24e821a5a6c9) by Lachlan Grose).
- Adding automatic splay fault/abut fault to map2loop model. ([2e2a4df](https://github.com/Loop3D/LoopStructural/commit/2e2a4df48d90e23eb03a04a6fa031869bc35f08e) by Lachlan Grose).

### Fixed
- Fixing bug for flake8 ([15afedc](https://github.com/Loop3D/LoopStructural/commit/15afedcab3d46645a4a601378511de2901cb1d73) by Lachlan Grose).
- Fixed bug in fault topology calculator ([a95ad85](https://github.com/Loop3D/LoopStructural/commit/a95ad8533bd121a7aac6cf04346d3c102be200a3) by Lachlan Grose).
- Fixing topology calculator ([f0038d1](https://github.com/Loop3D/LoopStructural/commit/f0038d1de218f947f4acdfe4030070226dca0cb1) by Lachlan Grose).


## [1.0.91](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.91) - 2021-03-30

<small>[Compare with 1.0.9](https://github.com/Loop3D/LoopStructural/compare/1.0.9...1.0.91)</small>

### Added
- Added opacity to isosurface ([a05385d](https://github.com/Loop3D/LoopStructural/commit/a05385da9425012e8d33d2dcbc879f83d5c55a70) by Lachlan Grose).
- Adding function to estimate fold hinge location ([f3b45fe](https://github.com/Loop3D/LoopStructural/commit/f3b45fe54466ddb147efd11e63f52b0989e27d0b) by Lachlan Grose).
- Added fault splay ([7f0d92f](https://github.com/Loop3D/LoopStructural/commit/7f0d92f3f7de45420493a60ff7b1c0925f36ae00) by Lachlan Grose).
- Adding splay fault to fault builder ([6e7ae27](https://github.com/Loop3D/LoopStructural/commit/6e7ae27932b14c6e7aba940f19b968979b2a8a41) by Lachlan Grose).
- Added ability to save current lavavu view ([1bcdf83](https://github.com/Loop3D/LoopStructural/commit/1bcdf83df7204d1fe08eb6158dd74557e982649f) by Lachlan Grose).
- Adding change log ([e1ca83e](https://github.com/Loop3D/LoopStructural/commit/e1ca83ef35d83f8ae3c5bdcd304bc1002bc8b61c) by Lachlan Grose).

### Changed
- Changed default dip annotation on map ([5019cf6](https://github.com/Loop3D/LoopStructural/commit/5019cf6a66c24c6fd293c0c9dc96c604265f1edf) by Lachlan Grose).

### Fixed
- Fix (hopefully) incompatibility with lavavu versions due to numpy ([ff0cda4](https://github.com/Loop3D/LoopStructural/commit/ff0cda443807489a88cd60cc5c5ef985b74105c8) by Lachlan Grose).
- Fixed nan comparison error for folds ([6bcafc7](https://github.com/Loop3D/LoopStructural/commit/6bcafc706772ea4202e72df5f388f843d5f25596) by Lachlan Grose).
- Fixed bug with empty dataframe ([64f6f1c](https://github.com/Loop3D/LoopStructural/commit/64f6f1cd2ac971c4d9fe2f9f3e5a461888b0b885) by Lachlan Grose).
- Fix - fixing up import issues ([bd4c0f3](https://github.com/Loop3D/LoopStructural/commit/bd4c0f36e74a5dd1c0331fce07d6d0f6237e67f5) by Lachlan Grose).


## [1.0.9](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.9) - 2021-03-16

<small>[Compare with 1.0.89](https://github.com/Loop3D/LoopStructural/compare/1.0.89...1.0.9)</small>


## [1.0.89](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.89) - 2021-03-14

<small>[Compare with 1.0.88](https://github.com/Loop3D/LoopStructural/compare/1.0.88...1.0.89)</small>


## [1.0.88](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.88) - 2021-03-14

<small>[Compare with 1.0.860](https://github.com/Loop3D/LoopStructural/compare/1.0.860...1.0.88)</small>

### Fixed
- Fixing numpy version ([d97c80f](https://github.com/Loop3D/LoopStructural/commit/d97c80ff91ccecf8322338bb1d98596d7b14f013) by Lachlan Grose).


## [1.0.860](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.860) - 2021-03-12

<small>[Compare with 1.0.86](https://github.com/Loop3D/LoopStructural/compare/1.0.86...1.0.860)</small>


## [1.0.86](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.86) - 2021-03-12

<small>[Compare with 1.0.85](https://github.com/Loop3D/LoopStructural/compare/1.0.85...1.0.86)</small>


## [1.0.85](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.85) - 2021-03-12

<small>[Compare with 1.0.84](https://github.com/Loop3D/LoopStructural/compare/1.0.84...1.0.85)</small>


## [1.0.84](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.84) - 2021-03-12

<small>[Compare with 1.0.83](https://github.com/Loop3D/LoopStructural/compare/1.0.83...1.0.84)</small>


## [1.0.83](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.83) - 2021-03-08

<small>[Compare with 1.0.82](https://github.com/Loop3D/LoopStructural/compare/1.0.82...1.0.83)</small>


## [1.0.82](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.82) - 2021-03-08

<small>[Compare with 1.0.81](https://github.com/Loop3D/LoopStructural/compare/1.0.81...1.0.82)</small>


## [1.0.81](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.81) - 2021-03-07

<small>[Compare with 1.0.8](https://github.com/Loop3D/LoopStructural/compare/1.0.8...1.0.81)</small>


## [1.0.8](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.8) - 2021-02-24

<small>[Compare with 1.0.75dev](https://github.com/Loop3D/LoopStructural/compare/1.0.75dev...1.0.8)</small>


## [1.0.75dev](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.75dev) - 2021-02-10

<small>[Compare with 1.74dev](https://github.com/Loop3D/LoopStructural/compare/1.74dev...1.0.75dev)</small>


## [1.74dev](https://github.com/Loop3D/LoopStructural/releases/tag/1.74dev) - 2021-02-10

<small>[Compare with 1.0.72dev](https://github.com/Loop3D/LoopStructural/compare/1.0.72dev...1.74dev)</small>


## [1.0.72dev](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.72dev) - 2021-02-10

<small>[Compare with 1.0.72-dev](https://github.com/Loop3D/LoopStructural/compare/1.0.72-dev...1.0.72dev)</small>


## [1.0.72-dev](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.72-dev) - 2021-02-10

<small>[Compare with 1.0.71-dev](https://github.com/Loop3D/LoopStructural/compare/1.0.71-dev...1.0.72-dev)</small>


## [1.0.71-dev](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.71-dev) - 2021-02-10

<small>[Compare with 1.0.7-dev](https://github.com/Loop3D/LoopStructural/compare/1.0.7-dev...1.0.71-dev)</small>


## [1.0.7-dev](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.7-dev) - 2021-02-10

<small>[Compare with 1.0.6-dev](https://github.com/Loop3D/LoopStructural/compare/1.0.6-dev...1.0.7-dev)</small>


## [1.0.6-dev](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.6-dev) - 2021-02-09

<small>[Compare with 1.0.6](https://github.com/Loop3D/LoopStructural/compare/1.0.6...1.0.6-dev)</small>

### Fixed
- Fix typo ([51a783c](https://github.com/Loop3D/LoopStructural/commit/51a783c07b4c0509b78ab4ee8e5b3c4fc3fe5e8c) by Vincent Picavet).


## [1.0.6](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.6) - 2021-01-13

<small>[Compare with 1.0.5](https://github.com/Loop3D/LoopStructural/compare/1.0.5...1.0.6)</small>

### Fixed
- Fixup formatting in line with flake8 requirements ([1821b91](https://github.com/Loop3D/LoopStructural/commit/1821b91b1b3372ad186e6d6ae82ba95f8b4ab477) by Vincent Fazio).


## [1.0.5](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.5) - 2020-11-23

<small>[Compare with 1.0.4](https://github.com/Loop3D/LoopStructural/compare/1.0.4...1.0.5)</small>

### Added
- Add missing package & improve exporter comments ([7da6f1d](https://github.com/Loop3D/LoopStructural/commit/7da6f1d488acebc4d07cb6e62a9c39ddb9ecfdd9) by Vincent Fazio).

### Fixed
- Fixed bugs with gradient norm does not seem to work as i expected ([d01f7e2](https://github.com/Loop3D/LoopStructural/commit/d01f7e2b252006215242deb4e65d07f336d3b706) by Lachlan Grose).


## [1.0.4](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.4) - 2020-09-02

<small>[Compare with latest](https://github.com/Loop3D/LoopStructural/compare/latest...1.0.4)</small>


## [latest](https://github.com/Loop3D/LoopStructural/releases/tag/latest) - 2020-09-02

<small>[Compare with 1.0.3](https://github.com/Loop3D/LoopStructural/compare/1.0.3...latest)</small>

### Added
- Add documentation ([2de74be](https://github.com/Loop3D/LoopStructural/commit/2de74be7d27d3ca6db95c94468f66251c5356b7a) by Vincent Fazio).
- Add export to file of 6-face cuboid volumes ([f9e4de9](https://github.com/Loop3D/LoopStructural/commit/f9e4de9e1796f166a2ae9baf3067ccbe541ea151) by Vincent Fazio).
- Added 'python3-dev' to ubuntu documentation ([77eab16](https://github.com/Loop3D/LoopStructural/commit/77eab16bc1b71f54c9504a3bae4b4eea276ac4c5) by Vincent Fazio).
- Adding numpy to conda build ([87860c1](https://github.com/Loop3D/LoopStructural/commit/87860c1d983d6cecf82f544690c9b95bc657103c) by Lachlan Grose).


## [1.0.3](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.3) - 2020-07-06

<small>[Compare with 1.0.2](https://github.com/Loop3D/LoopStructural/compare/1.0.2...1.0.3)</small>


## [1.0.2](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.2) - 2020-07-01

<small>[Compare with 1.0.1](https://github.com/Loop3D/LoopStructural/compare/1.0.1...1.0.2)</small>


## [1.0.1](https://github.com/Loop3D/LoopStructural/releases/tag/1.0.1) - 2020-06-24

<small>[Compare with v0.0.1-alpha](https://github.com/Loop3D/LoopStructural/compare/v0.0.1-alpha...1.0.1)</small>


## [v0.0.1-alpha](https://github.com/Loop3D/LoopStructural/releases/tag/v0.0.1-alpha) - 2020-06-23

<small>[Compare with 0.01](https://github.com/Loop3D/LoopStructural/compare/0.01...v0.0.1-alpha)</small>

### Added
- Added documentation and flag to skip faults. currently not adding unconformities... ([71da48f](https://github.com/Loop3D/LoopStructural/commit/71da48f081a8d3296655f93f7f4578f494d29ec0) by Lachlan Grose).
- Adding ability to add references to documentation ([b7b156a](https://github.com/Loop3D/LoopStructural/commit/b7b156a941011fdc4e3e54aec4f3fba049a69cee) by Lachlan Grose).
- Adding more examples ([02c52b1](https://github.com/Loop3D/LoopStructural/commit/02c52b183a3fcab7c0c10ef432ab47070e8ad263) by Lachlan Grose).
- Added stratigraphic column to geological model - means that the scalar field can be divided up into groups and then plotted using discrete blocks. ([9fac829](https://github.com/Loop3D/LoopStructural/commit/9fac829bdf0d724c11cca8e8767116867f5ee1d1) by Lachlan Grose).
- Added ability to plot model voxet ([1e32caa](https://github.com/Loop3D/LoopStructural/commit/1e32caa828fc78abffe5d26bc922789db15adfe8) by Lachlan Grose).
- Added tangent constraints ([6678478](https://github.com/Loop3D/LoopStructural/commit/66784785bf98e55df0c71666372bca318dc00d97) by Lachlan Grose).
- Adding fold geometry attributes to modelling output, removing pca warning ([404e522](https://github.com/Loop3D/LoopStructural/commit/404e522608961a44296066c267c8f93b17e100f2) by Lachlan Grose).
- Add model to feature so that the model geometry can be used to pick mean/max/min values ([06817bd](https://github.com/Loop3D/LoopStructural/commit/06817bd24297b08dff992f3c3890d041839f4017) by Lachlan Grose).
- Added calculate gradient misfit for geological feature, some changes for visualisation ([27be86c](https://github.com/Loop3D/LoopStructural/commit/27be86ce459887a04e4e775e7e064d57353cc9b2) by Lachlan Grose).
- Adding probabilistic geological model class ([400bff2](https://github.com/Loop3D/LoopStructural/commit/400bff22efeb58ca616172d1954850bf990488ba) by Lachlan Grose).
- Adding axis argument to map view allowing the mapview to be added to an existing axis e.g. subplot ([b56bf0a](https://github.com/Loop3D/LoopStructural/commit/b56bf0a4cb73f380998edc29b52679190ff9c917) by Lachlan Grose).
- Adding files ([0ba6245](https://github.com/Loop3D/LoopStructural/commit/0ba6245660a370c33095559949e382e0125ac461) by Lachlan Grose).
- Adding wheel ([8ec7d1b](https://github.com/Loop3D/LoopStructural/commit/8ec7d1b5a5099fba2315274e38c18dd8d1033f14) by Lachlan Grose).

### Changed
- Changed print statements to logger.info ([6cc3d32](https://github.com/Loop3D/LoopStructural/commit/6cc3d32b32b82a573a02b603a68e71bb25ec1aea) by Lachlan Grose).
- Changed nodes to @property so they aren't stored permanently ([0cdbe79](https://github.com/Loop3D/LoopStructural/commit/0cdbe79ca710c8821e139dc96af3c500f0553b09) by Lachlan Grose).
- Changed surfe default to single surfaces ([5716d79](https://github.com/Loop3D/LoopStructural/commit/5716d799bb92e715682d7a89e26741caccfb4b8a) by Lachlan Grose).

### Fixed
- Fixed typo in adding features ([8980972](https://github.com/Loop3D/LoopStructural/commit/8980972a3a218da64f9aa03739153def8499d9b7) by Lachlan Grose).
- Fixing typo ([bd09493](https://github.com/Loop3D/LoopStructural/commit/bd09493196b95cf26325124b4c2dd5b80775fa1b) by Lachlan Grose).
- Fixing up misc functions for probabilistic folds/faults ([e502efc](https://github.com/Loop3D/LoopStructural/commit/e502efc88a9f9f233ceb212f9f71f3baa4b9ef89) by Lachlan Grose).
- Fixing up some small errors. and ran tests ([1e82d92](https://github.com/Loop3D/LoopStructural/commit/1e82d92b1896ad11934655a905d208f17595e83e) by Lachlan Grose).
- Fixed bug in fdi ([d247fcd](https://github.com/Loop3D/LoopStructural/commit/d247fcd226e5e13f2b6584bd3cbf2af24cd8a5a6) by Lachlan Grose).

### Removed
- Remove metadata changes ([3e4185d](https://github.com/Loop3D/LoopStructural/commit/3e4185d74c38f4f4ef6ee105871429f1ff815b10) by Vincent Fazio).


## [0.01](https://github.com/Loop3D/LoopStructural/releases/tag/0.01) - 2020-04-01

<small>[Compare with first commit](https://github.com/Loop3D/LoopStructural/compare/c36f54d6419af8ac79c61b46e21ae5fd63c5a5be...0.01)</small>

### Added
- Adding surfe wrapper for surfe interpolator ([76571fd](https://github.com/Loop3D/LoopStructural/commit/76571fd52ef54ead701058592fc7231bcb3dac1b) by Lachlan Grose).
- Adding in lsqr as a solver ([774873c](https://github.com/Loop3D/LoopStructural/commit/774873ca25b88b146faa1cb3d14ba1b309083fd6) by Lachlan Grose).
- Adding wheel ([3a403e7](https://github.com/Loop3D/LoopStructural/commit/3a403e70a72e7068d4a4cf505051d6ddbf374256) by Lachlan Grose).
- Added lsqr and a few lines of documentation ([206fc5a](https://github.com/Loop3D/LoopStructural/commit/206fc5a1bd00978ae0ffce713c08ff7bebc02050) by Lachlan Grose).
- Adding duplex example ([ea55774](https://github.com/Loop3D/LoopStructural/commit/ea55774e33bfe88529375f9c1ac253b2b667a6bc) by Lachlan Grose).
- Adding deleted line back ([e2ebf31](https://github.com/Loop3D/LoopStructural/commit/e2ebf31a4e2a5b741d520a55c49383a84acc84d4) by Lachlan Grose).
- Adding notebook ([ab18ff0](https://github.com/Loop3D/LoopStructural/commit/ab18ff07fec59dedbb5da5c4df183ef39412c415) by Lachlan Grose).
- Adding per data point weighting to interpolation ([6fa8226](https://github.com/Loop3D/LoopStructural/commit/6fa8226f734ca13249daf003054a31bf1ad28f80) by Lachlan Grose).
- Adding banner to docs ([e5e80ed](https://github.com/Loop3D/LoopStructural/commit/e5e80ed297865e5033533abde77dd224b4aea369) by Lachlan Grose).
- Adding in unconformities ([a7a143c](https://github.com/Loop3D/LoopStructural/commit/a7a143c614892a5d7b1b2332594119bbdb806e29) by Lachlan Grose).
- Adding dsi back in ([14ce0e1](https://github.com/Loop3D/LoopStructural/commit/14ce0e145e5f6e96b50ed30384cf60105b6dd716) by Lachlan Grose).
- Adding compile directives for cython ([8b9c298](https://github.com/Loop3D/LoopStructural/commit/8b9c298e3c9975fef7d3b457abd335b56ca194e4) by Lachlan Grose).
- Adding some logging ([0602c1e](https://github.com/Loop3D/LoopStructural/commit/0602c1eb2669a8cec5ad6c12b25e503d4f6fd981) by Lachlan Grose).
- Adding test for refoded folds ([fbabb1c](https://github.com/Loop3D/LoopStructural/commit/fbabb1cd239cedebfb85fbcfdeedc32273d707b6) by Lachlan Grose).
- Added language level to cython ([0ad94ef](https://github.com/Loop3D/LoopStructural/commit/0ad94ef012a9c70b3a2e36a95232524d15198faf) by Lachlan Grose).
- Added a lot of logging functions ([83b7389](https://github.com/Loop3D/LoopStructural/commit/83b738954ae0484f5da4d1a26d5a7b4424cb0284) by Lachlan Grose).
- Adding polyfold model class ([0b9dac1](https://github.com/Loop3D/LoopStructural/commit/0b9dac1d85c97eafe1cb4499cfb1cd851cf0903f) by lachlan).
- Adding pandas to requirements ([d8896ef](https://github.com/Loop3D/LoopStructural/commit/d8896ef2d4c9daa85e586282b345b6a67cf78f1a) by Lachlan Grose).
- Added dummy test ([dccd804](https://github.com/Loop3D/LoopStructural/commit/dccd80402817ca792d77ca5e7fe13e58943c3b38) by Lachlan Grose).
- Adding apt installs to yml ([e3dd426](https://github.com/Loop3D/LoopStructural/commit/e3dd4260fed667f9effe4cec985c43b8f68e8afb) by lachlangrose).
- Added some helper functions for visualisation, graben example ([52a14bc](https://github.com/Loop3D/LoopStructural/commit/52a14bc202a6540c87c1eef21bf5c224623fdff2) by Lachlan Grose).
- Adding data for duplex ([3630c74](https://github.com/Loop3D/LoopStructural/commit/3630c745b80fdb8d70972a670cd4b04c392edba7) by Lachlan Grose).
- Adding missing nbs ([86189d3](https://github.com/Loop3D/LoopStructural/commit/86189d3a15f8f7701d1ed4b6d3c4dae5aacc0549) by Lachlan Grose).
- Added fold axis for fold rotation angle calculator back in ([5c20043](https://github.com/Loop3D/LoopStructural/commit/5c20043d0e28f1dab7587170c5bea7a94a36204c) by Lachlan Grose).
- Adding faulted intrusion dataset and example ([158ab4d](https://github.com/Loop3D/LoopStructural/commit/158ab4d3ddf3a4cf73634d98cfc1da2191b98862) by Lachlan Grose).
- Added faults to evaluate_gradient ([6adc0fd](https://github.com/Loop3D/LoopStructural/commit/6adc0fdc5fc4b0405f7419bfeeef2ec2d72da6ca) by lachlan).
- Added folded foliation to the model class, now need to add refolded folds ([5ebaf16](https://github.com/Loop3D/LoopStructural/commit/5ebaf168f40c3856b1c23f7ba45d57d62a81f796) by Lachlan Grose).
- Adding demo data for claudius model adding dataset module that can load the demo datasets adding manifest.in to specify which files to copy over when installing updated setup.py to copy module data ([da339c1](https://github.com/Loop3D/LoopStructural/commit/da339c1ac312c27acb9b4ac72f115e6aef2ed909) by Lachlan Grose).
- Added azimith and dip to dataframe parser ([64e8ba2](https://github.com/Loop3D/LoopStructural/commit/64e8ba262c8ab4123f64b84de0d9534454e83a78) by Lachlan Grose).
- Adding change summary ([afbadc7](https://github.com/Loop3D/LoopStructural/commit/afbadc73932724aacbb3d3b86d389627bce4a06d) by lachlan).
- Adding numpy tetmesh ([3ad7a35](https://github.com/Loop3D/LoopStructural/commit/3ad7a35c2e24c71ad7eca1714acf61f59d7638ab) by Lachlan Grose).
- Adding grad constraint to fdi + some docs ([645a525](https://github.com/Loop3D/LoopStructural/commit/645a52570878524fd34a4344adbf2a4b94f26d06) by Lachlan Grose).
- Added gradient direction constraint as well as gradient norm constraint ([8db5aea](https://github.com/Loop3D/LoopStructural/commit/8db5aea54033790b02b4f546957c6957bb981fee) by Lachlan Grose).
- Adding data ([ffda7f4](https://github.com/Loop3D/LoopStructural/commit/ffda7f46aef99da9ee82aabd1ae4e9e2aaeea29e) by Lachlan Grose).
- Adding returns to geological model ([c2d3d76](https://github.com/Loop3D/LoopStructural/commit/c2d3d761bd11e1aac06c61b75064eca3b48330f8) by Lachlan Grose).
- Added tutorial for external solver made faults so that the data points are faulted when the feaature is created ([34688d6](https://github.com/Loop3D/LoopStructural/commit/34688d679caa80c23f73d2e19abb6351d7aff0fe) by Lachlan Grose).
- Adding in some modifications to solvers ([6e3d081](https://github.com/Loop3D/LoopStructural/commit/6e3d08171e095b66d18ec9f949fec652f40878ca) by Lachlan Grose).
- Added logger rather than using print debugs ([d0563d9](https://github.com/Loop3D/LoopStructural/commit/d0563d9cdaaee1fd0e842a4ac1e0a8257d451909) by Lachlan Grose).
- Added more documentation to interpolators, also added ability for finite difference to use a dictionary of operators rather than just the default config ([93db9e1](https://github.com/Loop3D/LoopStructural/commit/93db9e14ee065c4c20855a8da76cf67afb303b93) by Lachlan Grose).
- Adding mapview ([044dc71](https://github.com/Loop3D/LoopStructural/commit/044dc714948ab35d976cfcd631a7297317a2a302) by lachlan).
- Adding missing files ([fd8a398](https://github.com/Loop3D/LoopStructural/commit/fd8a39871b358d64623072d00ecd2cdf59cb760c) by Lachlan Grose).
- Adding features for fold rotation and fault displacement ([1f2f5f1](https://github.com/Loop3D/LoopStructural/commit/1f2f5f181bba1cda2f08d3f0c438f5d8be398d4a) by lachlan).
- Adding tutorial nbs ([446623e](https://github.com/Loop3D/LoopStructural/commit/446623ec6e549ab8fd95f35ce8559505db638a5d) by Lachlan Grose).
- Added isosurfacing on another object for a geological feature ([4744055](https://github.com/Loop3D/LoopStructural/commit/47440553a88fa2ee1e856364fcc1677f70d62ef6) by Lachlan Grose).
- Adding constant fold axis and fixing tutorial for polydeformed fold ([60deda5](https://github.com/Loop3D/LoopStructural/commit/60deda5dcf3af1fa9cb0e7f6769606769b3e72af) by Lachlan Grose).
- Added scalar field plot option ([f60fe06](https://github.com/Loop3D/LoopStructural/commit/f60fe068e5c62e812d5b4a691b868833959a7955) by Lachlan Grose).
- Adding notebooks ([a18dff0](https://github.com/Loop3D/LoopStructural/commit/a18dff06e26ffdb1082a8a241d2ed0f2f24f9e35) by Lachlan Grose).
- Adding documentation figures ([a65f9ad](https://github.com/Loop3D/LoopStructural/commit/a65f9ad12d2d559f5849e05fa4e4a19d467a5e7f) by Lachlan Grose).
- Adding claudius exaple ([351d8e6](https://github.com/Loop3D/LoopStructural/commit/351d8e66a573ff12cbeb2fd299bd4a85879f0f22) by Lachlan Grose).
- Adding claudius data to examples and updating structured grid to move to origin ([2b92a8b](https://github.com/Loop3D/LoopStructural/commit/2b92a8beb2d84f652dc12daca8e40ab53b3af68c) by Lachlan Grose).
- Added mesh kwarg for structuralframe to be backwards compatible ([328fca8](https://github.com/Loop3D/LoopStructural/commit/328fca8604dbbc8ff01c2a8d7e72a6bc97382111) by Lachlan Grose).
- Added orthogonality constraint for finite difference ([291324f](https://github.com/Loop3D/LoopStructural/commit/291324fcae070078f22af406dbf9430c8829f59b) by Lachlan Grose).
- Added equality constraints for constrained least squares ([df22a2d](https://github.com/Loop3D/LoopStructural/commit/df22a2de5c79effd180fe984e4ae08ee6d443570) by Lachlan Grose).
- Adding grid ([5ddd113](https://github.com/Loop3D/LoopStructural/commit/5ddd11393cdd78d6a5957bae128d9f24edb5bf1c) by Lachlan Grose).
- Adding makefile ([7eaa2a1](https://github.com/Loop3D/LoopStructural/commit/7eaa2a13d6398bf657185ff1867d977100c2b5e9) by Lachlan Grose).
- Adding documentation source and map_viewer file ([5317220](https://github.com/Loop3D/LoopStructural/commit/53172204b9b8b4ba88b9693afe406c61df34500d) by Lachlan Grose).
- Added an update function to the geologicalfeature so that faulted features are correctly updated ([c0e10bd](https://github.com/Loop3D/LoopStructural/commit/c0e10bdbfae68bdef5a219c81e38bae9987e5424) by Lachlan Grose).
- Added update function for faulted geological feature so you can update the faulted feature if the parent feature or fault are reevaluated. example on noddy working well :) ([25f9da3](https://github.com/Loop3D/LoopStructural/commit/25f9da3a94dc4e4058c13232d645e9ef48f21581) by lachlan).
- Adding correct notebook ([9621a9c](https://github.com/Loop3D/LoopStructural/commit/9621a9c0e4106b1170d6cbceaab81ada88eeef19) by Lachlan Grose).
- Added check for fold frame coordinates during builder ([8e27197](https://github.com/Loop3D/LoopStructural/commit/8e27197fddc6140d269033764e54237c39c71bb7) by Lachlan Grose).
- Added from_strike_dip and from plunge constructors for geologicla points. also working on transforming points for fault interpolation ([56cae64](https://github.com/Loop3D/LoopStructural/commit/56cae6461b89647dad51ccfc6df1565d785c135a) by Lachlan Grose).
- Added object to define fault displacement using a combination of 1d functions that can be cubic splines or any other callable function ([c41ed5f](https://github.com/Loop3D/LoopStructural/commit/c41ed5f43fb9ef5eb0d10ba7cc878005cd36812c) by Lachlan Grose).
- Added option to plot surface normals ([090a1a4](https://github.com/Loop3D/LoopStructural/commit/090a1a4262341163bb2627d131be12f9cb831fd5) by Lachlan Grose).
- Adding faulted geological feature code, not working 100% yet ([64f0607](https://github.com/Loop3D/LoopStructural/commit/64f06078ff30d6019cbc1b43482e18a60abea0b5) by Lachlan Grose).
- Adding functions to geological feature ([47ef2e0](https://github.com/Loop3D/LoopStructural/commit/47ef2e0f6ee0cb6776e14a82b0e0c4a66aee3b33) by Lachlan Grose).
- Adding plot options for vector field ([e14afdf](https://github.com/Loop3D/LoopStructural/commit/e14afdf808349ab4e38e3476e1cd9818820b1361) by Lachlan Grose).
- Added marching tetrahedron cython functions ([0a88436](https://github.com/Loop3D/LoopStructural/commit/0a88436474ebf8497a2834615f25a585fb862d28) by Lachlan Grose).
- Added plotting method to tet_mesh ([f039ade](https://github.com/Loop3D/LoopStructural/commit/f039adee5a40e24196f43750f8b2b3a90640b9da) by Lachlan Grose).
- Added files for folding and an example ([e484583](https://github.com/Loop3D/LoopStructural/commit/e484583c72d484f31645cf7e330976e2ddd6df2b) by Lachlan Grose).

### Changed
- Changed tetra ordering as per roys code ([c221c5a](https://github.com/Loop3D/LoopStructural/commit/c221c5afeaf9aa711375eff76d7cf682cdf077c6) by Lachlan Grose).
- Changed model viewer isosurfacing to only compute scalarfield once. has changed api for vis a lot ([a247da4](https://github.com/Loop3D/LoopStructural/commit/a247da44b7098b1696593bfea99026ff59cecf41) by Lachlan Grose).
- Changed regions from a property on the mesh to a lambda function of a position or array of positions ([1d70618](https://github.com/Loop3D/LoopStructural/commit/1d7061805095030d3446860eb6d657f099682ad9) by Lachlan Grose).
- Changed some imports and debugging ([ea43f28](https://github.com/Loop3D/LoopStructural/commit/ea43f28b6710f13d7e2295fa09a44a0245f0d076) by Lachlan Grose).
- Changed plot to add in visualisation. ([a07bdef](https://github.com/Loop3D/LoopStructural/commit/a07bdef80a26ef468c2c2d15acbe0e41bf122611) by Lachlan Grose).
- Changed imports to relative ([826b83f](https://github.com/Loop3D/LoopStructural/commit/826b83fdffa46a5e98b8aa7c99068ed8ad224fcb) by Lachlan Grose).
- Changed how data are assigned to structural frame, gave structural frame [] operator to access feature ([ece58f5](https://github.com/Loop3D/LoopStructural/commit/ece58f5ab8c9c47dffb97f63a09d88cda7fc2715) by Lachlan Grose).
- Changed constant gradient to use normalised face normal. this should fix up problems with cg values being too high ([5c86bc5](https://github.com/Loop3D/LoopStructural/commit/5c86bc59ab16b73cd678c8a6c91880aa754bf695) by Lachlan Grose).
- Changed long to long long in cython definitions inline with int64 ([eb04cf3](https://github.com/Loop3D/LoopStructural/commit/eb04cf31d1d0131ffe48fe2fe3bed570708fee33) by Roy Thomson).
- Changed how pli weights are implemented - regularisation is weighted by the area of the shared face. points and gradient points are weighted by the volume of the tetrahedron. ([8bdff99](https://github.com/Loop3D/LoopStructural/commit/8bdff99e9d8f9ab68619aa87fccf5d5d1709f1b9) by Lachlan Grose).
- Changed boolean for faultedfeature gradient evaluation ([64730f3](https://github.com/Loop3D/LoopStructural/commit/64730f3bd76d5cb44258e8d79d4a6d83dfe85bc5) by Lachlan Grose).
- Changed storing of data for structural frame builder to use a list of lists rather than a dict, this can be passed to the fault segment more easily. in the future could consider removing the 'gx','gy' tags and just use frame index. ([06c2900](https://github.com/Loop3D/LoopStructural/commit/06c29005e2e13b93a6a8545575152999e768acf4) by Lachlan Grose).
- Changed interpolation weights to be a dictionary stored by the interpolation object that can be updated when calling setup interpolator. this means the interpolator can be rerun using the same weights ([8513693](https://github.com/Loop3D/LoopStructural/commit/8513693693eafbd5a8e4efedd9b2ae95eecc17a7) by Lachlan Grose).
- Changed structural frame so that the gz coordinte can be analytical ([5dda636](https://github.com/Loop3D/LoopStructural/commit/5dda6362371cebc5007be65d6e40bcc2960473ff) by Lachlan Grose).
- Changed name of tests folder to examples ([8b47c03](https://github.com/Loop3D/LoopStructural/commit/8b47c0340d19fd706f0180fb5addf241c0c06bfe) by Lachlan Grose).
- Changed docker file to use same base as lavavu and fme init now requires importing submodules manually ([c7c1cd2](https://github.com/Loop3D/LoopStructural/commit/c7c1cd220ccdb80c97c1f7de3a06feb219faf09e) by Lachlan Grose).

### Fixed
- Fixing import after refactor ([57b958e](https://github.com/Loop3D/LoopStructural/commit/57b958ef6bd5bbc2d5a811b44b1d8f9e8d2b117e) by Lachlan Grose).
- Fixing vector pairs. calculate strike vector and then dip vector is orthogonal to normal and strike. ([03954c8](https://github.com/Loop3D/LoopStructural/commit/03954c827b73ee3c9990de5f251886c03926b080) by Lachlan Grose).
- Fixed indexing for features and removed unneeded function from mapview ([0064fdf](https://github.com/Loop3D/LoopStructural/commit/0064fdfa7434dfd9325357ad2189358b88876302) by Lachlan Grose).
- Fixing up a few small bugs ([b3934b1](https://github.com/Loop3D/LoopStructural/commit/b3934b1d04e5419ac10d7e0287a43b7327d6a599) by Lachlan Grose).
- Fixing a few bugs ([f5eb7c4](https://github.com/Loop3D/LoopStructural/commit/f5eb7c48fef18e9284dde616d8bae485f02f11dc) by Lachlan Grose).
- Fixed fold frame for pli, fdi is still broken ([fc0edb3](https://github.com/Loop3D/LoopStructural/commit/fc0edb37017c7df941bd0455c1adb0072a003c74) by Lachlan Grose).
- Fixed regular tetmesh for interpolating single foliation - some issues with building a fold frame ([3e98feb](https://github.com/Loop3D/LoopStructural/commit/3e98feb3d93932c953545daf2ac6c391dad29fec) by Lachlan Grose).
- Fixed bug in masks ([1da3707](https://github.com/Loop3D/LoopStructural/commit/1da3707953b8c14e3accc549ce797ab7e07a3488) by Lachlan Grose).
- Fixed up out of bounds tetmesh issue ([0c9ba6d](https://github.com/Loop3D/LoopStructural/commit/0c9ba6dc7c7d1226915c8059df0b0d7ed7ca62c7) by Lachlan Grose).
- Fixing doc ([a535650](https://github.com/Loop3D/LoopStructural/commit/a535650c004ee10086aec8c48be8a4ac4c1fc966) by lachlan).
- Fixed bug in feature builder and adding ability for different neighbour masks for structured grid. ([52b9ad8](https://github.com/Loop3D/LoopStructural/commit/52b9ad876e39eb8964ebadbb22aab784bfb34fc5) by Lachlan Grose).
- Fixing up fold bugs ([ac4faad](https://github.com/Loop3D/LoopStructural/commit/ac4faad1966d1372fd20559aba0f50d533726978) by Lachlan Grose).
- Fixed lint issue, testing ideas for graben model with +ve and -ve displacements.... should work but still a bug ([7a7847e](https://github.com/Loop3D/LoopStructural/commit/7a7847ee338d6ba4a3317c0ec5c734d85bed2c2e) by lachlan).
- Fixing up some docstrings ([9f06535](https://github.com/Loop3D/LoopStructural/commit/9f06535b2be9510328a117550052ca177bbced4e) by lachlan).
- Fixing for flake8 ([7eead50](https://github.com/Loop3D/LoopStructural/commit/7eead506c76a5bed0a3c40ab3a4598728c3f8ff4) by Lachlan Grose).
- Fixing some bugs ([7e9f89b](https://github.com/Loop3D/LoopStructural/commit/7e9f89b348a208dcba6fc586404d8ce77da8c4c6) by Lachlan Grose).
- Fixed nodes array on structured grid to give actual nodes, meaning hard constraints work. ([d431a48](https://github.com/Loop3D/LoopStructural/commit/d431a4846d391a1fafb444e8c9b07ea7f4441d98) by lachlan).
- Fixing visualisation for faults ([559573e](https://github.com/Loop3D/LoopStructural/commit/559573ecdc69d9daa5a62024a121183d8a4a24c8) by Lachlan Grose).
- Fixed problem with fdi ([c498eb8](https://github.com/Loop3D/LoopStructural/commit/c498eb80f20e67316110914a528c3754493c2b38) by Lachlan Grose).
- Fixed fold norm problem and adding some figures ([e4674bf](https://github.com/Loop3D/LoopStructural/commit/e4674bfff11a2362db7d6e8864f9814dbc5846d0) by Lachlan Grose).
- Fixing up some warnings ([77e3a5d](https://github.com/Loop3D/LoopStructural/commit/77e3a5d0563762be3316e5983ec2b6f0c0484cb8) by Lachlan Grose).
- Fixing up faults for fdi ([8121208](https://github.com/Loop3D/LoopStructural/commit/8121208f1c9130076f17f344509d6ac8ae8fbae6) by Lachlan Grose).
- Fixing up some of the examples so they can be run in batch. most are working, ellipsoid faults need to be checked. ([ff2e180](https://github.com/Loop3D/LoopStructural/commit/ff2e180c4d6f57ead47735e8f9cdd71c33e32c1f) by Lachlan Grose).
- Fixing up some of the examples so they can be run in batch. ([8324073](https://github.com/Loop3D/LoopStructural/commit/8324073bee5e2da27931a2b7aed5fac31c87d746) by Lachlan Grose).
- Fixed bug with gradient not being applied ([666b774](https://github.com/Loop3D/LoopStructural/commit/666b7749bc0c7ed87fed1fd669e58e785709f7f3) by lachlan).
- Fixed pli constraint - point values were not added correctly ([b43ead3](https://github.com/Loop3D/LoopStructural/commit/b43ead3909a9bfff5ea86079934bc1965b85cf59) by lachlan).
- Fixing up fold constraints. resetting the interpolator when rerunning was not good for fold... ([33e6576](https://github.com/Loop3D/LoopStructural/commit/33e65764ae58a7d16778edf0fe8775d955ee536d) by Lachlan Grose).
- Fixed scalar field plot for cuboid ([7454350](https://github.com/Loop3D/LoopStructural/commit/7454350001207099f5cdb5aea2384ce8151c99d4) by lachlan).
- Fixing up finite difference interpolation ([e369ed1](https://github.com/Loop3D/LoopStructural/commit/e369ed15d39d20e7e5fb0cb047c79737f796f89d) by Lachlan Grose).
- Fixed up fdi so that it can be used to evaluate the function on pointsets. ([f92b293](https://github.com/Loop3D/LoopStructural/commit/f92b293236e41b459ed182d09f59887948f69e1e) by Lachlan Grose).
- Fixed examples for new syntax ([ae0152f](https://github.com/Loop3D/LoopStructural/commit/ae0152f14ea9da2af3a87b636d630a29773efea1) by lachlan).
- Fixed up docker and added a notebook dir back in for rabii ([f7b2829](https://github.com/Loop3D/LoopStructural/commit/f7b28296a826d1f724dd19a462173e4de72d0d36) by lachlan).
- Fixed triangle normals for visualisation ([97608dc](https://github.com/Loop3D/LoopStructural/commit/97608dcece2c3aeb36afc62fc6df149ad0438982) by Lachlan Grose).
- Fixed bug with gradient points and added visualisation of input data ([67c65a3](https://github.com/Loop3D/LoopStructural/commit/67c65a38c7339eef6fb65cdbb8602cdd9c76439a) by Lachlan Grose).
- Fixed problem with faults by adding buffer rather than selecting surrounding elements ([e9f1403](https://github.com/Loop3D/LoopStructural/commit/e9f14039aefadad50053745ec44a262464ce519c) by Lachlan Grose).
- Fixed isosurfaces - note api has changed basic example works ([144b560](https://github.com/Loop3D/LoopStructural/commit/144b5607c773a7accd78c8eed595d73c40b61111) by Lachlan Grose).
- Fixing isosurface ([200c972](https://github.com/Loop3D/LoopStructural/commit/200c97279e9119c84006278bc01bff7d8868cf73) by Lachlan Grose).
- Fixing up docker file, some issue with linking to x display for lavavu... maybe this is an issue with using docker? ([731ed67](https://github.com/Loop3D/LoopStructural/commit/731ed672831e5eeae6e72c954611059e49cb471e) by Lachlan Grose).
- Fix problem finding wrong element for point because of pca transformation ([820b3f0](https://github.com/Loop3D/LoopStructural/commit/820b3f06bc3a74f6a6c4f1786b9f80e850c8689d) by Lachlan Grose).

### Removed
- Removed cython from setup.py ([63642a9](https://github.com/Loop3D/LoopStructural/commit/63642a9358f633f0aec99220ae98370e8899bb28) by lachlan).
- Removed fold axis direction field from fold limb roation calc ([48a2c48](https://github.com/Loop3D/LoopStructural/commit/48a2c48e5f0fc270e0339d230bfa38e28b3e2696) by Lachlan Grose).
- Removed print from cython code, removing unused file. adding dataset for laurent2016 ([efde723](https://github.com/Loop3D/LoopStructural/commit/efde72332e0fb15cadc0e9d25a93e188609960f8) by Lachlan Grose).
- Removed accidental copy and paste ([98ff72c](https://github.com/Loop3D/LoopStructural/commit/98ff72c48f9cbc63c12d4010a245c86caf1d95e7) by Lachlan Grose).
- Removed old doc from structured grid, moving some imports to init files for easier use ([368e686](https://github.com/Loop3D/LoopStructural/commit/368e68695694cfd2853febdc070ee58d8b8355ff) by Lachlan Grose).
- Removed changes to cg, doesn't work. need to try weighting the gradient control points by a representative radius... but how to calculate this? is it the distance to the nearest neighbour? ([243fa7c](https://github.com/Loop3D/LoopStructural/commit/243fa7cad9717baa9fe1e07d3e84b431f897098a) by lachlan).
- Removed counter for cg/cgp ([f5f9cd9](https://github.com/Loop3D/LoopStructural/commit/f5f9cd90bc87108ff361905bc6de08a172b93785) by Lachlan Grose).
- Removed debug print ([90992c0](https://github.com/Loop3D/LoopStructural/commit/90992c03559d21515138487ea7ac4d68c4b1a864) by lachlan).
- Removed stray charater ([85fc424](https://github.com/Loop3D/LoopStructural/commit/85fc4243af4bd087685d1c3c49b4da00d1c4be42) by Lachlan Grose).
- Removed unused dsi functions ([0a35545](https://github.com/Loop3D/LoopStructural/commit/0a3554566ebe562254cc796d913b362d9d90d62a) by Lachlan Grose).


