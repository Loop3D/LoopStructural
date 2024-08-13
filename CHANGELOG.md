# Changelog

## [1.6.2](https://github.com/Loop3D/LoopStructural/compare/v1.6.1...v1.6.2) (2024-08-06)


### Bug Fixes

* extra import ([7d10434](https://github.com/Loop3D/LoopStructural/commit/7d10434eb11631fa501275c14d617ed014f092a7))
* tuple to Tuple ([e567810](https://github.com/Loop3D/LoopStructural/commit/e567810e4fafd36da7ccf5696dd9245a904d4462))
* updating type hints ([a064224](https://github.com/Loop3D/LoopStructural/commit/a0642243fac0bd7e90f28957b95d68e31bac0af7))

## [1.6.1](https://github.com/Loop3D/LoopStructural/compare/v1.6.0...v1.6.1) (2024-08-05)


### Bug Fixes

* adding inequality pairs ([ce33ac9](https://github.com/Loop3D/LoopStructural/commit/ce33ac9914d04550192fe621070e3b1b19e7038b))
* adding loopsolver optional depencency + admm solver option ([26edd3f](https://github.com/Loop3D/LoopStructural/commit/26edd3fa8ee7de429a7cd3682f030dfebf79f40b))
* adding offset to fault ([b75df73](https://github.com/Loop3D/LoopStructural/commit/b75df73410f2a2662ef90cfb0b15e2413acc7d00))
* adding omf export ([d03949e](https://github.com/Loop3D/LoopStructural/commit/d03949e109c2d2c3ed35c4b3dcd7893f18179494))
* linting ([0e75342](https://github.com/Loop3D/LoopStructural/commit/0e75342788624a691755a889db7293b39308eb3c))
* linting ([ad3bb55](https://github.com/Loop3D/LoopStructural/commit/ad3bb5540dd5a232fd014b8131ffdc4429b6d649))
* linting ([0d7a052](https://github.com/Loop3D/LoopStructural/commit/0d7a0522bc706f2018bdd642cb2f28db303fd057))
* rename properties_cell to cell_properties and properties_vertex to properties for the structured grid ([4368eb6](https://github.com/Loop3D/LoopStructural/commit/4368eb60f7e2f170782a91cc6003159fdd98875b))
* updating point/surface export and constructor ([6b74cfd](https://github.com/Loop3D/LoopStructural/commit/6b74cfd5d5edda5769653be9e3bc3d9dd6cb4415))

## [1.6.0](https://github.com/Loop3D/LoopStructural/compare/v1.5.13...v1.6.0) (2024-06-07)


### Features

* removing visualisation module. ([a51355d](https://github.com/Loop3D/LoopStructural/commit/a51355d688ebd682e0068f45a3eba6c155ee6bcb))


### Bug Fixes

* add a vtk grid from bb attribute ([4b50d91](https://github.com/Loop3D/LoopStructural/commit/4b50d9107aa559896e545292eccc72e45af8e819))
* add get data ([c5ec096](https://github.com/Loop3D/LoopStructural/commit/c5ec096b6d11ac293fc393e4a2ba7a331daae163))
* add get data to faultdisplacement feature ([43f05d8](https://github.com/Loop3D/LoopStructural/commit/43f05d8c9ee80364f193dca4c64d08f0a6ba56c4))
* add isosurfacing and scalar field method to base feature ([bbbfe0c](https://github.com/Loop3D/LoopStructural/commit/bbbfe0cca8bf80d634c61a9667ef47b0ada4578a))
* add name argument to p1 gradient orthogonal ([4f92ba4](https://github.com/Loop3D/LoopStructural/commit/4f92ba4acad233f180ed80a16870377648f5dcf1))
* add other plotters to visualisation ([eb89b10](https://github.com/Loop3D/LoopStructural/commit/eb89b104764f6891408b434667e5a5890d137337))
* add scale paramater to generated vector field ([45661cf](https://github.com/Loop3D/LoopStructural/commit/45661cffc4b7eb9114f52952bf4c650bdbdd8484))
* adding case when strat column hasn't been set ([53a49dd](https://github.com/Loop3D/LoopStructural/commit/53a49dda7af06d10648914e3ae86576939e1fb6f))
* adding cg as a solver and dict for solver params ([fba928b](https://github.com/Loop3D/LoopStructural/commit/fba928b2749c153a277b58385d4131dad7950ded))
* adding export methods ([ceeee02](https://github.com/Loop3D/LoopStructural/commit/ceeee0277e05e4db9f558d5fabbef33cf7120707))
* adding fold example back in ([86d3377](https://github.com/Loop3D/LoopStructural/commit/86d33777563dfbbd2408ef0808f42ad73fd5604c))
* adding geoh5 export for points and grid ([c8641a6](https://github.com/Loop3D/LoopStructural/commit/c8641a6ce1aa33340a41e6282ed3ad7d325e5979))
* adding get_data for structural frame ([90dbb94](https://github.com/Loop3D/LoopStructural/commit/90dbb94a9840961dcaacb2c19f0ae19a3aac26b2))
* adding get_data method ([66ac1c1](https://github.com/Loop3D/LoopStructural/commit/66ac1c1cb333d7fcf061605deb83361846231858))
* adding method to evaluate all stratigraphic feature gradients. ([07b6078](https://github.com/Loop3D/LoopStructural/commit/07b6078cf98dd6b0f405ce5f55ad709fa6aea957))
* adding name argumen to finite difference gradient orthogonal ([064ae65](https://github.com/Loop3D/LoopStructural/commit/064ae65905f75c1b6a51c46438d9469a1604304d))
* adding ones column for constraints when weights missing ([44cc8fb](https://github.com/Loop3D/LoopStructural/commit/44cc8fb02607b433ced1635ef6942536060290fa))
* adding option to not store vertices for 2d unstructured supports to save memory ([64bc744](https://github.com/Loop3D/LoopStructural/commit/64bc7448911cdde58d5e0dee02d203aae8d0dcd4))
* adding pyvista to docs build ([b3bdcd7](https://github.com/Loop3D/LoopStructural/commit/b3bdcd78d7cdfb137b98b39e730fa62c1a9341fd))
* adding structure grid data class ([5e5035f](https://github.com/Loop3D/LoopStructural/commit/5e5035f5315ded01c8089c5ca9b31dd133898aed))
* adding surface getters ([6c75d4b](https://github.com/Loop3D/LoopStructural/commit/6c75d4b8f07cb7e686ed36580c6b7a94576e5217))
* adding value point save. ([2093f4e](https://github.com/Loop3D/LoopStructural/commit/2093f4e83400ae5092707d8a04330322547d51e7))
* adding vector/value points ([382eb74](https://github.com/Loop3D/LoopStructural/commit/382eb74bb7e0165c8de8fb6d6dfcf74c6ff9bfcb))
* allow faults without trace to be built if centre exists ([436540b](https://github.com/Loop3D/LoopStructural/commit/436540b03f99b78344b34a9b8c4f83780625adae))
* atol for isclose ([7165790](https://github.com/Loop3D/LoopStructural/commit/7165790202459ebfbbd4b01b54bd03edfc5d5913))
* auto select support dimension ([aba41bd](https://github.com/Loop3D/LoopStructural/commit/aba41bd13865075437bb575a119696ab14cb14e4))
* calculate fault normal using rotation of strike ([8707d53](https://github.com/Loop3D/LoopStructural/commit/8707d5358910cf0325e6329bc65d8e31a3167f54))
* cast data to float to avoid deprecation warning ([1692d67](https://github.com/Loop3D/LoopStructural/commit/1692d679fee3e4ec289699afcc61655502a08ea7))
* change fault segement evaluate value to scaled displacement vector ([3b8f8c3](https://github.com/Loop3D/LoopStructural/commit/3b8f8c3ec6ab6f885db2e1ea241f66d213ed3167))
* change feature.scalar_field() from vtk type to structure grid ([edadfb3](https://github.com/Loop3D/LoopStructural/commit/edadfb3f63eceb98c7624310d65fe5c89075e4fd))
* Change to new thickness name in project file ([4e567b7](https://github.com/Loop3D/LoopStructural/commit/4e567b77b0831b025259c7cd5727adce4f35fba6))
* changing vtk to method instead of attribute to allow parameter overloading ([bf30047](https://github.com/Loop3D/LoopStructural/commit/bf3004787af26a36e904b1db1bf7c43f467057fe))
* check if strat col for evaluate model ([83f5e86](https://github.com/Loop3D/LoopStructural/commit/83f5e86bbfbffa4d41809d9f0fd1647d751dbf97))
* created structuredgrid data type to replace vtk regulargrid ([15fdb3c](https://github.com/Loop3D/LoopStructural/commit/15fdb3c4a5a4c75af454186339cbc02b9361e685))
* custom solver updates solved state ([6d7264c](https://github.com/Loop3D/LoopStructural/commit/6d7264c4083ba1f09b47f980e586f43be100ff41))
* disable type hint because of circular import ([c2df4f4](https://github.com/Loop3D/LoopStructural/commit/c2df4f4ca8f4d53d430601c9394366b05a031c95))
* don't add fault above an unconformity! ([a5fc543](https://github.com/Loop3D/LoopStructural/commit/a5fc54349d4031cccd055d969ee324e6990b773c))
* don't add unconformities to unconformities. ([ffa11f1](https://github.com/Loop3D/LoopStructural/commit/ffa11f118fb03c0d1b76594e08ae8950bce06308))
* enforce loopstructuralvisualisation version ([53bd28b](https://github.com/Loop3D/LoopStructural/commit/53bd28b754652dfb796ced70b58f454df7c6103e))
* feature gradient masked by unconformities ([ae5324a](https://github.com/Loop3D/LoopStructural/commit/ae5324a250379996a8725437e71db7ee2e969340))
* get data for intrusion feature ([19fdc40](https://github.com/Loop3D/LoopStructural/commit/19fdc403e18fd073fa883ea46f46655c1ccc194c))
* implementing model.save method ([dd54899](https://github.com/Loop3D/LoopStructural/commit/dd548999b3caca2d0ea7ce8dbbb3092a1eb01bd1))
* indexing error ([b3f67cf](https://github.com/Loop3D/LoopStructural/commit/b3f67cfac756a81370b4d6d3d834a6c0a8e4ad83))
* indexing for applying fault rotation ([99f48f2](https://github.com/Loop3D/LoopStructural/commit/99f48f23631f53f110eccbcd63a732fbfbdb27fc))
* interpolator map support map for all 2d supports ([547fea0](https://github.com/Loop3D/LoopStructural/commit/547fea024255d647e5f25a0aaf04e5e1d163e9d3))
* interpolator support is rescaled for fault displacement. ([d886e81](https://github.com/Loop3D/LoopStructural/commit/d886e81b756a9efd867db85358f48307b840e137))
* linting ([a99c178](https://github.com/Loop3D/LoopStructural/commit/a99c178387405faa09db52045a63bf7ddaa6005b))
* making allow bbox to do rescaling ([5150860](https://github.com/Loop3D/LoopStructural/commit/5150860427c4cb3681fb29ad3278075852be74b0))
* moving bbox test ([8357c79](https://github.com/Loop3D/LoopStructural/commit/8357c79fb14dc066ce24e3d4a1e6315bcdce5115))
* parameter names, and adding get_data method ([e770a45](https://github.com/Loop3D/LoopStructural/commit/e770a45786f18780258704c6392ddb6c98930dfc))
* pass vector tollerance through vtk method for vectorpoints ([8343bc6](https://github.com/Loop3D/LoopStructural/commit/8343bc6f7bffdea2660747ca48ddcba02ba3e9ec))
* put fault normal points on the trace ([2391b30](https://github.com/Loop3D/LoopStructural/commit/2391b30ef3a344cf8b8df398bece2cfc5d43732f))
* remove api for now, Isosurfacer in utils ([e33c9dc](https://github.com/Loop3D/LoopStructural/commit/e33c9dc1b2ce060ec49a2899f883eb92bea253a9))
* remove shuffle ([0581fb1](https://github.com/Loop3D/LoopStructural/commit/0581fb10c60a76ea5d37cb03193ecb2188680140))
* remove surfepy import ([5ba8f14](https://github.com/Loop3D/LoopStructural/commit/5ba8f1427d2761c190b34c252199c91c7a388ec0))
* removign vol weighting ([47ad944](https://github.com/Loop3D/LoopStructural/commit/47ad9448126a38c1e7bd1be403ee9f4f2121b43e))
* removing random shuffle from orthogonal constraints ([4667906](https://github.com/Loop3D/LoopStructural/commit/466790683a1dbe87f0e110a634ac69d17b5070b4))
* rename id to stratigraphy for model block ([640ac0b](https://github.com/Loop3D/LoopStructural/commit/640ac0b0235b35a41e8633a1a84acd50a0d31885))
* return isosurfaces as a list rather than dictionary ([baa083a](https://github.com/Loop3D/LoopStructural/commit/baa083a68215ac1b89835c59ede0df77d4ea605f))
* return no strat ids when strat column not defined ([e77d58c](https://github.com/Loop3D/LoopStructural/commit/e77d58cc11b97ca2f78512fbc29a7056d59bc81e))
* set support to be the same for all fault frame components ([b9c7500](https://github.com/Loop3D/LoopStructural/commit/b9c7500436a73be754def958cb112d11b85cf1bd))
* step_vector take into account dimensions=2 ([3d90210](https://github.com/Loop3D/LoopStructural/commit/3d90210cf5b926360c3f2e41339796e1aab3d410))
* step_vector take into account dimensions=2 ([e62e6ea](https://github.com/Loop3D/LoopStructural/commit/e62e6ea46078c053ec0a52c9d2d3c397cb3a42d9))
* step_vector use self.dimensions to set length ([72e0c32](https://github.com/Loop3D/LoopStructural/commit/72e0c32af91a4ae1386a096c88e55eab8e39ed0f))
* surface export for tsurf, obj, vtk, pickle ([5812d3b](https://github.com/Loop3D/LoopStructural/commit/5812d3bccbf889c0828abb5a9f5281da260e182d))
* update for project file changes ([0a64def](https://github.com/Loop3D/LoopStructural/commit/0a64defc19da74d3da3f31accd34af916a60c1a5))
* updating bbox for 2d case as well as 3d. ([31f58e0](https://github.com/Loop3D/LoopStructural/commit/31f58e057b919c0071c877d2ee7f86a8e3c290f4))
* updating bounding box dimensions to use size of origin array ([e56c868](https://github.com/Loop3D/LoopStructural/commit/e56c868b0d38ff6db1936c4ef57702e1701eacea))
* updating bounding box for exporters ([f7e2571](https://github.com/Loop3D/LoopStructural/commit/f7e25717bbcb6652f438fe1212beddab39ac215d))
* updating imports ([40673c0](https://github.com/Loop3D/LoopStructural/commit/40673c0fc7c05f6ab1b6e77cf9e6b09026432cfa))
* updating solve_system kwargs ([ccbacff](https://github.com/Loop3D/LoopStructural/commit/ccbacffc8bc64abbcd5c93fd8acc00bbe7132bc9))
* upgrading libtiff ([6484857](https://github.com/Loop3D/LoopStructural/commit/648485761cb3e9669ebf57ba18138c31021279d5))
* Use ThicknessMedian instead of ThicknessMean as Mean isn't populated ([5f489b9](https://github.com/Loop3D/LoopStructural/commit/5f489b99ce5ba196fd97c17a35619a5d8db9266d))
* weights for vector constraints is optional ([72e2f06](https://github.com/Loop3D/LoopStructural/commit/72e2f06daee43c2aa41ba9751c43f291057d2f1b))


### Documentation

* adding geoh5py into docs build ([3274299](https://github.com/Loop3D/LoopStructural/commit/3274299b0434afaf1bbfce5c8b20fe13d49128ce))
* adding/updating docs ([fa8c564](https://github.com/Loop3D/LoopStructural/commit/fa8c5647c7df3feca6094d79b8ef781ebca13cc4))
* model.save ([2fb945e](https://github.com/Loop3D/LoopStructural/commit/2fb945eac062df616fc8785f5b4c9ea78a9a9d54))
* removing api from docs ([9170b3e](https://github.com/Loop3D/LoopStructural/commit/9170b3ec6cfbd2dac3c45e358f4bba4a6690b579))
* removing export example ([b8a45da](https://github.com/Loop3D/LoopStructural/commit/b8a45dad28dae90050e29e68a9b3bb1d98d387b5))
* updating docs for new visualisation ([0f55b5d](https://github.com/Loop3D/LoopStructural/commit/0f55b5d2680988cbf11a91a2afc18e8ead0ae8f8))
* upgrading visualisation version ([47f10d9](https://github.com/Loop3D/LoopStructural/commit/47f10d9a188c5588e360341a12a62871ce13b615))

## [1.5.13](https://github.com/Loop3D/LoopStructural/compare/v1.5.12...v1.5.13) (2024-03-07)


### Bug Fixes

* issue with 2d supports ([79bf9c0](https://github.com/Loop3D/LoopStructural/commit/79bf9c0ad5cfcde551d00ae76e1462349420f316))
* pyproject.toml includes submodules ([55c0e0a](https://github.com/Loop3D/LoopStructural/commit/55c0e0abeca8fa320fc44571417b16a4b08742ee))
* remove element volume scaling. ([4350c0e](https://github.com/Loop3D/LoopStructural/commit/4350c0e821a19871af6a6d4f60e3a6abe32e0d7f))
* update version of pypi upload action ([#180](https://github.com/Loop3D/LoopStructural/issues/180)) ([852f062](https://github.com/Loop3D/LoopStructural/commit/852f06250baeac190dd2961d3c5a1fece7506abc))

## [1.5.12](https://github.com/Loop3D/LoopStructural/compare/v1.5.11...v1.5.12) (2024-03-05)


### Bug Fixes

* adding imports back into init file ([b1d2ba9](https://github.com/Loop3D/LoopStructural/commit/b1d2ba97d5e4c7878165dc930d683df25fa61c49))
* adding pyproject.toml ([474101c](https://github.com/Loop3D/LoopStructural/commit/474101c00c8be510dd8e4e57152ff7b0a95a4cb1))
* bounding box can be defined from max, origin or nsteps + step vector ([a8d364c](https://github.com/Loop3D/LoopStructural/commit/a8d364ca5735fbd5efbe9fc2f535a5c6698b4d08))
* change pyvista property to vtk ([8277304](https://github.com/Loop3D/LoopStructural/commit/8277304e932cc2ebe6581a8409a9d0a2375698bc))
* exposing interpolator api in module import ([59e4f7c](https://github.com/Loop3D/LoopStructural/commit/59e4f7c196c6b803d6a642ea7976e4c733b5581a))
* if interpolator results in 0, still set up to date ([6f18788](https://github.com/Loop3D/LoopStructural/commit/6f187888a69844a271c16125d75f1c19ce8a2aa3))
* make rotation work with array of axis and angle ([72f7744](https://github.com/Loop3D/LoopStructural/commit/72f77441760e254a4e14e0a2738b041b2633499f))
* making code compatible with linter ([716038e](https://github.com/Loop3D/LoopStructural/commit/716038e6c69031f6f6d374e10e0b8964b09cd523))
* mapping interpolators to cythonless version ([e755c29](https://github.com/Loop3D/LoopStructural/commit/e755c29f851b5c67f74b8bd1dd26f278db349e18))
* moving interpolator api ([2e1b008](https://github.com/Loop3D/LoopStructural/commit/2e1b0087e12696b359ea8cc8d9643d8b39afff60))
* pyproject.toml missing modules ([4a4bf92](https://github.com/Loop3D/LoopStructural/commit/4a4bf92452fdb32ee3cf01f2155e11a17e832da5))
* Removing cython  dependency ([#168](https://github.com/Loop3D/LoopStructural/issues/168)) ([21e5572](https://github.com/Loop3D/LoopStructural/commit/21e5572a02329b8a8a18a328bd76a8cfb449344b))
* rename grad stepness to remove _ ([e8a7877](https://github.com/Loop3D/LoopStructural/commit/e8a78770d1ce6cc504a21f946363b739e52cff09))
* return axis and angle for movement of fault points ([48e6261](https://github.com/Loop3D/LoopStructural/commit/48e62616581cda2da89e0d79d323f8332d19a376))
* scaling fault normal by minor axis ([ac99448](https://github.com/Loop3D/LoopStructural/commit/ac994480b559d1f58c9c8befa54dc7016b4e2cb9))
* under constrained faults now work ([50a04af](https://github.com/Loop3D/LoopStructural/commit/50a04afee891be9f6b80108d9baf6fab71c8ffba))
* unstructured tetra only take first 4 elements ([4673bbc](https://github.com/Loop3D/LoopStructural/commit/4673bbc14d2517bf10138a741e72db910c42d2fe))


### Documentation

* updating design and contributors guide ([e10c7bc](https://github.com/Loop3D/LoopStructural/commit/e10c7bcb9872a6776a4270dd643b70c14fe49008))

## [1.5.11](https://github.com/Loop3D/LoopStructural/compare/v1.5.10...v1.5.11) (2023-12-03)


### Bug Fixes

* :art: moving get interpolator out of geologicalmodel ([3f52950](https://github.com/Loop3D/LoopStructural/commit/3f5295050bff2cbd0ecb769e703678550abcc632))
* :bug: add default solver to discrete interpolator ([3b93f7f](https://github.com/Loop3D/LoopStructural/commit/3b93f7f594ba519368ebde8d3ed6e3a83979ce14))
* :bug: evalute changed to evaluate ([4e93caf](https://github.com/Loop3D/LoopStructural/commit/4e93caf37f6c917f530bb0d26ffafdbff1279cb3))
* :bug: return a as array not transpose ([46e9edf](https://github.com/Loop3D/LoopStructural/commit/46e9edfa69442747bcdcd83af283a7584a8bfbf8))
* :sparkles: adding interpolator factory ([245517b](https://github.com/Loop3D/LoopStructural/commit/245517bdf3cfd2d0dce2f5b6961d54193e5f7a1b))
* :sparkles: adding support factory ([c3ec7a8](https://github.com/Loop3D/LoopStructural/commit/c3ec7a85dadb4e19bb65be6fd8d21852f54ce16c))
* :sparkles: adding support type enum ([82c0bd5](https://github.com/Loop3D/LoopStructural/commit/82c0bd59a16f7ccc0843141dc90fd14def19c1f4))
* :sparkles: bounding box object ([75941f0](https://github.com/Loop3D/LoopStructural/commit/75941f02d6cdfeaf29306ff6d55428cb77933dde))
* :sparkles: new api for accessing interpolation for a single scalar field ([5bca235](https://github.com/Loop3D/LoopStructural/commit/5bca2350e88f8e0f857e1b30da823de4a490f974))
* :sparkles: transformation object ([48e6c32](https://github.com/Loop3D/LoopStructural/commit/48e6c32b0149227b6c9efa80cf3a5cd6ef6c7b05))
* add logger to surface ([8fcf89c](https://github.com/Loop3D/LoopStructural/commit/8fcf89c79a51d26382ff160cbb43dd4c95b9cf27))
* adding a factory method to create an interpolator ([8aba4bd](https://github.com/Loop3D/LoopStructural/commit/8aba4bdfdecd1ef5dfd976c9c8145d7a171f6c3a))
* adding api module ([7ab7c07](https://github.com/Loop3D/LoopStructural/commit/7ab7c077ea983a8569613ab1fc50555bd4276421))
* adding datatypes module ([2102c72](https://github.com/Loop3D/LoopStructural/commit/2102c72c3f8a83aec416ad6c241728b16ac95fa6))
* adding example execution time to git ignore ([7ac99d3](https://github.com/Loop3D/LoopStructural/commit/7ac99d378b236f2fa4a48295cb0108852b9c7029))
* adding exception if nsteps=0 or <0 ([0940c55](https://github.com/Loop3D/LoopStructural/commit/0940c552b5c1ba931c8d6cfe02702c2e07d8020e))
* adding fixture for bb and interpolatortype ([2a3e84e](https://github.com/Loop3D/LoopStructural/commit/2a3e84ee88b12cccf2e28a8bd7ed7e9a585bba3a))
* adding isovalue to surfaces created by isosurfacer ([754541f](https://github.com/Loop3D/LoopStructural/commit/754541fc4010bc3632a101e4ca4f895325baa49b))
* adding mesa ([7f0fb94](https://github.com/Loop3D/LoopStructural/commit/7f0fb9476348fdf7dcf8402420793737a3da2dd5))
* adding nsteps to bb ([82e4dac](https://github.com/Loop3D/LoopStructural/commit/82e4dac7685882323ed65bfb47cb8fa4b4d7f69b))
* adding placeholder for interpolate api ([f9709b5](https://github.com/Loop3D/LoopStructural/commit/f9709b54515b97cee9eed429110cbb5101edd24c))
* adding surface data type for storing triangular surfaces ([0d383cf](https://github.com/Loop3D/LoopStructural/commit/0d383cf5d472852b5d2690b71efa9682b0b84fe4))
* adding surfacer ([47c390d](https://github.com/Loop3D/LoopStructural/commit/47c390d0abcb4fc825e48c1c344f6309212ed785))
* adding test to bb ([485f084](https://github.com/Loop3D/LoopStructural/commit/485f08485a0d602fc7479378faaf92ad698da255))
* bug with size of constraints array ([66a9825](https://github.com/Loop3D/LoopStructural/commit/66a982554358b6a629aa8d7c2a68469626bd30da))
* change to interpolation factory ([91bcbed](https://github.com/Loop3D/LoopStructural/commit/91bcbed956d2ca80f59c03c730ea53ba50c1436b))
* Check for empty features and escape early ([083a195](https://github.com/Loop3D/LoopStructural/commit/083a195004f45cc52086b503441f93f97cb3f21d))
* create folded fold frame using updated code ([7e1db31](https://github.com/Loop3D/LoopStructural/commit/7e1db31ed275099cd830b695dd416ff7aa7cf994))
* Ensure modifications to data frame are on a copy of that frame ([ed61243](https://github.com/Loop3D/LoopStructural/commit/ed6124380bc95ba91f55c88916c40e2a29463747))
* fixing lavavu version ([57f649b](https://github.com/Loop3D/LoopStructural/commit/57f649b252b9f2ca37dca42d82f97724a26be271))
* fixing python 3.10 for docs ([6ea4d84](https://github.com/Loop3D/LoopStructural/commit/6ea4d847de2b737d9ac9279fc5861aa848bfb0e5))
* flake8 ([9ef0081](https://github.com/Loop3D/LoopStructural/commit/9ef00811979f41eac624508eb436112d3060f43c))
* flake8 error ([a37fec2](https://github.com/Loop3D/LoopStructural/commit/a37fec28c1946fbd99bfe71eb8dd01dc0cc1966b))
* formatting ([8175437](https://github.com/Loop3D/LoopStructural/commit/817543773f649c5a03143029fcb6f81acf1fe56e))
* init file for datatypes ([7f1f62a](https://github.com/Loop3D/LoopStructural/commit/7f1f62ab5e6cd62077cbf594a45cee724feb7716))
* isinside and init with list ([0ff6735](https://github.com/Loop3D/LoopStructural/commit/0ff6735025ca6eaefc8494861a4b32f5f6af7463))
* manually incrementing version ([4229d6a](https://github.com/Loop3D/LoopStructural/commit/4229d6a95691e66b0b0fe63df3eb469140b22189))
* move create interpolator to feature builders ([1c01bdd](https://github.com/Loop3D/LoopStructural/commit/1c01bdde350e23e4018c1a61135d16ed309f766c))
* moving get interpolator to separate function ([2790d76](https://github.com/Loop3D/LoopStructural/commit/2790d761ec27f00be9a6beee7305a0827f4617b3))
* osmesa ([54595cb](https://github.com/Loop3D/LoopStructural/commit/54595cb5ebcc80de9bd09bb8bdd694685feb1e39))
* pypi only on release ([619557b](https://github.com/Loop3D/LoopStructural/commit/619557b3ff76c9dd5996fc6f8136de82369ca248))
* relocating bounding box ([8b18eb2](https://github.com/Loop3D/LoopStructural/commit/8b18eb20fbdea8f84bd87e9832c7f63cc2761613))
* removing calculate pairs from dsi ([e35003c](https://github.com/Loop3D/LoopStructural/commit/e35003cab2ca27fdc504557696f8e43cd3b9e2ed))
* removing print property values from lavavu ([306793c](https://github.com/Loop3D/LoopStructural/commit/306793ca4a997752a1e479f9c7ea50ac99695189))
* removing wavelength variable from refolded fold ([64278b7](https://github.com/Loop3D/LoopStructural/commit/64278b788bb74cad0e2dbf804f097e573a473f9d))
* return surface object from isosurfacer ([6f28d67](https://github.com/Loop3D/LoopStructural/commit/6f28d67cf11d3e6c14068ee0a6d7426ad9efb944))
* set min number of elements to 3 ([7a7e9ba](https://github.com/Loop3D/LoopStructural/commit/7a7e9ba62825cc03f5db176d009258c42e097c39))
* specify argument names for get_interpolator ([e203b4d](https://github.com/Loop3D/LoopStructural/commit/e203b4d5cfacbd6786b4aeda42ed44f45d10a664))
* trying lavavu without osmesa ([f41b7e7](https://github.com/Loop3D/LoopStructural/commit/f41b7e7caddcbc5936b898e1fe426d3fcbf8b1bd))
* trying to use tini for x ([4b38221](https://github.com/Loop3D/LoopStructural/commit/4b38221a6a09c4c4d9746dac379d3e9f23c1260c))
* updating dockerfile for documentaiton build ([1bd6d38](https://github.com/Loop3D/LoopStructural/commit/1bd6d387443280bb43de46a6c2e6bcde3c8785b7))
* updating docs to use LS version ([f8aac3a](https://github.com/Loop3D/LoopStructural/commit/f8aac3ab577b89f3b29de104df1a7c76840575cc))
* updating examples to not use pyamg ([66239cb](https://github.com/Loop3D/LoopStructural/commit/66239cb480d88c7f4d3a6530e841ae3524a65a8a))
* updating lavavu to use new bounding box ([96306db](https://github.com/Loop3D/LoopStructural/commit/96306db7ef7c41017183f67046ce32ac4f3ab5bb))
* updating product to prod ([5d050ca](https://github.com/Loop3D/LoopStructural/commit/5d050cab8c63591bf41291fe6c973a37d3849bee))
* updating so change work for building model ([08adb56](https://github.com/Loop3D/LoopStructural/commit/08adb5624e9d42a889c1def8fb1c0e126c75cffd))
* updating to new bbox ([fb450e0](https://github.com/Loop3D/LoopStructural/commit/fb450e02ef52b4959cf4e8e94af157bc9516a93a))
* Use fault names if present before labelling fault Fault_XXX ([7bf34fc](https://github.com/Loop3D/LoopStructural/commit/7bf34fc3e841e243b1a34dbaf36294b5ce17306f))
* use property array rather than store properties on grid ([2d95b85](https://github.com/Loop3D/LoopStructural/commit/2d95b852eb0f327daddbf6a786f354853c645242))
* value changed to point ([4144c59](https://github.com/Loop3D/LoopStructural/commit/4144c59d134b420875bb4247a361e26be1547c9c))

## [1.5.10](https://github.com/Loop3D/LoopStructural/compare/v1.5.9...v1.5.10) (2023-03-14)


### Bug Fixes

* adding 3D fault displacement function ([438e699](https://github.com/Loop3D/LoopStructural/commit/438e699842d4c97879fb5ab9a89fd2746845d41e))
* major change for base grid ([0a49817](https://github.com/Loop3D/LoopStructural/commit/0a4981785abf671a0153cce2c4c64c1a8650602c))
* trying to use bash to recognise anaconda command ([1ce8e32](https://github.com/Loop3D/LoopStructural/commit/1ce8e320d1fa03b5c8bda8a80a411a7361e009c6))
* updating tetmesh class for new indexing ([e53369a](https://github.com/Loop3D/LoopStructural/commit/e53369a42a1b4b990ed5c1a88bcf7440163200e7))

## [1.5.9](https://github.com/Loop3D/LoopStructural/compare/v1.5.8...v1.5.9) (2023-03-06)


### Bug Fixes

* adding clean function to interpolator ([d0de548](https://github.com/Loop3D/LoopStructural/commit/d0de548f6cbc756c04b029d294c991b1fc6004e5))
* adding geopandas to dockerfile for docs ([32637ec](https://github.com/Loop3D/LoopStructural/commit/32637ec5dea4f20334e56387e39c8d81005403c2))
* adding type hints ([916972c](https://github.com/Loop3D/LoopStructural/commit/916972c2f8da2169bca6f0a71a53e73fcfcbe8d6))
* all tests passing ([02a4d70](https://github.com/Loop3D/LoopStructural/commit/02a4d7015dcefd89caf6ed8852def42375fcfbb8))
* changing build script for docs ([c0c5ba2](https://github.com/Loop3D/LoopStructural/commit/c0c5ba203389ece188648d1a3c7d57e5b12948fc))
* check interpolator type for installing grad constraints ([bfa90ac](https://github.com/Loop3D/LoopStructural/commit/bfa90ac084e5880bd1e0a1ed1bcf5b7eb2b51c2f))
* code changes for merge ([7aeec5c](https://github.com/Loop3D/LoopStructural/commit/7aeec5c0166b743d6fef44b8b56ebdd382adc2da))
* code clenaing and lateral data update ([d202656](https://github.com/Loop3D/LoopStructural/commit/d20265673eb0efd1783434dea590fd8cea95eea9))
* dockerfile python version ([49de51e](https://github.com/Loop3D/LoopStructural/commit/49de51e86ecf699d6de6aa68e579920241dd0942))
* docs/docker building in python 3.10 ([26bfaa2](https://github.com/Loop3D/LoopStructural/commit/26bfaa2634f2c7567be9c6ce0bd67ea708238b6e))
* fault overprinting not added to model ([750ccce](https://github.com/Loop3D/LoopStructural/commit/750ccce0c53c01cce5b994750b6402ad8ad828e6))
* if w == 0 don't add constraint ([3aea341](https://github.com/Loop3D/LoopStructural/commit/3aea3410dfa551ef64e9772dc30fcba86c8c30dd))
* intrusion builder uses basebuilder ([bc0004a](https://github.com/Loop3D/LoopStructural/commit/bc0004a807e61f095f30141e64b435ca559bcff3))
* intrusions test updated ([e13a9aa](https://github.com/Loop3D/LoopStructural/commit/e13a9aaaa9bee5a7fda4800a8701906bb86a566a))
* kwargs in interpolator ([06f39d7](https://github.com/Loop3D/LoopStructural/commit/06f39d732807366469f83f0c0be1a6b07d843f4f))
* lateral constraint update ([c830149](https://github.com/Loop3D/LoopStructural/commit/c830149c9da44959cb5966b20bdf7a3a3b894694))
* remove of old SGS variables ([ab87b08](https://github.com/Loop3D/LoopStructural/commit/ab87b08040f1efb502c5290d1ef5acd407b8c9fb))
* remove temp variable ([d95394e](https://github.com/Loop3D/LoopStructural/commit/d95394ec79cc1af762d77a8b9c4dcfb4a03faac3))
* removing global_indices ([4960168](https://github.com/Loop3D/LoopStructural/commit/49601680479de64e67a402fe54b70078e56e3dd1))
* removing of #codes ([32c2b7e](https://github.com/Loop3D/LoopStructural/commit/32c2b7eba62ebc47381a162a8593545641c2c58f))
* removing unused imports ([11085fe](https://github.com/Loop3D/LoopStructural/commit/11085fe988754db68aaaaccc5af5ee5007783776))
* removing unused packages ([de08249](https://github.com/Loop3D/LoopStructural/commit/de082497756625711c0f3fb0e80c02b70f70343f))
* slicing for weights ([75890ca](https://github.com/Loop3D/LoopStructural/commit/75890ca4ce4cf258888330776453e2528d654417))
* svariogram parameter not being passed ([18d9fe2](https://github.com/Loop3D/LoopStructural/commit/18d9fe2f6ef27e896e4680eeac2e5ef8dd6027be))
* type hints and formatting ([2b9dd41](https://github.com/Loop3D/LoopStructural/commit/2b9dd4134cafda56ccb3ae8b598777a7a9b6f529))
* update of conceptual model functions ([b61e9d3](https://github.com/Loop3D/LoopStructural/commit/b61e9d356c0f226f9eb1e9a9047170d08aa2383c))
* update of intrusion interpolation parameters ([0140cd4](https://github.com/Loop3D/LoopStructural/commit/0140cd47ea4b19d21f150a606a608f1dfa87e807))
* updating ci action versions ([478a401](https://github.com/Loop3D/LoopStructural/commit/478a401925f1af919841080113226dc011e2be77))


### Documentation

* making dockerfile build ([0c91968](https://github.com/Loop3D/LoopStructural/commit/0c919684fc80239b5f2e0fdc7e901145f4ebfb88))

## [1.5.8](https://github.com/Loop3D/LoopStructural/compare/v1.5.7...v1.5.8) (2023-01-24)


### Bug Fixes

* force version bump again ([a38f60d](https://github.com/Loop3D/LoopStructural/commit/a38f60d79a95c12e5a4f45cb9359569e710781cf))

## [1.5.7](https://github.com/Loop3D/LoopStructural/compare/v1.5.6...v1.5.7) (2023-01-23)


### Bug Fixes

* triggering release ([c63f12f](https://github.com/Loop3D/LoopStructural/commit/c63f12f09d20e8b5cbf3f1a317486836ee6d057f))

## [1.5.6](https://github.com/Loop3D/LoopStructural/compare/v1.5.5...v1.5.6) (2023-01-19)


### Bug Fixes

* force bump version ([b5cb2ad](https://github.com/Loop3D/LoopStructural/commit/b5cb2ad8ff43696c844eb3f5dd9f20930335a41d))

## [1.5.5](https://github.com/Loop3D/LoopStructural/compare/v1.5.4...v1.5.5) (2023-01-17)


### Bug Fixes

* adding numpy and cython to sdist ([c6ff2e1](https://github.com/Loop3D/LoopStructural/commit/c6ff2e11794ebe45f4ce393d67b309ffae1d4868))
* ci formatting ([44d434f](https://github.com/Loop3D/LoopStructural/commit/44d434f68706cdefa359194a189f03c8db52f280))
* fix numpy to 1.21 ([13e0796](https://github.com/Loop3D/LoopStructural/commit/13e0796a64d22b485168ba6fd198b68e0e70c6ec))

## [1.5.4](https://github.com/Loop3D/LoopStructural/compare/v1.5.3...v1.5.4) (2023-01-12)


### Bug Fixes

* absolute import to relative import ([0ff0ba9](https://github.com/Loop3D/LoopStructural/commit/0ff0ba95f5efa4a543232a901b8fcc1ad542a216))
* adding meshio to docs docker ([12ca5d3](https://github.com/Loop3D/LoopStructural/commit/12ca5d3e477e10cbf5737b112acba1c039358f80))
* disabling fault_network example ([c390b7c](https://github.com/Loop3D/LoopStructural/commit/c390b7cd2a823b8445c38fdc7561681647749d07))
* fixing bug with fault network example ([283cdfe](https://github.com/Loop3D/LoopStructural/commit/283cdfe6b01129fa741e0005168a0f10372a8e89))
* incorrect variable in process data ([1776269](https://github.com/Loop3D/LoopStructural/commit/177626955ed3cd9ec6221eed2b774afccc989930))
* lock python to 3.8 for lavavu ([f38b6a1](https://github.com/Loop3D/LoopStructural/commit/f38b6a12e31da188278d020470dd57cee6cd223f))
* making dockerfile use pypi code ([6462c41](https://github.com/Loop3D/LoopStructural/commit/6462c4101741f32389d1484f6f4c466ec2c576c8))
* removing old create dtm function ([b559674](https://github.com/Loop3D/LoopStructural/commit/b559674d46580706bc4136bd561db5f14dfc076f))
* removing probability module ([5b0b6a2](https://github.com/Loop3D/LoopStructural/commit/5b0b6a2cca303066439c6e9d634798fdfeb93535))
* typo in feature exporter ([a735e9f](https://github.com/Loop3D/LoopStructural/commit/a735e9ffcfe2d42dff06a28752194aa4b44829ec))


### Documentation

* fixing lavavu version ([a993f60](https://github.com/Loop3D/LoopStructural/commit/a993f6053d0814b36ac128ffceef7dffda9fb719))
* rename builder to builders ([d332399](https://github.com/Loop3D/LoopStructural/commit/d33239921623dbd909cc94bc31353e4b2ff11efd))

## [1.5.3](https://github.com/Loop3D/LoopStructural/compare/v1.5.2...v1.5.3) (2022-11-02)


### Bug Fixes

* adding copy method to geological feature ([01651cb](https://github.com/Loop3D/LoopStructural/commit/01651cb1a685462c375c9b379de8be345f6241f7))
* adding debugging page to doc index ([a09aa0d](https://github.com/Loop3D/LoopStructural/commit/a09aa0d276b797436a91a8d6197d171db7909b58))
* adding interpolation back to builder ([1a537aa](https://github.com/Loop3D/LoopStructural/commit/1a537aa6ed06264de0ad40aae28d4cc672cdab6c))
* allow regularisation to be set for faults using processor ([07920dc](https://github.com/Loop3D/LoopStructural/commit/07920dc076c0147e962da24378188abe28a2931e))
* model set for fault segment ([3eb0328](https://github.com/Loop3D/LoopStructural/commit/3eb03281c7e343cff781402bcc0828a7a4b63285))
* remove prints while solving ([3aae8a9](https://github.com/Loop3D/LoopStructural/commit/3aae8a9e6c2240dc2917bae10aa3e43b9194bbf0))
* removing theano from docs docker ([a706575](https://github.com/Loop3D/LoopStructural/commit/a706575fc8c9828bf97cb1a7e419f3b5ad833303))
* tolerance geological modle ([fbc7e4a](https://github.com/Loop3D/LoopStructural/commit/fbc7e4a95c9618b7f8bacea4929e01da4ce83c19))
* updating dockerfile ([70ea4f7](https://github.com/Loop3D/LoopStructural/commit/70ea4f72334ca1dcd31b83129a315cf9aabdba0f))
* updating dockerfile so that it builds ([0e5482a](https://github.com/Loop3D/LoopStructural/commit/0e5482a9e1b7b353827c2db53920efe5df33a080))


### Documentation

* adding more documentation on structural frames and faults ([2fab248](https://github.com/Loop3D/LoopStructural/commit/2fab2483b16207576052ffa1e7f809e1bf08b92e))
* path for fault frame image ([921cd19](https://github.com/Loop3D/LoopStructural/commit/921cd195ec492571ac6fe68abca51159ce7258df))

## [1.5.2](https://github.com/Loop3D/LoopStructural/compare/v1.5.1...v1.5.2) (2022-10-11)


### Bug Fixes

* adding BaseBuilder to builders init ([7f56fdd](https://github.com/Loop3D/LoopStructural/commit/7f56fdd72c8edeb188217f4314c02179c1ba790e))
* adding checks for data size in interpolator ([ea8ccca](https://github.com/Loop3D/LoopStructural/commit/ea8cccacfa7b91f92b08eb697dadd48740bd7c65))
* adding old auto_examples directory ([cefe790](https://github.com/Loop3D/LoopStructural/commit/cefe790e1e93bfe6c727800707fcb8ecdc0bdeaf))
* code review for intrusions ([debc391](https://github.com/Loop3D/LoopStructural/commit/debc391501d6ccb86769ab7109c80cf8a3c32ffa))
* creating basebuilder ([3480e3d](https://github.com/Loop3D/LoopStructural/commit/3480e3dec5c4260fdb96f42ea68499445be42d6b))
* feature builder sets interpolator interpolation__region ([520cc01](https://github.com/Loop3D/LoopStructural/commit/520cc013bc1ceb4cbace9afcfa224a6239c3e726))
* formatting ([9cb79cb](https://github.com/Loop3D/LoopStructural/commit/9cb79cb91327ccf777630dadfe782e2434475bc3))
* starting unit id @ 1 means basement =0 ([ba627bb](https://github.com/Loop3D/LoopStructural/commit/ba627bbe5d251efc6c2adf7036930aea477588a4))


### Documentation

* adding basic testing guideline ([7c6ccef](https://github.com/Loop3D/LoopStructural/commit/7c6ccefbf04fb293453a7f13b04b2371e956ece3))
* adding docstrings ([ffe9614](https://github.com/Loop3D/LoopStructural/commit/ffe9614bcfecf1606fb1130e976df7264987e597))
* adding geopandas to documentation docker ([ab443ab](https://github.com/Loop3D/LoopStructural/commit/ab443ab7c9dd8e14038782319ae82fa70896ae58))
* adding links to kwargs in docstring ([6e3a6e2](https://github.com/Loop3D/LoopStructural/commit/6e3a6e2b809842c13bf0700192a098d0f46f6c40))
* fixed bug with documentation not being generated ([994d5dc](https://github.com/Loop3D/LoopStructural/commit/994d5dc20d4e44840271b4bcf1df6e5c5325ffa1))

## [1.5.1](https://github.com/Loop3D/LoopStructural/compare/v1.5.0...v1.5.1) (2022-09-15)


### Bug Fixes

* bump ([d175b99](https://github.com/Loop3D/LoopStructural/commit/d175b99f925734b32d62b6b50d2ed2064ebfedc5))

## [1.5.0](https://github.com/Loop3D/LoopStructural/compare/v1.4.13...v1.5.0) (2022-09-14)


### Features

* bump version to test release-please ([08d380d](https://github.com/Loop3D/LoopStructural/commit/08d380d6f3d050f752848a566f62a6d4e849ffc1))


### Bug Fixes

* setting versionfile for version.py ([6fc564e](https://github.com/Loop3D/LoopStructural/commit/6fc564e38a996888bd71804213d0495a5a6714e0))


### Documentation

* adding example for local weight change ([9a8f773](https://github.com/Loop3D/LoopStructural/commit/9a8f7732774a4ff89f4e8c20846e5bb9d6ad71b1))

### [1.4.13](https://www.github.com/Loop3D/LoopStructural/compare/v1.4.12...v1.4.13) (2022-08-29)


### Bug Fixes

* removing example that crashes docs ([63a9ce9](https://www.github.com/Loop3D/LoopStructural/commit/63a9ce907339e8a12395b4f7b1bf150e5b537e18))

### [1.4.12](https://www.github.com/Loop3D/LoopStructural/compare/v1.4.11...v1.4.12) (2022-08-27)


### Bug Fixes

* enabling fault network example ([5560d3f](https://www.github.com/Loop3D/LoopStructural/commit/5560d3f03aee7213dec7088b0e07f8a4150935f4))
* failfast in strategy ([2f3304e](https://www.github.com/Loop3D/LoopStructural/commit/2f3304ee40220ac1b5f49ddf6b29d6a2e524520b))
* path to copy change log ([c061755](https://www.github.com/Loop3D/LoopStructural/commit/c0617555acb1f1f57247c1fa7138e9a412f06207))

### [1.4.11](https://www.github.com/Loop3D/LoopStructural/compare/v1.4.10...v1.4.11) (2022-08-10)


### Bug Fixes

* add data for structural frame, adds data for all coordinates ([6ebdf84](https://www.github.com/Loop3D/LoopStructural/commit/6ebdf842071d8822053008c6e7a19eb8c0131414))
* adding new dataset for fault, map in installer ([b32edf3](https://www.github.com/Loop3D/LoopStructural/commit/b32edf3a739b51347db9233d8b35a8f92fa5bbb5))
* basal unconformities cropping lowest surface ([5c47ef6](https://www.github.com/Loop3D/LoopStructural/commit/5c47ef6020b2d2c1da53e771db5dac1308d870db))
* faults not added to features because using string not enum ([faf71ee](https://www.github.com/Loop3D/LoopStructural/commit/faf71ee6584d3da63e77343847dfbcfafefb4b21))
* faults were not faulting unconformities ([ec520fa](https://www.github.com/Loop3D/LoopStructural/commit/ec520fa88f7b0c0efeb1c80516b6b77967cd7982))
* geological map example stratigraphic order loaded incorrectly ([b98397b](https://www.github.com/Loop3D/LoopStructural/commit/b98397b6adafc279f0bfa9701a77b8b79f779aed))
* single group surfaces not plotting, changed min strat column to 0 ([44481a3](https://www.github.com/Loop3D/LoopStructural/commit/44481a3155b0dafcbe009f148625362dbcda1770))

### [1.4.10](https://www.github.com/Loop3D/LoopStructural/compare/v1.4.9...v1.4.10) (2022-07-20)


### Bug Fixes

* :bug: unconformities weren't working. ([304335d](https://www.github.com/Loop3D/LoopStructural/commit/304335d0421297ca026fa26a2e6fed2e4c16332b))
* adding loopjson encoder ([ec1b84e](https://www.github.com/Loop3D/LoopStructural/commit/ec1b84e7881977d7ec1fec7ecf00b45770f5d832))
* adding more informative errors to p2 tetra ([c63cae7](https://www.github.com/Loop3D/LoopStructural/commit/c63cae7e9c03be10b566d4b5b85987956d164a78))
* adding python alternatives for cg and fold cg ([5ffef3b](https://www.github.com/Loop3D/LoopStructural/commit/5ffef3b482e1fedaaf2ec6a6c8c44a392d6def6f))
* catch all exceptions, raise import ([a071416](https://www.github.com/Loop3D/LoopStructural/commit/a071416cb602db135c9fa9d67d99ceb79ddb0407))
* catch dfi import error ([f9c8aa0](https://www.github.com/Loop3D/LoopStructural/commit/f9c8aa0b04996dd5ad887871b950b2851c4310bb))
* catch modulenotfound raise importerror ([836ee84](https://www.github.com/Loop3D/LoopStructural/commit/836ee843f418956fc9bd578868030eba537960e2))
* don't add unconformity if base feature is None ([be4d8ac](https://www.github.com/Loop3D/LoopStructural/commit/be4d8ac220d1fe9a58738652ccf823631cab607e))
* faults enabled bug, ~ does not flip boolean, changing to not ([0c68788](https://www.github.com/Loop3D/LoopStructural/commit/0c6878837e8eeb24b2ca9009cced0a8ecb6750ac))
* feature type enums renaming ([616b554](https://www.github.com/Loop3D/LoopStructural/commit/616b55440944055d284484381ace2d068e5259c1))
* featuretype enum ([c256746](https://www.github.com/Loop3D/LoopStructural/commit/c2567465855b9a8ae6a97c1ded5972ec18a89b30))
* fold cg ([a35f41a](https://www.github.com/Loop3D/LoopStructural/commit/a35f41a0a1f2be50e9d15ae8d4e2d92877526eff))
* lavavu.camera bug when keys are missing ([b554c3c](https://www.github.com/Loop3D/LoopStructural/commit/b554c3c143b27e592672f87d1982ae455779aade))
* missing enum ([bc34bb4](https://www.github.com/Loop3D/LoopStructural/commit/bc34bb4e8cb07a006529c26a7982b3e2b3ea0241))
* model plotter uses BaseFeature to check type ([a2bd0f0](https://www.github.com/Loop3D/LoopStructural/commit/a2bd0f01a7c297ecf01fc0c612dc0d123ac16702))
* pli import is optional ([6347742](https://www.github.com/Loop3D/LoopStructural/commit/6347742847460821a11ce840af51972f0b5df43b))
* raise error if cython can't be imported ([f5b82eb](https://www.github.com/Loop3D/LoopStructural/commit/f5b82ebf96e7b93052693f99af783ea18511aa48))
* recursive error fix ([f3180f3](https://www.github.com/Loop3D/LoopStructural/commit/f3180f3beaae36a424be1de7f2512800fbe39022))
* renaming enums to be consistent with calls ([57ed443](https://www.github.com/Loop3D/LoopStructural/commit/57ed44349c5b8f04e967d049f55aa35f1440e324))
* unconformities causing nan slicing due to recurvsive error ([728fc5f](https://www.github.com/Loop3D/LoopStructural/commit/728fc5f053345ae6901631bf60bd1d2457574d66))


### Documentation

* :memo: added example for plotting unconformities ([bad53a0](https://www.github.com/Loop3D/LoopStructural/commit/bad53a01e21e39dcdbdc3a5c218e9c9957050f9b))
* fixing new documentation so that data can be loaded ([7ac5a0d](https://www.github.com/Loop3D/LoopStructural/commit/7ac5a0d1752c6aad2e8f5112ec4d6efeaa93dd20))
* missing variable ([7df88a2](https://www.github.com/Loop3D/LoopStructural/commit/7df88a2111c259de367b734a095452a747b8a019))
* renaming cython to _cython so its not included in docs ([27e88a1](https://www.github.com/Loop3D/LoopStructural/commit/27e88a1d4f3c8f3461481972fc1b13fda7513f0d))
* trying to ignore cython ([6582d0c](https://www.github.com/Loop3D/LoopStructural/commit/6582d0c8381588e8f59c623641c7d8ce451891fc))
* updating documentation ([493beea](https://www.github.com/Loop3D/LoopStructural/commit/493beeac3ade0b35c49895d91c64aa8926d9a3b7))
* using docker to build docs ([f499aae](https://www.github.com/Loop3D/LoopStructural/commit/f499aae76107e09e15d2b9ce56e75fa9b2e0560b))

### [1.4.9](https://www.github.com/Loop3D/LoopStructural/compare/v1.4.8...v1.4.9) (2022-05-05)


### Documentation

* fixing docs ci ([2e0d472](https://www.github.com/Loop3D/LoopStructural/commit/2e0d47292ccc20c071256cbb51ac8736c6851f9e))

### [1.4.8](https://www.github.com/Loop3D/LoopStructural/compare/v1.4.7...v1.4.8) (2022-05-05)


### Bug Fixes

* :bug: faults where feature name given not fault_name were crashing ([dae9e92](https://www.github.com/Loop3D/LoopStructural/commit/dae9e92f4b2f94ebcaefe85dfc419bfdef69b43a))


### Documentation

* adding geopandas to doc test ([b6613a9](https://www.github.com/Loop3D/LoopStructural/commit/b6613a91fd158405ac01d8c9675618bcb6e6af12))
* fixing fault network example ([6f4f2fc](https://www.github.com/Loop3D/LoopStructural/commit/6f4f2fcba8e1e708c17acfa569d6592f0a33c9d6))
* refactoring documentation + adding new example ([d15ba6c](https://www.github.com/Loop3D/LoopStructural/commit/d15ba6c391af889e0a0e57f4eced126fe6521461))

### [1.4.7](https://www.github.com/Loop3D/LoopStructural/compare/v1.4.6...v1.4.7) (2022-05-04)


### Bug Fixes

* :bug: fault function attribute is now a property accessor ([3b9cc2e](https://www.github.com/Loop3D/LoopStructural/commit/3b9cc2e29d10d5b05af12da23ac8ccc6c2dca11a))
* :zap: flag not uptodate when build arguments change ([6f80029](https://www.github.com/Loop3D/LoopStructural/commit/6f8002911af8d52def836bc555193cfe20485d9c))
* changing import error message to warning ([de22530](https://www.github.com/Loop3D/LoopStructural/commit/de22530702a55ce322ae5cbfeeaa8c0ed6dd2875))
* fold frame bug for refolded folds ([bf272d1](https://www.github.com/Loop3D/LoopStructural/commit/bf272d19c55338421d24edda9255d184b5ce60e3))
* if downthrow not provided estimate abutting direction ([5a7c8d5](https://www.github.com/Loop3D/LoopStructural/commit/5a7c8d5ee16b1f8d080e496894507f1ef09ed102))
* temp disable fault centre from project file ([44d813d](https://www.github.com/Loop3D/LoopStructural/commit/44d813d4f4b66fffcbce97099b92e2c04ceb4294))


### Documentation

* :memo: updating docs for pydata, actually committing changes ([be33a6f](https://www.github.com/Loop3D/LoopStructural/commit/be33a6fe2d8bd1d8f69117203c7bc460706eda9c))
* :memo: updating documentation to use pydata theme ([c5b28d5](https://www.github.com/Loop3D/LoopStructural/commit/c5b28d5238663d581c661bc4220487253cba4eca))
* adding fault network example ([2bcbfac](https://www.github.com/Loop3D/LoopStructural/commit/2bcbfac75cc06312a7bb55e9abd5c06a01abbaeb))
* adding more to design document ([8322c60](https://www.github.com/Loop3D/LoopStructural/commit/8322c6036f95cb5e382e54176e6d73200d55bf27))

### [1.4.6](https://www.github.com/Loop3D/LoopStructural/compare/v1.4.5...v1.4.6) (2022-04-26)


### Bug Fixes

* :ambulance: quick fix for projectfile import error ([34fb1d6](https://www.github.com/Loop3D/LoopStructural/commit/34fb1d695b101cfbefbcb5f8e8bbb0449e170c18))
* :ambulance: try 2, catching loopimporterror as well as import error ([5db6a81](https://www.github.com/Loop3D/LoopStructural/commit/5db6a8113d5dd0f6f2c3c59078cbd812f4aecafb))

### [1.4.5](https://www.github.com/Loop3D/LoopStructural/compare/v1.4.4...v1.4.5) (2022-04-14)


### Bug Fixes

* :art: adding LoopProjectFileProcessor to modeling init ([fa1c8f7](https://www.github.com/Loop3D/LoopStructural/commit/fa1c8f7333d4c387988ddf3969e51958df5d123a))
* :bug: adding origin and maximum to the data processor ([bb7b811](https://www.github.com/Loop3D/LoopStructural/commit/bb7b811bf773ffbeeafb270188ed1ee4992afa01))
* :bug: adding scaling to splay faults to reduce bubbles ([7d618fe](https://www.github.com/Loop3D/LoopStructural/commit/7d618fedfa6ca9adb8f3c2040cc2ac93630ce2cd))
* :bug: can add faults even when stratigraphic column hasn't been set ([33366f9](https://www.github.com/Loop3D/LoopStructural/commit/33366f926924f73ba8139c9d87e59c7ba09d06a2))
* :bug: fault geometry weren't being set if they were calculated from the fault trace ([ddec5cb](https://www.github.com/Loop3D/LoopStructural/commit/ddec5cb45259186aef4bde2bbf98e0bd3f518079))
* :bug: fixing warning message with divide by 0 ([ebd12ff](https://www.github.com/Loop3D/LoopStructural/commit/ebd12ff6a0af0a312076327a2ab26a03a5333661))
* bug with fold constraints causing crash ([3abb65e](https://www.github.com/Loop3D/LoopStructural/commit/3abb65e47d2f90f9f2e2500ddbe0d362bec136fb))
* bug with fold interpolator ([8d94ca4](https://www.github.com/Loop3D/LoopStructural/commit/8d94ca4cd8aa2ff7de17cae6b53f77a3abc06af4))
* temp disable splay fault ([721270e](https://www.github.com/Loop3D/LoopStructural/commit/721270e942258c9c220aaaaaf01391e36d9af958))

### [1.4.4](https://www.github.com/Loop3D/LoopStructural/compare/v1.4.3...v1.4.4) (2022-03-21)


### Bug Fixes

* :ambulance: import failing when surfepy not install ([5041ad6](https://www.github.com/Loop3D/LoopStructural/commit/5041ad6692cbcd5f3a58b9edb2fc8168e73905ec))
* :ambulance: number of steps was actually number of cells causing issues with indexing ([88b8fee](https://www.github.com/Loop3D/LoopStructural/commit/88b8fee380cf50a35801247a3c430d84eb4e14f4))
* :ambulance: surfe import was failing ([63d61a4](https://www.github.com/Loop3D/LoopStructural/commit/63d61a4253a9b10cbc7e4e07ae73d42fbd1aa84e))
* :bug: fixed bug with fault clipping, uses normal direction to infer hanging wall side away from interpolated surfaces ([7ea9beb](https://www.github.com/Loop3D/LoopStructural/commit/7ea9beb33993a82f87686c530f7f4d4851a6b8f9))
* :bug: small fixes for making sure model update is called when params are changed ([a22ed61](https://www.github.com/Loop3D/LoopStructural/commit/a22ed6112e5cde81577523fdb142bfee446128ad))
* :bug: stratigraphic column was not covering scalar field range. ([990e0f4](https://www.github.com/Loop3D/LoopStructural/commit/990e0f48cd02af1965418f3371d2e267e6affdb8))
* :bug: weight for p2 interpolator wasn't being applied ([5002a8a](https://www.github.com/Loop3D/LoopStructural/commit/5002a8a69cdf1d69d52a36f877b81a49c8e86fa1))
* :fire: changes fault slip warning to info ([16d16d7](https://www.github.com/Loop3D/LoopStructural/commit/16d16d7442f2f2a96ca457ec82ee3ebe5b7e1929))
* :rotating_light: fixig pylint error ([d5c99cd](https://www.github.com/Loop3D/LoopStructural/commit/d5c99cdaffc94e4d6883de395923ae0278f23fe0))
* :sparkles: adding config class for storing module wide configs ([9407ff7](https://www.github.com/Loop3D/LoopStructural/commit/9407ff7b523160786430e390491a49b62dd631d0))
* :sparkles: faults/fault builder now store the fault geometry parameters ([f92749b](https://www.github.com/Loop3D/LoopStructural/commit/f92749b22494bb418f219a19f4f880f6278ee946))
* :sparkles: flagging interpolator not up to date when data changes ([96c8ddb](https://www.github.com/Loop3D/LoopStructural/commit/96c8ddbe29d7f6e61efac89d6d0d9fca0d750089))
* :sparkles: intrusionbuilder-->intrusionframebuilder, intrusionbuilder will build intrusionfeature ([7d5935c](https://www.github.com/Loop3D/LoopStructural/commit/7d5935c53c751096c1ca019dd0df31a71991f767))
* adding debugging code for indexing structured grid ([59e130a](https://www.github.com/Loop3D/LoopStructural/commit/59e130a422db5a547d8fca4a228811074d7b5af6))
* adding intrusion dependencies to setup.py ([08f6a22](https://www.github.com/Loop3D/LoopStructural/commit/08f6a224d2479814fc00cedb70013299f0617bdb))
* adding logger and remove unused files ([bd3f428](https://www.github.com/Loop3D/LoopStructural/commit/bd3f428d65c0ce369b4b36970e7abce11a68abeb))
* change logs ([b8000db](https://www.github.com/Loop3D/LoopStructural/commit/b8000db5caa6638360d4e1d9ace5ca59f8e3397a))
* change pandas append to concat ([51bbf5e](https://www.github.com/Loop3D/LoopStructural/commit/51bbf5e69fc62f902d2be7bb6706edab879c4401))
* changing ckmeans to sklearn kmeans ([60f3ec5](https://www.github.com/Loop3D/LoopStructural/commit/60f3ec5d6761e4bde94ea8dd1a56e16afa5506d6))
* fixing import warnings ([3b5c0f3](https://www.github.com/Loop3D/LoopStructural/commit/3b5c0f3a3b48d39cc5f7ddac664aec99b3cf91c0))
* improvements in modified scalar field ([d1de75d](https://www.github.com/Loop3D/LoopStructural/commit/d1de75d01b8131f150fc936f1ff4d74b71ef588d))
* Intrusion Builder class ([f697ee6](https://www.github.com/Loop3D/LoopStructural/commit/f697ee62e4c948c74083f5762bb0b9333a8d8086))
* intrusion feature containing pre IBody fxs ([cd6fe61](https://www.github.com/Loop3D/LoopStructural/commit/cd6fe6103505fb7f36eeacf91fcead01de17b951))
* merge origin/intrusion to local intrusion ([2507b33](https://www.github.com/Loop3D/LoopStructural/commit/2507b33b5a45267eddd0a7c944a5250a3d457db1))
* moving config to own file ([3dc2a19](https://www.github.com/Loop3D/LoopStructural/commit/3dc2a199bbe8cae904255c575c32c718004580c1))
* remove of ckwrap library ([e5ddf5e](https://www.github.com/Loop3D/LoopStructural/commit/e5ddf5ee9f96e995809cfb0426c2a78beb506c9f))
* temp removing rotate from structured support ([2e2d1c6](https://www.github.com/Loop3D/LoopStructural/commit/2e2d1c6774eb597f099488d3f83f0af4dddc7702))
* unit tests for intrusions ([85b9b4b](https://www.github.com/Loop3D/LoopStructural/commit/85b9b4b034ffb13a811bc30ea4ccf7487b7e257e))
* update of builder ([cdf7ebb](https://www.github.com/Loop3D/LoopStructural/commit/cdf7ebb5e70ea238a426d6161928d86cd818b38a))
* update of conceptual models ([948beef](https://www.github.com/Loop3D/LoopStructural/commit/948beef1cf8cadf516ec20610fb89d878bddcc0f))
* update of intrusion module ([bba0163](https://www.github.com/Loop3D/LoopStructural/commit/bba0163fc4d93def34a5fbb38d948e9c13280802))
* update skfmm import ([ea855f0](https://www.github.com/Loop3D/LoopStructural/commit/ea855f03902e093855f3c1046ceabb9243b5dd90))
* variogram parameters now in set_sgs_parameters ([26ca7b3](https://www.github.com/Loop3D/LoopStructural/commit/26ca7b37ac298ab0c99446820324e56e758f8223))


### Performance Improvements

* :zap: reducing memory for unstructured tetra ([1ee5dd6](https://www.github.com/Loop3D/LoopStructural/commit/1ee5dd6934602aff82c280429e62f18c699cef2e))

### [1.4.3](https://www.github.com/Loop3D/LoopStructural/compare/v1.4.2...v1.4.3) (2022-02-21)


### Bug Fixes

* :ambulance: added scale parameter for calculate_topo matrix. ([839e0e0](https://www.github.com/Loop3D/LoopStructural/commit/839e0e07fd8333c70ca71c0d8cf356ccd2e9e604))
* :bug: adding experimental flag for loopstructural ([f810479](https://www.github.com/Loop3D/LoopStructural/commit/f810479f6d057d21f9229ae645d39257b051d906))
* :bug: changes support barycentre to property ([d2e1f0b](https://www.github.com/Loop3D/LoopStructural/commit/d2e1f0b54d71eb2cdae3a92e9fcbce9a7e189d26))
* :bug: fixing bugs for flake8 ([9dce30c](https://www.github.com/Loop3D/LoopStructural/commit/9dce30c012ba6543a7936275a4452c25d30412fc))
* :bug: small changes to pass tests. ([14c7580](https://www.github.com/Loop3D/LoopStructural/commit/14c758044d36962b6b63d7af7dd6783638539023))
* :zap: changing aabb to sparse matrix ([2cbdab7](https://www.github.com/Loop3D/LoopStructural/commit/2cbdab7feb1cafe4491766f056b3fe58e32ac706))
* adding isnan to inequality check ([d780300](https://www.github.com/Loop3D/LoopStructural/commit/d7803005971e8e72663ebe2538251407079cb30b))
* adding node inequality constraints ([de33914](https://www.github.com/Loop3D/LoopStructural/commit/de3391403a7012895f7bf1bf170862ca2cfa6ed7))
* adding option to use mkl ([22883ba](https://www.github.com/Loop3D/LoopStructural/commit/22883ba59c13476ff77f91265268d10e56cd8a84))
* adding p2 interpolator for 3d :) ([5d0acfb](https://www.github.com/Loop3D/LoopStructural/commit/5d0acfb99458f2cacfbbfe099164e32ebe2a6d3b))
* bug with equality constraints ([b192cff](https://www.github.com/Loop3D/LoopStructural/commit/b192cff05987b6ad9350ade53208943afca18823))
* changing base unstructred support to work ([50ca16b](https://www.github.com/Loop3D/LoopStructural/commit/50ca16b70ddafd29ead061e153fe2b0652ca0e2b))
* correct ordering of shape functions ([6f2f113](https://www.github.com/Loop3D/LoopStructural/commit/6f2f113abb42a14c4b5278a635066f1696bf062d))
* minor fix to variogram saving wavelength guess as attribute ([2dbc29e](https://www.github.com/Loop3D/LoopStructural/commit/2dbc29e0696fe25bc3f55afb6cee1bea7bb9d9f1))
* modifying for 3d/generic case. ([c6a8ea3](https://www.github.com/Loop3D/LoopStructural/commit/c6a8ea3dcfdb7098d9e77df83f87e0f50036b670))
* norm constraints p2 ([93d3165](https://www.github.com/Loop3D/LoopStructural/commit/93d31654a9403f285a5ce3e4cfc3c6b4f91becf5))
* region funct type incorrect causing ([a5d27f1](https://www.github.com/Loop3D/LoopStructural/commit/a5d27f1569e0337971cc4997c673333ff1842d83))
* storing constraints as numpy array ([6c75379](https://www.github.com/Loop3D/LoopStructural/commit/6c753799019e19bb69d8621d5020787d5c42339f))


### Documentation

* :memo: adding warning that analysis module is experimental ([0574e57](https://www.github.com/Loop3D/LoopStructural/commit/0574e57bc639fe809b3e9b02531ce0e10742ea36))
* :memo: updating documentation for evaluate model ([1dcc716](https://www.github.com/Loop3D/LoopStructural/commit/1dcc716c6ae5cb0c77a54bf794ee14e077a874ce))

### [1.4.2](https://www.github.com/Loop3D/LoopStructural/compare/v1.4.1...v1.4.2) (2022-02-07)


### Bug Fixes

* actually fixing divide by zero ([ef4d7d0](https://www.github.com/Loop3D/LoopStructural/commit/ef4d7d036a24730a5c296f8e7779d76e405b98ce))
* added function for plotting structural frames ([69dbc6d](https://www.github.com/Loop3D/LoopStructural/commit/69dbc6d12999952d6266d003e8e7bda37551670c))
* adding mpl syntax for 3d view points ([cc7788f](https://www.github.com/Loop3D/LoopStructural/commit/cc7788f62f82b4bdc4bba242c394212ed94e6d5d))
* bug where data is not float ([e6f7e1f](https://www.github.com/Loop3D/LoopStructural/commit/e6f7e1fb6c09d4702bc49218d5b7b4e73385de7a))
* bug with calculating average fault ori ([a2099cd](https://www.github.com/Loop3D/LoopStructural/commit/a2099cd0477a18a6448d5e556795b0e40a1f05a6))
* commenting out osqp ([c2fd2da](https://www.github.com/Loop3D/LoopStructural/commit/c2fd2da4c5f6b84424cbd0fa0f3e940b7cadeb4f))
* crashing when all constraints of a type outside of box ([7907f58](https://www.github.com/Loop3D/LoopStructural/commit/7907f5814a4e2a0672c3fb73616325e533549809))
* extra checks ([ba6f391](https://www.github.com/Loop3D/LoopStructural/commit/ba6f39134e02daf996e59fd321d2698ce1dee040))
* fault interpolation taking a long time ([78fa916](https://www.github.com/Loop3D/LoopStructural/commit/78fa9163bc0b5fd84ef93302541dbf2da474749b))
* faults interpolation mesh was not being adjusted ([33722f6](https://www.github.com/Loop3D/LoopStructural/commit/33722f6c103701a7c7bf5dadb88c3bd422752622))
* faults not rescaling when incorrect norm vectors used ([9706f25](https://www.github.com/Loop3D/LoopStructural/commit/9706f255c612af4ea702c89f870826eed380167b))
* faults not visible in model ([32cbb2a](https://www.github.com/Loop3D/LoopStructural/commit/32cbb2a849c7f0c87da4adc2485013973eb2380b))
* faults with no value data not interpolating ([bf55ec6](https://www.github.com/Loop3D/LoopStructural/commit/bf55ec684cd4635f7396a4245a79b2e9adf81697))
* interpolator type is now an enum ([8197552](https://www.github.com/Loop3D/LoopStructural/commit/81975524f9ee487d8de489501b5fd2455150d20f))
* norm constraint out of model cause crash ([20692e2](https://www.github.com/Loop3D/LoopStructural/commit/20692e21dd7b97911e1b38f66ba405dfc152cb5f))
* number of constraints warning message ([9abcc04](https://www.github.com/Loop3D/LoopStructural/commit/9abcc042947b06f80db121aa48f9a9c1f0648b81))
* removed divide by 0 error ([9f17489](https://www.github.com/Loop3D/LoopStructural/commit/9f17489cf8f233651778e817f7fd9c3258bc7d9f))
* updating logging ([f77f915](https://www.github.com/Loop3D/LoopStructural/commit/f77f9156dfae6a58ea4a31b4b577dd2f207d2b8a))

### [1.4.1](https://www.github.com/Loop3D/LoopStructural/compare/v1.4.0...v1.4.1) (2021-12-06)


### Bug Fixes

* updating github workflow  ([84b4017](https://www.github.com/Loop3D/LoopStructural/commit/84b4017030089b51e81b34e2faa52f9526ad6c74))
* updating version number manually  ([1f11571](https://www.github.com/Loop3D/LoopStructural/commit/1f11571cc743f6b366898471ce4dbd4dbffd337b))

## [1.4.0](https://www.github.com/Loop3D/LoopStructural/compare/v1.3.12...v1.4.0) (2021-12-06)


### Features

* bumping version ([e50a35e](https://www.github.com/Loop3D/LoopStructural/commit/e50a35eaef77d873794a1d0efa362601882152ff))

### [1.3.12](https://www.github.com/Loop3D/LoopStructural/compare/v1.3.11...v1.3.12) (2021-12-05)


### Bug Fixes

* adding polarity to LS datastructure ([61842f0](https://www.github.com/Loop3D/LoopStructural/commit/61842f0b2f4bf72fd53e0d5834a89e2a996e4f1a))
* error creating mesh for faults ([7efe910](https://www.github.com/Loop3D/LoopStructural/commit/7efe910e3afc716d6cd7de69eb2c534cf9d8a1b0))
* gradient norm used for folds when fold ([36f9cff](https://www.github.com/Loop3D/LoopStructural/commit/36f9cff04ad6dde5fdefdf3feead3984a21e3fc4))
* splot not working for overturned folds ([2e185b5](https://www.github.com/Loop3D/LoopStructural/commit/2e185b52a5544723aa8d34f5e87d6449acf39668))

### [1.3.11](https://www.github.com/Loop3D/LoopStructural/compare/v1.3.10...v1.3.11) (2021-11-29)


### Bug Fixes

* added fault stratigraphy relation to processor ([ac864d8](https://www.github.com/Loop3D/LoopStructural/commit/ac864d86b4cd670f00dc60a59782d144536268d4))
* added project file processor ([855a690](https://www.github.com/Loop3D/LoopStructural/commit/855a69088fb49d18407a2a4404e775b71c2dc0bb))
* bugfix crash when evaluating feature ([d24f656](https://www.github.com/Loop3D/LoopStructural/commit/d24f656f75f955c46c706619daad5212a3a48160))
* fault support resized for all faults ([63ea5ec](https://www.github.com/Loop3D/LoopStructural/commit/63ea5ec89e210cc986c02d526cedd4590d374679))
* formatting ([a2f4fee](https://www.github.com/Loop3D/LoopStructural/commit/a2f4fee13cfd825e9ffd0800f61f6ee452e5bb61))
* formatting ([06d7e32](https://www.github.com/Loop3D/LoopStructural/commit/06d7e327bf59de23340018d1a842ba1b21deb0e4))
* missing argument from test + typo ([dee158a](https://www.github.com/Loop3D/LoopStructural/commit/dee158aed45624d8a711a5efb665c9555f1b4712))
* multiple faults not applied to stratigraphy data ([d7f4689](https://www.github.com/Loop3D/LoopStructural/commit/d7f468944c6300424caefbaf3e59911dde71cae5))
* only host version built on action ([8c7d058](https://www.github.com/Loop3D/LoopStructural/commit/8c7d058f9526a648e06b2f25c32ecddcc0233bb8))
* typo in variable names, AAT instead of ATA ([6a1a98f](https://www.github.com/Loop3D/LoopStructural/commit/6a1a98ff7fbf1e16c28cc71f9e167633af429937))
* updating github ci for conda build ([3525ec1](https://www.github.com/Loop3D/LoopStructural/commit/3525ec1a602d2b3a8b735f2b8d9c8b380f87331f))

### [1.3.10](https://www.github.com/Loop3D/LoopStructural/compare/v1.3.9...v1.3.10) (2021-11-10)


### Bug Fixes

* conda build for all python ([04757b5](https://www.github.com/Loop3D/LoopStructural/commit/04757b506a095a3d53b914f951305c9f34da1a0d))


### Documentation

* fixing documentation autobuild ([54e0bc9](https://www.github.com/Loop3D/LoopStructural/commit/54e0bc934a12cdda56c8d5a2112f7f61b2c68412))

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
