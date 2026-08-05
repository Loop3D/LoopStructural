# Changelog

## [3.3.1](https://github.com/Loop3D/map2loop/compare/v3.3.0...v3.3.1) (2026-01-31)


### Bug Fixes

* use first line from shapely and add all dip/thickness estimates ([#246](https://github.com/Loop3D/map2loop/issues/246)) ([381d13f](https://github.com/Loop3D/map2loop/commit/381d13fb2b281568f81fd959139f6b3ac292df18))

## [3.3.0](https://github.com/Loop3D/map2loop/compare/v3.2.10...v3.3.0) (2026-01-22)


### Features

* AlongSection Thickness Calculator ([#240](https://github.com/Loop3D/map2loop/issues/240)) ([ba8f82a](https://github.com/Loop3D/map2loop/commit/ba8f82af0b1663c032bfb80babf8f7785a1bf794))

## [3.2.10](https://github.com/Loop3D/map2loop/compare/v3.2.9...v3.2.10) (2026-01-20)


### Bug Fixes

* add ability to specify fault field to topoloy ([92f754b](https://github.com/Loop3D/map2loop/commit/92f754b2511a18f41c0f990ca6f6391774858ad1))
* add call methods to allow for generic debug exporting ([3c0e0e0](https://github.com/Loop3D/map2loop/commit/3c0e0e075b9d91bc32341cc8e2f249a828b00ee4))
* add option to set logger level and handler. ([5fb0720](https://github.com/Loop3D/map2loop/commit/5fb0720d795a7123094c39f171c3f497f53eda07))
* adding safeguarding and fixing merge ([a675be2](https://github.com/Loop3D/map2loop/commit/a675be2d6cef5646a682dbcd5c706681830d1ad9))
* applying copilot suggestions ([4021640](https://github.com/Loop3D/map2loop/commit/40216402634b5b9349719b681fb48525080b48a8))
* beartype issues ([6e2ee69](https://github.com/Loop3D/map2loop/commit/6e2ee69f08d913ce3a55df3f08144ded01922829))
* change to no default arguments and add unit_name_column as a required ([65a7b86](https://github.com/Loop3D/map2loop/commit/65a7b86ae94c44c098aceb17b13bf3fa624716b9))
* only store geodataframe if layer is not empty ([46fade8](https://github.com/Loop3D/map2loop/commit/46fade8bc8a673c0b8c067115970da70fed209bc))
* project import for test ([6eb1518](https://github.com/Loop3D/map2loop/commit/6eb151882369e3e7a062f60722aac9b7e3436736))
* replace args for kwargs ([e1cb1e6](https://github.com/Loop3D/map2loop/commit/e1cb1e642ffc9a3c206d51a8ccf78a064d645b5f))
* save location tracking as a line geodataframe for qgis compatibility ([17e14ca](https://github.com/Loop3D/map2loop/commit/17e14cae5cd466e8e3185e2569c311a611e0bed4))
* set age based sorter args ([5226b49](https://github.com/Loop3D/map2loop/commit/5226b49cfeda5592d878fc07e5ff960ba125f51e))
* update maximise contacts to compute length and check for empty contacts ([f2cf211](https://github.com/Loop3D/map2loop/commit/f2cf211ce61e9547ac9a00e4d70f00cba24f3d4d))
* updating incorrect type hints ([71f473c](https://github.com/Loop3D/map2loop/commit/71f473cd7e18d984b474de3a1b4bfd91bac98d81))
* updating logic ([86b658c](https://github.com/Loop3D/map2loop/commit/86b658c582c8ca6102729e8f8fd7370f28a113f3))
* updating requirements for alpha sorter plus add check for empty contacts ([b0f12b6](https://github.com/Loop3D/map2loop/commit/b0f12b68cfea2e0d16935720f01240230d912595))

## [3.2.9](https://github.com/Loop3D/map2loop/compare/v3.2.8...v3.2.9) (2025-12-01)


### Bug Fixes

* move extra sort arguments to sorter init ([#241](https://github.com/Loop3D/map2loop/issues/241)) ([3414fec](https://github.com/Loop3D/map2loop/commit/3414fec29969407f7b1b07ee59bfbedfe90a63b1))

## [3.2.8](https://github.com/Loop3D/map2loop/compare/v3.2.7...v3.2.8) (2025-10-29)


### Bug Fixes

* add is_strike parameter in thickness calculator ([#237](https://github.com/Loop3D/map2loop/issues/237)) ([1423a62](https://github.com/Loop3D/map2loop/commit/1423a62d4ccf54f7c860422b47223247af9b418c))

## [3.2.7](https://github.com/Loop3D/map2loop/compare/v3.2.6...v3.2.7) (2025-09-19)


### Bug Fixes

* DTM is optional ([#234](https://github.com/Loop3D/map2loop/issues/234)) ([92e0d24](https://github.com/Loop3D/map2loop/commit/92e0d24e7c1fa39c75a5fa6baf31704f91096bf0))

## [3.2.6](https://github.com/Loop3D/map2loop/compare/v3.2.5...v3.2.6) (2025-09-17)


### Bug Fixes

* refactor sorter ([#233](https://github.com/Loop3D/map2loop/issues/233)) ([3666bcb](https://github.com/Loop3D/map2loop/commit/3666bcb99e9a38c30ba87e94b3753f51a85b5c7a))

## [3.2.5](https://github.com/Loop3D/map2loop/compare/v3.2.4...v3.2.5) (2025-08-15)


### Bug Fixes

* handle optional columns 'INTRUSIVE' and 'SILL' in geology data extraction ([#230](https://github.com/Loop3D/map2loop/issues/230)) ([3cbc9e3](https://github.com/Loop3D/map2loop/commit/3cbc9e384068b9afb1530dae686398d4575dd064))

## [3.2.4](https://github.com/Loop3D/map2loop/compare/v3.2.3...v3.2.4) (2025-08-13)


### Bug Fixes

* trigger release for v3.3.0 ([54eab7c](https://github.com/Loop3D/map2loop/commit/54eab7ca462eb27cc1a6320ece8d1ecf2ebba2f5))
* trigger release for v3.3.0 ([54eab7c](https://github.com/Loop3D/map2loop/commit/54eab7ca462eb27cc1a6320ece8d1ecf2ebba2f5))
* trigger release for v3.3.0 ([8a4e093](https://github.com/Loop3D/map2loop/commit/8a4e093fedaefa4755db393e008bfc6484be0ba2))

## [3.2.3](https://github.com/Loop3D/map2loop/compare/v3.2.2...v3.2.3) (2025-06-20)


### Bug Fixes

* accept pathlib.Path for save_path ([4e2e799](https://github.com/Loop3D/map2loop/commit/4e2e7990eebaf0b47146632ac1311d1f96d9e9c3))
* allow str and pathlib.Path ([ff23f45](https://github.com/Loop3D/map2loop/commit/ff23f451b244d8bce05dafbd27a91c64fd9dcbe0))
* fold filtering to support multiple search terms ([c6d0b7a](https://github.com/Loop3D/map2loop/commit/c6d0b7aa2773bbe551fd71ec02d47f31417696f4))
* remove duplicate functions ([f629b2d](https://github.com/Loop3D/map2loop/commit/f629b2df62bc9da7ea1409e619107e580a8d052d))
* support multiple search terms for faults ([ad0289f](https://github.com/Loop3D/map2loop/commit/ad0289fa6de5af9e311a9e21efa13e949ff8867a))
* use pathlib.Path for save_path ([628040d](https://github.com/Loop3D/map2loop/commit/628040d42780086f4ee1ff39e258545ef53171b9))


### Documentation

* comment typo ([1c42d2e](https://github.com/Loop3D/map2loop/commit/1c42d2e450793c40abd7ad3b0c60843089f9a7e8))

## [3.2.2](https://github.com/Loop3D/map2loop/compare/v3.2.1...v3.2.2) (2025-01-13)


### Bug Fixes

* add featureId when parsing fault_orientations ([#177](https://github.com/Loop3D/map2loop/issues/177)) ([924c2cf](https://github.com/Loop3D/map2loop/commit/924c2cf696688abe3dc1c8579753daeb2f1b45e4))

## [3.2.1](https://github.com/Loop3D/map2loop/compare/v3.2.0...v3.2.1) (2025-01-12)


### Bug Fixes

* include dependencies in site-packages  - issue [#169](https://github.com/Loop3D/map2loop/issues/169) ([#170](https://github.com/Loop3D/map2loop/issues/170)) ([b33532b](https://github.com/Loop3D/map2loop/commit/b33532b56473148433fd192e182aadee028dc875))

## [3.2.0](https://github.com/Loop3D/map2loop/compare/v3.1.13...v3.2.0) (2024-12-16)


### Features

* v3.2 ([#153](https://github.com/Loop3D/map2loop/issues/153)) ([7978841](https://github.com/Loop3D/map2loop/commit/7978841b7106faf478492fe20770f17d9e244fbb))
* Thickness calculators can now be used simultaneously with `project.set_thickness_calculator([StructuralPoint(), InterpolatedStructure()])`
* Minimum fault length is a key in the config dictionary, and if not provided, it is calculated from the bbox area
* `ignore_lithology_codes()` and `ignore_fault_codes()` are new parameters in the config_dictionary
*  legacy config files with legacy keys are not allowed
*  `map2loop.utils.update_from_legacy_file` function may convert the hjson files into an updated config dictionary

### Bug Fixes
* all bug fixes can be checked in ([PR#153](https://github.com/Loop3D/map2loop/pull/153))

## [3.1.13](https://github.com/Loop3D/map2loop/compare/v3.1.12...v3.1.13) (2024-11-19)


### Bug Fixes

* bump with vtag ([5ebc07b](https://github.com/Loop3D/map2loop/commit/5ebc07b0014dcbd6c2e053113e0f39c99d2ee0a6))
* change docs build to combine all dependencies in one file ([fb51025](https://github.com/Loop3D/map2loop/commit/fb510250526fc35ab265e8273949a2c51ba52359))
* lowercase m2l and add vtag to rp manifest ([7311e40](https://github.com/Loop3D/map2loop/commit/7311e402b0cb7d4c5efcb3272b8caebca00a901b))
* remove vtag setting temporarily to see if this is what's stopping version bump PR ([c73e6fd](https://github.com/Loop3D/map2loop/commit/c73e6fda01f42d76b21267f2feae76c1e5e799ce))
* revert and double check typo on release-please ([d455350](https://github.com/Loop3D/map2loop/commit/d4553504909a25935d2c7115b1830ab56a93a071))
* rp manifest update ([e0bb75b](https://github.com/Loop3D/map2loop/commit/e0bb75b9e70acb98d0706b7320ea87f6cfd72c80))
* test trigger release ([2d87633](https://github.com/Loop3D/map2loop/commit/2d876337b826730cd7246f86b79d85b73acfe33c))
* typo in rp config ([9b1ba12](https://github.com/Loop3D/map2loop/commit/9b1ba122da3ff720e8ecc2354d6031d2df8cb7c0))

## [3.1.12](https://github.com/Loop3D/map2loop/compare/v3.1.11...v3.1.12) (2024-09-11)


### Bug Fixes

* removing additional indents ([cdd3eab](https://github.com/Loop3D/map2loop/commit/cdd3eabecf6c22fb88cdf64a9f2bad72b324eba1))

## [3.1.11](https://github.com/Loop3D/map2loop/compare/v3.1.10...v3.1.11) (2024-08-09)


### Bug Fixes

* added check for collocated points ([#125](https://github.com/Loop3D/map2loop/issues/125)) ([2bef09f](https://github.com/Loop3D/map2loop/commit/2bef09fa3dbfecbbb012a1e393c6851e6c3bbced))
* issue 122 ([c643721](https://github.com/Loop3D/map2loop/commit/c6437215125aea4d3a33b4150611ecc55eeb97c8))
* linting issues ([94af4d9](https://github.com/Loop3D/map2loop/commit/94af4d9fe1bf689ddf269dfafa07a40bc5df545f))
* remove pip pin for LPF because LPF has incorrect version tags ([dd99e15](https://github.com/Loop3D/map2loop/commit/dd99e15d167154bd902c005e012360b5005f20cf))
* update str format ([c3339a6](https://github.com/Loop3D/map2loop/commit/c3339a63161e7945140f55e5adf394e585fa590a))

## [3.1.10](https://github.com/Loop3D/map2loop/compare/v3.1.9...v3.1.10) (2024-08-02)


### Documentation

* corrected title ([86148a5](https://github.com/Loop3D/map2loop/commit/86148a59de1b96cddcae22c4d470c096c5fdbd5b))

## [3.1.9](https://github.com/Loop3D/map2loop/compare/v3.1.8...v3.1.9) (2024-08-02)


### Bug Fixes

* adapt floats to new np & unpin ([7bdd346](https://github.com/Loop3D/map2loop/commit/7bdd346fc8a565f107ce6bcfe8ec687701baa028))
* adapt floats to new np & unpin ([eac8ab2](https://github.com/Loop3D/map2loop/commit/eac8ab2e49142dd85e3bb4e9a09fb98d03482bab))
* adding poppler dependency for gdal.... why? ([a716933](https://github.com/Loop3D/map2loop/commit/a71693380a72bc71445295ac392171a0727ebe08))
* linting issues ([b743de4](https://github.com/Loop3D/map2loop/commit/b743de4ab1a8765cea703d041a65b1248abad1fc))
* linting issues ([1bf6600](https://github.com/Loop3D/map2loop/commit/1bf6600a5fe12a65eb8b866376f01b5f87c8d710))
* pin numpy version ([8a55469](https://github.com/Loop3D/map2loop/commit/8a554693b0a968e578560735f5bbb02db13d8b36))
* pin numpy version ([7dfd91c](https://github.com/Loop3D/map2loop/commit/7dfd91cbf20d1022b00bf007e87a8b99cb3836b3))
* update np version ([4acb184](https://github.com/Loop3D/map2loop/commit/4acb1849e7d51b3e878e59661e98d14e58187291))
* update np version ([5614446](https://github.com/Loop3D/map2loop/commit/5614446c616479ab753b72ca1132fcd99344c94f))
* update numpy version ([0fa6537](https://github.com/Loop3D/map2loop/commit/0fa6537841f1fa464adddae663351795e3158db1))
* update numpy version ([b312176](https://github.com/Loop3D/map2loop/commit/b31217681d4f21ca557c9aabd5f1756b9e8e57ed))
* use numpy type checker to catch new numpy typing ([48ac7cb](https://github.com/Loop3D/map2loop/commit/48ac7cb06d69f8549c484af6222ca6a96ab849d3))
* use numpy type checker to catch new numpy typing ([678f25e](https://github.com/Loop3D/map2loop/commit/678f25ec77e6d8ab207b98e7cd32be706a191c35))

## [3.1.8](https://github.com/Loop3D/map2loop/compare/v3.1.7...v3.1.8) (2024-06-20)


### Bug Fixes

* pin LPF to 0.1.3 ([f3766a8](https://github.com/Loop3D/map2loop/commit/f3766a83d0f09f293429cc975933cbf79b01c5a2))

## [3.1.7](https://github.com/Loop3D/map2loop/compare/v3.1.6...v3.1.7) (2024-06-18)


### Bug Fixes

* adapt floats to new np & unpin ([1a08401](https://github.com/Loop3D/map2loop/commit/1a08401ad65d3897452b223de8717db8b4e5b4ab))
* pin numpy version ([16130c6](https://github.com/Loop3D/map2loop/commit/16130c68b5c5a25c78f44ba8d11242c97b0cd8c4))
* update np version ([615985c](https://github.com/Loop3D/map2loop/commit/615985ce2fab530354585a1715e0edea55f2ba27))
* update numpy version ([0479ea1](https://github.com/Loop3D/map2loop/commit/0479ea1789380895789fca6d25334fe6a1638f12))

## [3.1.6](https://github.com/Loop3D/map2loop/compare/v3.1.5...v3.1.6) (2024-06-13)


### Bug Fixes

* add invalid hex input handling ([40028a8](https://github.com/Loop3D/map2loop/commit/40028a87b06c7610095edb0b78ec451838ff85a6))
* add the suggestions ([dc85b2e](https://github.com/Loop3D/map2loop/commit/dc85b2e5d89475511904275fc482d3bfd3f33fd8))
* comment the code ([6e064b1](https://github.com/Loop3D/map2loop/commit/6e064b1a3cb53af1fa92ea51faa4917d2c50663c))
* fix for duplicate units ([4052eca](https://github.com/Loop3D/map2loop/commit/4052eca85825a8ddd9d3c2f6c57c3f838ae2dd3e))
* make sure all colours are unique ([f90159a](https://github.com/Loop3D/map2loop/commit/f90159a83691278ee7cfca07c3aba64a144816c8))
* small fixes - double check for integer & add test information at the beginning of the test scripts ([a916343](https://github.com/Loop3D/map2loop/commit/a916343fba5ec59602fab9892b7d67f86702723b))
* update the tests for new function names (hex_to_rgb) ([e1778b7](https://github.com/Loop3D/map2loop/commit/e1778b734ca13d4c422daacd29de41c505f92fea))

## [3.1.5](https://github.com/Loop3D/map2loop/compare/v3.1.4...v3.1.5) (2024-06-06)


### Bug Fixes

* add 2 more links to try from GA WCS if timing out & prints out where the DEM was downloaded from ([92f73a5](https://github.com/Loop3D/map2loop/commit/92f73a55ad998af76fe02fd09702d98492c6e431))
* add catchall exception to mapdata_parse_structure & comment code ([59c677c](https://github.com/Loop3D/map2loop/commit/59c677c5988c1e411723c61f3482dbd792ed19a7))
* add test for all structures less than 360 ([f08c42a](https://github.com/Loop3D/map2loop/commit/f08c42a7e455657c0d042e1346fa5ce5fdb98775))
* allow 2 minutes for server connection & add add available ga server link ([900a50d](https://github.com/Loop3D/map2loop/commit/900a50d34f44e50c508466d617e4ce8d85c0c316))
* avoid double imports ([3347751](https://github.com/Loop3D/map2loop/commit/334775120d8db88a82e70980f68940d024524687))
* create the tmp file function missing ([55688fe](https://github.com/Loop3D/map2loop/commit/55688fe85f7aca2678a15160a73a972427a3bf34))
* ensure all dipdir vals &lt; 360 ([cf21a6b](https://github.com/Loop3D/map2loop/commit/cf21a6ba8cc0c48dfb7cc442d59a6cd63678ab83))
* fix for altitude not being saved properly in LPF ([b2c6638](https://github.com/Loop3D/map2loop/commit/b2c663866f478752670aca41a188aa195df9c90f))
* make the templates a bit easier to fill out ([81f54fe](https://github.com/Loop3D/map2loop/commit/81f54fe9f7163dc997219b4fe092c81b6c7b3ca3))
* move location of the PR template ([9209c25](https://github.com/Loop3D/map2loop/commit/9209c251520bfd993a51b04124ad62a831a58922))
* return just column with new vals instead of inplace modification ([722b98c](https://github.com/Loop3D/map2loop/commit/722b98cdbcf1f820b9d1a480f09ac59d287f1a04))
* revert back to original ([6ed00c9](https://github.com/Loop3D/map2loop/commit/6ed00c999901c3236e53fc039930510be9d38570))
* revert back to original timeout ([301242f](https://github.com/Loop3D/map2loop/commit/301242f931b92ac01c4251ecf9039119d9da49d5))
* revert back to original url ([26e4971](https://github.com/Loop3D/map2loop/commit/26e497198744c591879fea96c10935d6b79b7d9d))
* update question template ([8182e50](https://github.com/Loop3D/map2loop/commit/8182e50361852de817e1440adcb278b5f46f4796))
* update tests that depend on server to skip in case server is unavailable ([646be9e](https://github.com/Loop3D/map2loop/commit/646be9ea06dade92f16208da7e681d6d2b0e6c65))
* use gpd.points_from_xy ([2f931c5](https://github.com/Loop3D/map2loop/commit/2f931c59025997133aeea9c10f725858613b7f6b))
* verbose dipdir 360 test name ([6ffe6bf](https://github.com/Loop3D/map2loop/commit/6ffe6bf31ab860ee054bb175a63c00fee7b8a7ff))


### Documentation

* update libtiff version ([0a99ac8](https://github.com/Loop3D/map2loop/commit/0a99ac8aaf6980196216137f15060c44b2f77cd9))

## [3.1.4](https://github.com/Loop3D/map2loop/compare/v3.1.3...v3.1.4) (2024-05-29)


### Bug Fixes

* actions updated to prevent linting issues ([9a3e86a](https://github.com/Loop3D/map2loop/commit/9a3e86a4e5c7e89f93e627081cc51df21c75ba5b))
* add pull request template ([f8ccac7](https://github.com/Loop3D/map2loop/commit/f8ccac74d5d778be0f675c683bf4e29d59027001))
* added issue_templates ([06413ad](https://github.com/Loop3D/map2loop/commit/06413ade7660f33da5245b00493119d77649a681))
* few typos grammar ([97bfe16](https://github.com/Loop3D/map2loop/commit/97bfe16cdf2653e957946e8366a1dc801bbc4b2f))
* fix for issue [#73](https://github.com/Loop3D/map2loop/issues/73) ([39c65a7](https://github.com/Loop3D/map2loop/commit/39c65a721af57038e73a26c0522b303c5a40c16b))
* keep things light ([50751a9](https://github.com/Loop3D/map2loop/commit/50751a940c78a981ec4c423e3c914615cede962f))
* merge PR templates into one single PR template ([dcbad85](https://github.com/Loop3D/map2loop/commit/dcbad8503174c9700651bb5c491323336d10c8b8))
* remove tini from docker ([fc86b89](https://github.com/Loop3D/map2loop/commit/fc86b89994b578b6d603a514ec1c889623092df7))
* remove value from issue templates ([11a005d](https://github.com/Loop3D/map2loop/commit/11a005db964bcf4d82749b5c334dea9a8e394cad))
* revert back to having strike allowance as discussed ([48a19e0](https://github.com/Loop3D/map2loop/commit/48a19e091e22bb0772c22bfa4ebc2b9ff2f0bffe))
* StructuralPoint actually prints warnings ([ab6f276](https://github.com/Loop3D/map2loop/commit/ab6f276a991f538dd25bca467e3869b283a2df10))
* tests should be with pytest ([a04fb2a](https://github.com/Loop3D/map2loop/commit/a04fb2a569ba756e5ab2317782fd96bc8d15fd0c))
* update wording in templates for clarity ([ecb1092](https://github.com/Loop3D/map2loop/commit/ecb10926b25472672d7929d94aec0fec222ad09f))


### Documentation

* remove lavavu; no need for 3D vis for the examples ([6b42cdc](https://github.com/Loop3D/map2loop/commit/6b42cdcbc9ca7fd9e10ceee31444196f8088ee4f))

## [3.1.3](https://github.com/Loop3D/map2loop/compare/v3.1.2...v3.1.3) (2024-05-14)


### Bug Fixes

* remove tini from docker ([89fd60f](https://github.com/Loop3D/map2loop/commit/89fd60f35445f96894c153d3f8c5f6451bcbf0ed))


### Documentation

* remove lavavu; no need for 3D vis for the examples ([e8b1ada](https://github.com/Loop3D/map2loop/commit/e8b1ada0345cfebef114be555e770fd594155b16))

## [3.1.2](https://github.com/Loop3D/map2loop/compare/v3.1.1...v3.1.2) (2024-05-09)


### Bug Fixes

* update doc deploy ([cf0de5b](https://github.com/Loop3D/map2loop/commit/cf0de5b937c251ac71522488392f7ba7d65d78c6))

## [3.1.1](https://github.com/Loop3D/map2loop/compare/v3.1.0...v3.1.1) (2024-05-09)


### Bug Fixes

* trigger build ([f1879ce](https://github.com/Loop3D/map2loop/commit/f1879ced58154610ccadd46488cdf7f109ca1810))

## 3.1.0 (2024-05-09)


### Features

* Add new stratigraphic column sorters based on contact length and unit relationships fix: transpose dtm for height calculations, export geotiff using gdal ([36ddcca](https://github.com/Loop3D/map2loop/commit/36ddcca650d7f12c102b77e5df23d7a2b7710ee3))
* added utils function ([3c01767](https://github.com/Loop3D/map2loop/commit/3c01767c0f3adde4ccae59d84e01e65849b30993))
* Allow for empty dataframes if no filename specified. Mainly for fault and fold files ([f81d9b4](https://github.com/Loop3D/map2loop/commit/f81d9b4496a9d863864f72d2c1d2813871bccb24))
* new ThicknessCalculatorBeta ([f703832](https://github.com/Loop3D/map2loop/commit/f7038326fb2ff291ab23e3c3e52e7ad653b68e1a))
* thickness calculator gamma ([6926a45](https://github.com/Loop3D/map2loop/commit/6926a4540198f8fb25dd2c8905be5c2072ea511a))
* two interpolators for thickness calculations ([e118f6a](https://github.com/Loop3D/map2loop/commit/e118f6afb7202595e09c6ef3798c28c4ea7035f6))
* two interpolators for thickness calculators ([11f3a04](https://github.com/Loop3D/map2loop/commit/11f3a0455103b464285b47878dc226916b7bc233))


### Bug Fixes

* 3.12 ([264b4c0](https://github.com/Loop3D/map2loop/commit/264b4c069fc069cd84d30c43d4931f5815391fc4))
* actually remove the rank & type lines ([16801d2](https://github.com/Loop3D/map2loop/commit/16801d23c52d74c291161a52ac2fbb541a133529))
* add back the project import ([01dc0bb](https://github.com/Loop3D/map2loop/commit/01dc0bbacfa605aa5aac12dca94aed4d158420fa))
* add comments to the gamma calculator ([936b495](https://github.com/Loop3D/map2loop/commit/936b495d10ea28c52349f66cb35638d79ed329ff))
* Add data files as package_data so that it installs with .py files ([09b0cb8](https://github.com/Loop3D/map2loop/commit/09b0cb8a72d9eed06c6a3de7c62a0931ceb6b846))
* add descriptive variable names ([672bece](https://github.com/Loop3D/map2loop/commit/672bece19ddd3e0dac775a47a4577a7d5c163d45))
* Add exception handling around hjson access ([#24](https://github.com/Loop3D/map2loop/issues/24)) ([f627f83](https://github.com/Loop3D/map2loop/commit/f627f83079ba0c1315fc59d98fc2bb2fe27a0191))
* Add exception handling around hjson access, throw with better message ([d8804de](https://github.com/Loop3D/map2loop/commit/d8804de89a9ecc8705914f959b3c4a4858f968d9))
* add feature to count segnum of polygons with holes ([03bcbcc](https://github.com/Loop3D/map2loop/commit/03bcbcc00b6df28a0dc0879d6fa41ae0bc355cf8))
* Add layer id to structural data using a spatial join ([#38](https://github.com/Loop3D/map2loop/issues/38)) ([09c53d6](https://github.com/Loop3D/map2loop/commit/09c53d6fbb6b44c46b9e88cf9044f7ccd1ebd555))
* add lazy function to load files from state name ([0e0bfa7](https://github.com/Loop3D/map2loop/commit/0e0bfa7f72f80e89aba5a548a7282743dc0bd04a))
* add LPF 0.1.0 ([a4165b3](https://github.com/Loop3D/map2loop/commit/a4165b3f81c1f9d2625f1268561583e8b967bd7f))
* Add multiple tries to online config file access ([8a137ce](https://github.com/Loop3D/map2loop/commit/8a137cedc298596ea4ce858a5ff4e256eb169b7c))
* add option for user defined line length ([23aeb06](https://github.com/Loop3D/map2loop/commit/23aeb06c0f422598e054727ab0d8357da76a7705))
* add options/warnings to overwrite existing lpf files ([165bc17](https://github.com/Loop3D/map2loop/commit/165bc17f7a8181f9d1bab717120a4f3e86dab81a))
* Add owslib as dependency for testing ([0439da8](https://github.com/Loop3D/map2loop/commit/0439da869856b187defbed631c592a848b3684a1))
* add Pathlib to aus state ([e9d7add](https://github.com/Loop3D/map2loop/commit/e9d7addcbf545d67a0547c0cd0bf5dc73b9748d8))
* add remaining thickness fields ([130af98](https://github.com/Loop3D/map2loop/commit/130af988df47ca27ad5aa3ecba67655517697430))
* add some comments throughout the script and remove unnecessary comments ([390b864](https://github.com/Loop3D/map2loop/commit/390b864545223e01bb7a475429d4d6a7298732d2))
* Add special case for SA config file to lowercase the column names.  Also check for empty map2model relationship files and deal appropriately with them ([82054df](https://github.com/Loop3D/map2loop/commit/82054dfec9f82b353a63b42e5f3d815fdeff5fee))
* add str and Pathlib to beartypes ([16a0247](https://github.com/Loop3D/map2loop/commit/16a02470b520b3a118c69e29a63e3bb7eeda25a2))
* add strike attribute to the class ([80e4445](https://github.com/Loop3D/map2loop/commit/80e4445bcfbe73f3a430c77dce0a96de56c913c6))
* add test datafiles to gitignore ([951e7fe](https://github.com/Loop3D/map2loop/commit/951e7fea8ea358a00426fd7cdcedf552e865acc2))
* Add test for featureId ([37f5a96](https://github.com/Loop3D/map2loop/commit/37f5a96dce30ee3115105701f1c2c96f58d3d9fe))
* add thickness fields to proj and thickness_calculators ([f6ec60f](https://github.com/Loop3D/map2loop/commit/f6ec60faf0d075fe8fa8fedcc60884840b9dc0f9))
* added function to find the datasets folder after install. ([9623d48](https://github.com/Loop3D/map2loop/commit/9623d48643f664683e9403d14208a1e6a991f949))
* Added new stratigraphic sorters, a placeholder for throw calculators and updated thickness calculator ([e08e010](https://github.com/Loop3D/map2loop/commit/e08e010f38070c175412f0b92946b5ce705c5409))
* added option for grid resolution ([a0ef130](https://github.com/Loop3D/map2loop/commit/a0ef130bb12148bab8cb2bb7bddf7fb33db4f085))
* added segNum to geology samples ([0ad6983](https://github.com/Loop3D/map2loop/commit/0ad69839dc25d7400e80a4ef2deed5e2ace7dcf3))
* adding beartype to requreiments ([4d66b8f](https://github.com/Loop3D/map2loop/commit/4d66b8fd76d8ba5a64abdc0838b9797dfcc75dca))
* adding empty samples to project class ([bf100b2](https://github.com/Loop3D/map2loop/commit/bf100b2c23a284072d11d54e080d882da0460942))
* adding fault orientation calculation ([e35f0b6](https://github.com/Loop3D/map2loop/commit/e35f0b6a78ae85482a66f6636b54bc079432dd51))
* adding fault orientation calculator ([08c9f5e](https://github.com/Loop3D/map2loop/commit/08c9f5e943ee31b52ca613c4da44b3521018eb22))
* adding fault orientation data type ([f0e82e2](https://github.com/Loop3D/map2loop/commit/f0e82e20553f7f930ea4237c7516a18c61907090))
* adding orientation type to fault config. ([6359937](https://github.com/Loop3D/map2loop/commit/635993739bac4acf5c6e6ad83a884c391602285c))
* adding random colour generator ([e6315c4](https://github.com/Loop3D/map2loop/commit/e6315c47b9f3c18140da62f3d124021e79ae403b))
* adding to_dict for config ([968504e](https://github.com/Loop3D/map2loop/commit/968504e27e56b343feee58a8a1ed297f493c5d87))
* adding to_dict for config ([968504e](https://github.com/Loop3D/map2loop/commit/968504e27e56b343feee58a8a1ed297f493c5d87))
* adding to_dict for config ([258254e](https://github.com/Loop3D/map2loop/commit/258254e7fc4e2628c035a400083c04dc4f81866b))
* adding version and ignoring init file for f401 ([db6ac5b](https://github.com/Loop3D/map2loop/commit/db6ac5b48f7ade6d57b2a6a094811c083c996bd0))
* adding version attribute ([1a2ed67](https://github.com/Loop3D/map2loop/commit/1a2ed677443ae2497aee9bc6b0a929e7566632b7))
* avoid future warnings from geopandas ([e877d93](https://github.com/Loop3D/map2loop/commit/e877d936426fca7b58f9825e8f82687660edf5c2))
* can access sampled data for thickness calculations ([1ae2d2d](https://github.com/Loop3D/map2loop/commit/1ae2d2d199806c2d15e828543f890a73b56909ee))
* change config file name convention for clarity ([15123ac](https://github.com/Loop3D/map2loop/commit/15123ac68a590b85010cd766f588a4fa859b4d05))
* change docker to see if documentation builds ([10f4fde](https://github.com/Loop3D/map2loop/commit/10f4fdea4a9a856ce30941cb5f1878e63d226870))
* change make html statement in docs build ([26aeba6](https://github.com/Loop3D/map2loop/commit/26aeba65d10bc4251f97d584a15ded9c04ad5b74))
* change segNum to featureId ([aab8c8c](https://github.com/Loop3D/map2loop/commit/aab8c8c90c42dbe3a85431ed669f386a187a4349))
* change segNum to featureId ([2689a3a](https://github.com/Loop3D/map2loop/commit/2689a3a7307ae49a7f9198be7365471c43ccdad2))
* changing boolean comparison ([d90df18](https://github.com/Loop3D/map2loop/commit/d90df183cb9042a994525cfd698297aa588ea04f))
* changing isinstance to issubclass ([0e9ecc8](https://github.com/Loop3D/map2loop/commit/0e9ecc82b37a0a979aefc8fb3ae83b406f30a37a))
* Check for empty dataframe before concat ([c475513](https://github.com/Loop3D/map2loop/commit/c475513c84c356d1a2b09bcd331c861f124b3317))
* correct typo on actions & reverse test path ([a5fbdbc](https://github.com/Loop3D/map2loop/commit/a5fbdbcc394629a4d77bc04b2e7f65c4e0aa05e8))
* data works without init files on subfolders ([cc833bb](https://github.com/Loop3D/map2loop/commit/cc833bbee451f1b540237cc9fa6abc3c3747979d))
* default fault dip/dipdir to nan ([917cee6](https://github.com/Loop3D/map2loop/commit/917cee6bc7c27ed18a067403e7780d54656f92f5))
* default fault values to nan not 0 to prevent unexpected behaviours ([87863c7](https://github.com/Loop3D/map2loop/commit/87863c76654da4fee87fed45ac5015e30c49defa))
* Default temp path changed to work in linux ([5bd773b](https://github.com/Loop3D/map2loop/commit/5bd773bdd27268aa8b516b799a16dca141056c0f))
* Default to not showing file access debug info ([fbcf10d](https://github.com/Loop3D/map2loop/commit/fbcf10df02eb98fdeb5c2dfafffe5624835eac60))
* delete runtime from repo ([36d59bc](https://github.com/Loop3D/map2loop/commit/36d59bce5bd7a8dfedd279c877a1d9b7ee6bf18b))
* Dependencies issue with new version of geopandas breaking read_file on remote zip files. Also, fixed errant line that installed the old version of map2model ([4eb720a](https://github.com/Loop3D/map2loop/commit/4eb720a2d766b821d0cf31e46bc8b3f571ea6fc7))
* edit path of config filename ([6864d7d](https://github.com/Loop3D/map2loop/commit/6864d7d1efa204ff286b9ff4533103d984d5ec9b))
* endpoints have correct length ([c36995d](https://github.com/Loop3D/map2loop/commit/c36995dc2a31ce16d64e3ef46a721c53ee8c602c))
* fix typo ([f69b42f](https://github.com/Loop3D/map2loop/commit/f69b42f8640827df2b7ca7e35629aac729d34684))
* fixed NaN values and NaN warnings ([f7cbd78](https://github.com/Loop3D/map2loop/commit/f7cbd78cae78043916b79ada02d2c6677408db1b))
* Forgot that env varibles can't have hyphons ([c31bb1c](https://github.com/Loop3D/map2loop/commit/c31bb1cfe33b992ad9f6bebb3d406bd1f5563d45))
* formatted using ruff ([06ce94d](https://github.com/Loop3D/map2loop/commit/06ce94d6d4926eea7743341bbb17e766e3156bfa))
* init file recovered ([f9f56cc](https://github.com/Loop3D/map2loop/commit/f9f56cccab3eb310d6a80058305b4ca19746264c))
* install pytest with gh actions ([7b25639](https://github.com/Loop3D/map2loop/commit/7b25639310e03186c4db3df7312070825b034a2b))
* Issue with conda build test env not working, skipping ([2fffb69](https://github.com/Loop3D/map2loop/commit/2fffb696556f48f0620dfbac356971301fb15e89))
* letting docs build run on docker ([8b4ade0](https://github.com/Loop3D/map2loop/commit/8b4ade0675aa5d497cbd469a85f59b1d19b09239))
* linting ([c1bf893](https://github.com/Loop3D/map2loop/commit/c1bf893ee6d2316d153f4e0fb66f57f90a7dc9af))
* linting ([4cf51b7](https://github.com/Loop3D/map2loop/commit/4cf51b73192972c826e74373181e11883455b2a0))
* make sure Pathlib builds the right path ([f10270d](https://github.com/Loop3D/map2loop/commit/f10270d2666c75b38f936a94800532edc630ff53))
* minor adds/delete ([e479875](https://github.com/Loop3D/map2loop/commit/e4798750c6b13ac74214d639add556dce6111e9f))
* minor changes to the test for gamma thickness calc ([9211826](https://github.com/Loop3D/map2loop/commit/921182667babe99bca6e1bed3df1ef7dac073504))
* more updates to agree with Pathlib ([7a19f34](https://github.com/Loop3D/map2loop/commit/7a19f345db25a2079b9768722619b72a08b42f99))
* Move clut and config files to _datasets directory and use local files for doc-test case ([23e6b8b](https://github.com/Loop3D/map2loop/commit/23e6b8b2d85c6bc5e29d0d766d23d981279e8a0e))
* neat output details ([085d4eb](https://github.com/Loop3D/map2loop/commit/085d4ebbb213a800bb9850e1f496aa4c3993c7e8))
* new flake8 compliance ([6776e25](https://github.com/Loop3D/map2loop/commit/6776e25a4c22faf483ffd70010202211b0b874b2))
* Pass featureId of multi lines to project file, also unpin geopandas ([1a4299e](https://github.com/Loop3D/map2loop/commit/1a4299e16f616a8286f241262dbccfea0f069586))
* Path typo ([883b497](https://github.com/Loop3D/map2loop/commit/883b4978be6f0479c17ee50038fa3f2577628ca9))
* Remove deprecation ignore warnings and then deprecations ([79cd333](https://github.com/Loop3D/map2loop/commit/79cd33317625c96091690c9b97309be74cba2fdd))
* Remove j in loop as it was not referenced ([3314828](https://github.com/Loop3D/map2loop/commit/33148287a8c978c775d016eeb293923907c9f3c6))
* remove print statement ([ecd1e44](https://github.com/Loop3D/map2loop/commit/ecd1e4449fe3736f555e9cffe6dbae64e61017ba))
* remove segNum_test ([cbd05f2](https://github.com/Loop3D/map2loop/commit/cbd05f22256ce1b2d115112b3921a8bbf25443a3))
* remove this test for now and see if everything else is ok. ([16a3b88](https://github.com/Loop3D/map2loop/commit/16a3b882d82ebef47fc814dd5e9bfaf5729d0082))
* remove unnecessary else stat ([bea7e7a](https://github.com/Loop3D/map2loop/commit/bea7e7aa26f28691b325617b8439704dae093d30))
* remove unnecessary print statements ([58535bf](https://github.com/Loop3D/map2loop/commit/58535bf07ec3d172dbf0881a9c0cc3be95b042f8))
* remove unnecessary variable declarations ([fb44fef](https://github.com/Loop3D/map2loop/commit/fb44fef387ef521d0b5e0f48f205144a35190ba4))
* remove unsused variable ([fa3e4bb](https://github.com/Loop3D/map2loop/commit/fa3e4bbc71e939ff0836a55166268d32ebc5b866))
* removed float type hint ([87713d6](https://github.com/Loop3D/map2loop/commit/87713d62dbcd086b7662cb1d6145b048729bae3d))
* removed rank & type from LPF obs tables ([750acbd](https://github.com/Loop3D/map2loop/commit/750acbd505f901c74b892b1df261893fb414d681))
* removing commented code ([22db166](https://github.com/Loop3D/map2loop/commit/22db1668448e088fa02caa4f5bd2e17980774dc2))
* rename in thicknessCalc names in test ([fda68e6](https://github.com/Loop3D/map2loop/commit/fda68e686bc36195e95e24cb733d3a2475b90b54))
* Rename segNum to featureId ([274beca](https://github.com/Loop3D/map2loop/commit/274becaad09bb8500f29edc097c1d877bd7487e4))
* renaming compute to calculate ([f5a0ffb](https://github.com/Loop3D/map2loop/commit/f5a0ffb3beeb6c59dc2c3758bcc7251a74cea9a1))
* replace beta and gamma nomenclature ([a9040b0](https://github.com/Loop3D/map2loop/commit/a9040b0f365b7035d67af7846ad5e3ca4a54574a))
* revert hjson back to loop3d version as hjson-py does not work ([84b846a](https://github.com/Loop3D/map2loop/commit/84b846ae8efd293428acb2bb1c4515f0f48e6a31))
* revert naming conversion to original ([5e7e5ec](https://github.com/Loop3D/map2loop/commit/5e7e5ec43d5cb557590205893db1953be8b74545))
* Reworked dataframe to avoid chained assignment ([77b322f](https://github.com/Loop3D/map2loop/commit/77b322f5daccbf74d316bd6adfd3e1dd70e47d5e))
* ruff linter and docs fix ([c7e42a6](https://github.com/Loop3D/map2loop/commit/c7e42a68e4824ee1eff7a0c8c8797e1c2b7365db))
* run ci on all branches but only release from master ([9e0f977](https://github.com/Loop3D/map2loop/commit/9e0f977a4d4dbee51ef357db638b339810863b25))
* sampled_contacts are in map_data now ([9d57d95](https://github.com/Loop3D/map2loop/commit/9d57d95f6bca4a2f5e6ce4f0ab9fbcd5fed647a7))
* setup_grid docstring updated ([67b453b](https://github.com/Loop3D/map2loop/commit/67b453b5a32239004936ff1d31dc6a094f068975))
* simplify output generation; assignment of -1 to uncalculated thickness ([6acd2d5](https://github.com/Loop3D/map2loop/commit/6acd2d5eea382fa9fcf72aaefb0c628014681bce))
* simplify the procedure to get the relative path of the library ([416d101](https://github.com/Loop3D/map2loop/commit/416d101d73cf850381dc2169d5f6569f0f175b69))
* Still forcing versioning for release-please ([08e647b](https://github.com/Loop3D/map2loop/commit/08e647bac1dded0c807bbb9a3571a741505bd488))
* store xyz end points correctly ([af29ab6](https://github.com/Loop3D/map2loop/commit/af29ab6b6bb6937baa836cfe73e7003ebecd11a4))
* test changing dirs to see if docs build ([98098fc](https://github.com/Loop3D/map2loop/commit/98098fcd3c829a075de884bb613e331c6e9b0e1e))
* testing change tests directory to map2loop folder, because of json file location ([371a2a2](https://github.com/Loop3D/map2loop/commit/371a2a2e54aeb4b5619c7daf16b85123e71ee74e))
* testing test location again ([c775ed9](https://github.com/Loop3D/map2loop/commit/c775ed9ce45e6653e9b752b59a8592246df848fb))
* This environment variable has to go somewhere ([2495a89](https://github.com/Loop3D/map2loop/commit/2495a89093f5808cd75c184b442576bfdc1083ee))
* Try again setting environment variable for documentation test case ([4e6348b](https://github.com/Loop3D/map2loop/commit/4e6348bfebd6131db471cd9273c34330165f0e99))
* try getting colour data, if file doesn't load print warning ([dd4f365](https://github.com/Loop3D/map2loop/commit/dd4f36517d6fdeb5f04d7aee571260129e69d019))
* typo ([025de10](https://github.com/Loop3D/map2loop/commit/025de10e07fac9e045e516a2f6311a51b2fed1b1))
* typo ([ea6d0cd](https://github.com/Loop3D/map2loop/commit/ea6d0cd692fb47adec6403d4e16edb9ce6a1eb00))
* Typo fix for config dictionary converter ([15b543b](https://github.com/Loop3D/map2loop/commit/15b543bc077550c683c94d7bb97a74389b1de29e))
* typo in calculate ([c0ec36b](https://github.com/Loop3D/map2loop/commit/c0ec36b2444fb136b80e0356debf5f13ad454258))
* Typo in intrusion flag, ignoring intrusions for contact calculations, remove errant sorter from take_best option ([527de6c](https://github.com/Loop3D/map2loop/commit/527de6c339e56ead036ab5ed5d94d06d77673ae6))
* upd v3.1 ([182079b](https://github.com/Loop3D/map2loop/commit/182079ba13693d4bd3b6a8373c5954076f36ec8a))
* update (thanks RT & LG) ([f6f59cf](https://github.com/Loop3D/map2loop/commit/f6f59cf2dc2cbcf40351745cba1359db9c8db552))
* update all to overwrite LPF ([d0e38db](https://github.com/Loop3D/map2loop/commit/d0e38db735984103a4e95d43ffd5dc974ae6755f))
* update aus clut files with documentation links ([55bc845](https://github.com/Loop3D/map2loop/commit/55bc84516224001e29fcc79b7663cffdb90d5655))
* update clut file names ([c2c020f](https://github.com/Loop3D/map2loop/commit/c2c020fa0db8a625cc3b9095c231a9ebd22975d2))
* update code as per Roy's review ([450b372](https://github.com/Loop3D/map2loop/commit/450b372cfc43063a2de2eebf2b6e05458114c425))
* update docs dockerfile build ([deefe49](https://github.com/Loop3D/map2loop/commit/deefe4933db82b8a6633deb154eab6ed26582a5c))
* update for python 3.12 ([23c3d76](https://github.com/Loop3D/map2loop/commit/23c3d76f06a7d195b083f1eafb793fa147554ed9))
* update gh actions to revert to roiginal test location ([9b341ad](https://github.com/Loop3D/map2loop/commit/9b341ad0c4658cfe909a981e70c80b05358ce5c7))
* update gitignore ([d2df267](https://github.com/Loop3D/map2loop/commit/d2df2670cd4e3153def613c0c56e294326495a8e))
* update LPF dependencies ([cbb5489](https://github.com/Loop3D/map2loop/commit/cbb54899368ae5950b29a101bfca648b900d6527))
* update test to check segnum & geology ID and spatial match ([7cbf143](https://github.com/Loop3D/map2loop/commit/7cbf143b84a4c7cb3f34b563d2f9ef3202b0433e))
* update tests to run from upper folder ([1a8210d](https://github.com/Loop3D/map2loop/commit/1a8210dc41942164b202b0cfbea044091c67e5b0))
* update the rest of the variable names to thicknessMedian except in LPF ([2f8816b](https://github.com/Loop3D/map2loop/commit/2f8816becbb090219717cd8dde0ad6a390753c0f))
* update the StructuralPoint test ([dc7c6da](https://github.com/Loop3D/map2loop/commit/dc7c6dafab01cef2dd215a5df52f6e7af9391463))
* update to new thickness variable name ([5f79706](https://github.com/Loop3D/map2loop/commit/5f797064cf557f0cc47c32bdd6e73dd60b20853e))
* update to pathlib in documentation test environment ([fd263c2](https://github.com/Loop3D/map2loop/commit/fd263c20e86b9b58a0f44a2ed84cca912e8660e9))
* update to solve the Pathlib ([5b0ace7](https://github.com/Loop3D/map2loop/commit/5b0ace70e96f6ec3c5f70873579d474c1232761b))
* update variable names for thicknesscalculator StructPoint ([f6a2207](https://github.com/Loop3D/map2loop/commit/f6a22076c174e4eed0c17db0a38a62d386b91e23))
* updated docstring ([93d47bf](https://github.com/Loop3D/map2loop/commit/93d47bfac84397cf842a9398d1b0bd0dae873c5b))
* updated type hint for basal_constacts ([c125c63](https://github.com/Loop3D/map2loop/commit/c125c6340f822fbe13f5839f8fda3e8463da44b9))
* use correct variable rbf ([79791f1](https://github.com/Loop3D/map2loop/commit/79791f125dfe0748b80e243ef79f86ad10197752))
* Use geology sampler for contact sampling ([615d19a](https://github.com/Loop3D/map2loop/commit/615d19adefddeb736871ed4b4ad5329e8b1ac9e8))
* Use MANIFEST.in to include clut and config files in installation ([d49a53f](https://github.com/Loop3D/map2loop/commit/d49a53f526ba7b4d14b7f88f5091986054717841))
* use Nx3 array instead of vector components ([ae6739e](https://github.com/Loop3D/map2loop/commit/ae6739eaceb0d0b4aabb2650a593e799f11c4a24))
* use shapely.object to avoid clashes ([73dd77f](https://github.com/Loop3D/map2loop/commit/73dd77f4f911eb10237819c862382bb36481fc34))


### Documentation

* add 3.11 back on so docs build (LPF needs 3.12 as well) ([9f817be](https://github.com/Loop3D/map2loop/commit/9f817beba63e581a85145fea2521bb95c6942c29))
* add back the .sh to build docs locally ([a5ccb0f](https://github.com/Loop3D/map2loop/commit/a5ccb0fdf6557dba3b38245c49833128b3093f0c))
* add docstrings for ThicknessCalculatorGamma ([39187f7](https://github.com/Loop3D/map2loop/commit/39187f773f2e5d4647992630d4bf25c7eee2897d))
* Add import modules and example parameters to sorter and thickness calculator docs. ([3ed1ea7](https://github.com/Loop3D/map2loop/commit/3ed1ea73128bcce142fbd84a38947fd99b24332e))
* Add import modules and example parameters to sorter and thickness calculator docs. ([#33](https://github.com/Loop3D/map2loop/issues/33)) ([b4e6247](https://github.com/Loop3D/map2loop/commit/b4e6247117baf012009c7a81c4b923beb38bef5c))
* add m2l verson ([c5cf102](https://github.com/Loop3D/map2loop/commit/c5cf10243f0824731649623c423078cb159e61c7))
* add space to avoid warning ([d615498](https://github.com/Loop3D/map2loop/commit/d615498c6356fdeb239f0bd631ceb8767998c9c5))
* added thickness to the base class. ([e263389](https://github.com/Loop3D/map2loop/commit/e2633892dec2489609bc79ec7992ee99f597fd99))
* adding aylas code template ([48c1093](https://github.com/Loop3D/map2loop/commit/48c10935654b35aef191c160aa7c2cb69dc38332))
* adding clut presets to docs ([65b23fd](https://github.com/Loop3D/map2loop/commit/65b23fdaf935711c7ed6f775176840eaaa474a7c))
* adding doc specific requirments ([b26d4dd](https://github.com/Loop3D/map2loop/commit/b26d4dd9558792db939e4212b7b6bd2ada36a35b))
* adding docker file ([0acf2e4](https://github.com/Loop3D/map2loop/commit/0acf2e44bfc68b313d7be7c249583d74c8067b66))
* adding documentaion to map2loop ([2a30da4](https://github.com/Loop3D/map2loop/commit/2a30da44329d0adc0eb8026e861c7cc8fa23ffac))
* adding example for hamersley ([3e7688d](https://github.com/Loop3D/map2loop/commit/3e7688d1053225b12ab2ed3b897900f3bc4b4648))
* adding hjson file ([1d03505](https://github.com/Loop3D/map2loop/commit/1d03505debdf67ce2e7bc313e935173ed6868f38))
* adding m2l template ([b2f8ea0](https://github.com/Loop3D/map2loop/commit/b2f8ea0dd92455d5106a5d8260e56053ac5e2199))
* adding readme to docs ([e43968c](https://github.com/Loop3D/map2loop/commit/e43968ce25535fe6fb943b4d3fef8ce6e4b165dd))
* adding sphinx ([6e9b068](https://github.com/Loop3D/map2loop/commit/6e9b068dd46aed3dff40a19d4059f126cf265338))
* adding the remaining classes to docs ([c652d84](https://github.com/Loop3D/map2loop/commit/c652d842616091f00a4d38f291be1776b798e763))
* adding user guide ([47c8d9b](https://github.com/Loop3D/map2loop/commit/47c8d9b939d303ae7f68832c1550bd8da1735bca))
* adding user guide for map2loop  ([57c2263](https://github.com/Loop3D/map2loop/commit/57c22633d46d598f9a013e09984287a2a9953535))
* change docker run folder location ([54d5969](https://github.com/Loop3D/map2loop/commit/54d596949dae8815487037e802e7af4e39d362e6))
* changing map2loop-3 to map2loop ([1fe8c40](https://github.com/Loop3D/map2loop/commit/1fe8c4051d8e7c3489fa60f96d5245c3b0163467))
* finding the right json path may be the problem ([85ddcb3](https://github.com/Loop3D/map2loop/commit/85ddcb396cae2864f2f5b329ea0ce9513c03a71f))
* fix channel typo on readme ([feee724](https://github.com/Loop3D/map2loop/commit/feee724c9e7afdbdb388521229fe2b4641162fbe))
* fix link formatting for Brockman Syncline ([0dba31f](https://github.com/Loop3D/map2loop/commit/0dba31f50ee7f953e947f051c0b6d61be1945009))
* fix short underline ([9ecd5b4](https://github.com/Loop3D/map2loop/commit/9ecd5b40f646d4604ea78901e36080320530b804))
* fix the path for clut with PathLib too ([5c525af](https://github.com/Loop3D/map2loop/commit/5c525af1ecae4490ff8bd4615274329852f52922))
* Fixed indentation of code block ([5be5208](https://github.com/Loop3D/map2loop/commit/5be5208d0e3121ce3328dcc788309d7e81d7c199))
* fixing path to images ([6b311d6](https://github.com/Loop3D/map2loop/commit/6b311d6a6f3ce53e896c39e620823285b0ff79db))
* geographical --&gt; geological ([b368c19](https://github.com/Loop3D/map2loop/commit/b368c196e31484c023a21b2b61b6c6b2c8c9f2e0))
* ignore files ([9de740f](https://github.com/Loop3D/map2loop/commit/9de740f274a155b87e80a50776383a70b12e572a))
* keep consistent levels ([052de43](https://github.com/Loop3D/map2loop/commit/052de43496decf1c1b99f3b32990e06701f678af))
* make sure m2l is updated with local docker build ([dc4d8cc](https://github.com/Loop3D/map2loop/commit/dc4d8ccb3ee70fac9d9f3e489cd1510b9f945047))
* making code block work ([d0ed967](https://github.com/Loop3D/map2loop/commit/d0ed967eede6cab6373da3d92dfcdf7d004b678c))
* Move examples dir for doc generation ([241444f](https://github.com/Loop3D/map2loop/commit/241444fabf3bba69c8fdfbee279cef6d97b01f04))
* moving images ([a800455](https://github.com/Loop3D/map2loop/commit/a800455496e92c34853e82bef618ee541072358f))
* no need for lavavu ([efe4b9d](https://github.com/Loop3D/map2loop/commit/efe4b9d9263bb04437c33f882af6443052aa7f71))
* overwrite LPF ([62e75b0](https://github.com/Loop3D/map2loop/commit/62e75b0e2974e5af47646be64736c3dd384875d1))
* paths fixed to work with pathlib within docker ([5b5b503](https://github.com/Loop3D/map2loop/commit/5b5b503e1c69018ed02a44d34d97ef4d91f02380))
* remove conda channel from dependencies installation ([59fe775](https://github.com/Loop3D/map2loop/commit/59fe775f9c4307352e7e5320c9648c984927ffc2))
* remove python3.11 bind. m2model updated to 3.12 ([5e84484](https://github.com/Loop3D/map2loop/commit/5e84484420db6845159d2bb56a0ec5fdf99a5db0))
* removing LS import for config ([bcee039](https://github.com/Loop3D/map2loop/commit/bcee03998035ae92cb48dd2d02f921ce39bed91b))
* revert back to original ([f8d6d2c](https://github.com/Loop3D/map2loop/commit/f8d6d2caa79d54c64827d4ce3172c1173a95d096))
* revert docs to update m2l build locally ([e46c45e](https://github.com/Loop3D/map2loop/commit/e46c45e983f3481268a1fe294770cb5adb0a2651))
* test to see if python3.11 allows docs to build (map2model requirement) ([9d15efd](https://github.com/Loop3D/map2loop/commit/9d15efde25407512da1240cffcba7acb42528728))
* typo ([bdbe907](https://github.com/Loop3D/map2loop/commit/bdbe907bb16695865f42f74ad95f786a9976d0af))
* update docstring for thickness calc gamma class ([bfebedc](https://github.com/Loop3D/map2loop/commit/bfebedccb0802b56a54c60c9889ddc3f6802f7e4))
* update hamersley example to not use the legacy ([d85ffb9](https://github.com/Loop3D/map2loop/commit/d85ffb93d86159b21aa2e74948e64d8684cf485e))
* update installation steps for small clarifications ([c975f6f](https://github.com/Loop3D/map2loop/commit/c975f6f6a9b484d709c007bdaea214751e9fd649))
* update installation steps with channels ([76d7897](https://github.com/Loop3D/map2loop/commit/76d7897160f869ccb74945ae3487cce07d5c7bec))
* update link ([184e522](https://github.com/Loop3D/map2loop/commit/184e522f45d6acc0a4e9fd14fbd16503adbc8f7f))

## [3.0.6](https://github.com/Loop3D/map2loop/compare/3.0.5...3.0.6) (2024-01-25)


### Bug Fixes

* Remove deprecation ignore warnings and then deprecations ([79cd333](https://github.com/Loop3D/map2loop/commit/79cd33317625c96091690c9b97309be74cba2fdd))

## [3.0.5](https://github.com/Loop3D/map2loop/compare/3.0.4...3.0.5) (2024-01-25)


### Bug Fixes

* Typo in intrusion flag, ignoring intrusions for contact calculations, remove errant sorter from take_best option ([527de6c](https://github.com/Loop3D/map2loop/commit/527de6c339e56ead036ab5ed5d94d06d77673ae6))

## [3.0.4](https://github.com/Loop3D/map2loop-3/compare/3.0.3...3.0.4) (2024-01-11)


### Bug Fixes

* Issue with conda build test env not working, skipping ([2fffb69](https://github.com/Loop3D/map2loop-3/commit/2fffb696556f48f0620dfbac356971301fb15e89))

## [3.0.3](https://github.com/Loop3D/map2loop-3/compare/3.0.2...3.0.3) (2024-01-11)


### Bug Fixes

* Add owslib as dependency for testing ([0439da8](https://github.com/Loop3D/map2loop-3/commit/0439da869856b187defbed631c592a848b3684a1))

## [3.0.2](https://github.com/Loop3D/map2loop-3/compare/3.0.1...3.0.2) (2024-01-11)


### Bug Fixes

* Still forcing versioning for release-please ([08e647b](https://github.com/Loop3D/map2loop-3/commit/08e647bac1dded0c807bbb9a3571a741505bd488))
