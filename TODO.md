# To-do list

Here we keep track of work that is either planned to be done in the long-term future or in between two consecutive releases.
In particular, one can refer to the [Past work section](#past-work) to get some insights about the changes that will be included in the next release.
However, such a section is not intended to be a complete list and, indeed, you should keep in mind that this file is more for the developers than for the users and some statements in the following may sound cryptic.
Rather refer to the [CHANGELOG](https://github.com/AG-Philipsen/cl2qcd/blob/master/CHANGELOG.md) to read about notable changes to the project in a more user-friendly form.

### Legend

* :new: Feature to be added
* :fire: Bug/inconsistency fixing
* :recycle: Refactoring
* :mag: Work about code testing
* :cl: Work implies device coding
* :question: Decision making
* :memo: Documentation

----

## Future work

### High priority

 - [ ] :fire: :recycle: At the moment the information written in the header of every LIME file is basically meaningless and actually misleading. Parameters written there are not correct and they should be consistent with those of the run that produced such a configuration. Moreover, the hmcversion should be called according to the algorithm and its value should be a string containing the code name and the real version (e.g. `CL2QCD_v1.0`).
 - [ ] :new: Require at least version `1.60` of `boost` in order to freely use new unit tests features about [command line options](https://www.boost.org/doc/libs/1_60_0/libs/test/doc/html/boost_test/change_log.html).
             The code base has been already prepared for that and it should be matter of few changes in the main `CMakeLists.txt` file.
             However, tests should be run to check that everything is really fine.
 - [ ] :new: :question: [It seems](https://stackoverflow.com/a/42124857) that for more recent versions of `boost`, also recent versions of **CMake** are needed to correctly detect dependencies in `boost`.
                        It might be good to implement some kind of check in the main `CMakeLists.txt` file (note that for `boost 1.62` also `CMake 3.7` [is needed](https://stackoverflow.com/a/40929234)).
 - [ ] :mag: Write meaningful tests for executables, avoiding, if possible, any non-deterministic behaviour.
             Adjust the `CMakeLists.txt` file, using our `add_unit_test` **CMake** macro.
 - [ ] :mag: Add tests to the hardware/lattices folder, where part of the functionality that was in physics/lattices has been moved to.
             Remove consistently tests from physics/lattices which will be moved.
 - [ ] :mag: Refactor all tests in `physics`, relying to analytically determined values (or on reference values calculated in Mathematica).
             In this refactoring, mockups for the parameters should be created.
             Afterwards, the dependence on `meta` of `system` can be finally removed (i.e. remove `System(meta::Inputparameters& parameters)` constructor).
 - [ ] :mag: :question: The `--nDevices=1` option is not specified in all tests.
                        When not specified, the test runs on all available devices.
                        Should this be the case, or is it better to act differently?
 - [ ] :mag: :question: Most of the tests at the physics level contain the initialisation of a C-array of `char*` to be passed as kind of `argv` to the `Inputparameters` class.
                        This implies that the size of such an array is hard-coded as `argc` and this aspect is definitely not nice to have.
                        It is in general easier to build a `std::vector<std::string>` and then convert it to `std::vector<const char*>`, whose methods `size()` and `data()` give `argc` and `argv`, respectively and dynamically.
                        In the future, maybe, once refactored all the tests in physics with mockups for the parameters, this `argc/argv` approach might not be needed.
 - [ ] :recycle: :fire: At the level of `hardware/code`, all the classes with `cl_kernels` as private members should initialize them to `0` in the constructor initializer list.
                        Then, in the `fill_kernels` function, some of them should be built depending on the parameters.
                        In this way, those not built would automatically be set to `0`.
                        At the moment, this is not consistently done in all modules and often some `cl_kernels` are set to `0` in the `fill_kernels` function.
 - [ ] :recycle: :cl: Names of functions, variables, etc. in the `ocl_kernel` directory are definitely bad and it is very hard to read the code because of that.
                      Some work should be done to improve the situation, but tests should be run carefully to check that nothing was broken.

### Normal priority

 - [ ] :fire: :question: The compilation with the **CMake** variable `LAZY_HALO_UPDATES` switched off has been fixed in commit [88affa4](https://github.com/AG-Philipsen/cl2qcd/commit/88affa4ff4f68eddfd27c93dd951471511d376e6), but in general one should guarantee the correctness of multi-GPU support independently from the state of the **CMake** variables `LAZY_HALO_UPDATES` and `ASYNC_HALO_UPDATES`.
                         Specific tests, if possible, should be developed and this could be a good opportunity to make arbitrary the direction in which the lattice is split.
 - [ ] :recycle: In every class, all members should be initialised to meaningful values in the initialisation list.
                 This is not done everywhere consistently and it should be checked and fixed.
 - [ ] :recycle: Around in the code base there are `boost::lexical_cast<std::string>` which should be all replaceable by `std::to_string`, avoiding then the external dependency.
                 In particular, there is no [necessity](https://stackoverflow.com/a/29399444) of using this cast at the moment.
 - [ ] :recycle: Tests should probably in general only accept some given command line options like `--useCPU` which are then passed in the `CMakeLists.txt` files.
                 However, no check is ever performed and, if the user/developer manually pass other options, they are sometimes silently ignored.
                 This should be probably improved.

### Low priority


----

## Past work
