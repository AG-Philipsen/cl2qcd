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

### Normal priority

 - [ ] :recycle: In every class, all members should be initialised to meaningful values in the initialisation list.
                 This is not done everywhere consistently and it should be checked and fixed.
 - [ ] :recycle: Around in the code base there are `boost::lexical_cast<std::string>` which should be all replaceable by `std::to_string`, avoiding then the external dependency.
                 In particular, there is no [necessity](https://stackoverflow.com/a/29399444) of using this cast at the moment.
 - [ ] :recycle: Tests should probably in general only accept some given command line options like --useCPU which are then passed in the `CMakeLists.txt` files.
                 However, no check is ever performed and, if the user/developer manually pass other options, they are sometimes silently ignored.
                 This should be probably improved.

### Low priority


----

## Past work
