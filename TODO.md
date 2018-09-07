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

### Normal priority


### Low priority


----

## Past work
