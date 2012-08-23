# Build specifications for the contents of the Meta package
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

#
# Collect all code in this part of the source tree
#

set_abs_paths(META_FILES
	inputparameters.cpp
)
add_library(meta ${META_FILES})

#
# Definition of tets
#

add_executable(inputparameters_test inputparameters_test.cpp)
target_link_libraries(inputparameters_test
		optimal)
configure_file(test_input_1 test_input_1 COPYONLY)
configure_file(test_input_2 test_input_2 COPYONLY)
add_test(meta::Inputparameters inputparameters_test)

