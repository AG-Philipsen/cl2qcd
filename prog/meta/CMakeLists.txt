# Build specifications for the contents of the Meta package

add_executable(inputparameters_test inputparameters_test.cpp)
target_link_libraries(inputparameters_test
		optimal)
configure_file(test_input_1 test_input_1 COPYONLY)
configure_file(test_input_2 test_input_2 COPYONLY)
add_test(Inputparameters inputparameters_test)

add_library(meta inputparameters.cpp)
