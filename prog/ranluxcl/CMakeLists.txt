add_executable(ranluxcltest ranluxcltest.cpp)
target_link_libraries( ranluxcltest
	${OPENCL_LIBRARIES}
)
add_test(NAME Ranluxcl_CPU COMMAND ranluxcltest 1 0 1 1)
add_test(NAME Ranluxcl_GPU COMMAND ranluxcltest 1 1 1 1)
set_tests_properties(Ranluxcl_CPU Ranluxcl_GPU PROPERTIES PASS_REGULAR_EXPRESSION "Correctness check PASSED!")
