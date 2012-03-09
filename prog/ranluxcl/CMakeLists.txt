add_executable(ranluxcltest ranluxcltest.cpp)
target_link_libraries( ranluxcltest
	${OPENCL_LIBRARIES}
)
