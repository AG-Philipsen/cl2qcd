# Build configuration for some microbenchmarks for OpenCL based HMC
#
# Copyright 2011 Matthias Bach <bach@compeng.uni-frankfurt.de>

if(PERFORM_PROFILING)

	find_package(Boost 1.33.0 REQUIRED COMPONENTS program_options)
	include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

	add_executable(bandwidth bandwidth.cpp)
	target_link_libraries( bandwidth
		optimal_profile
		${Boost_LIBRARIES}
	)

endif(PERFORM_PROFILING)
