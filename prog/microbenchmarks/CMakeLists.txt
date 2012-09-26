# Build configuration for some microbenchmarks for OpenCL based HMC
#
# Copyright 2011 Matthias Bach <bach@compeng.uni-frankfurt.de>

if(PERFORM_PROFILING)

	add_executable(bandwidth bandwidth.cpp)
	target_link_libraries( bandwidth
		optimal
	)

endif(PERFORM_PROFILING)
