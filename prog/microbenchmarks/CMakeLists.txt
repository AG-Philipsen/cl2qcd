# Build configuration for some microbenchmarks for OpenCL based HMC
#
# Copyright 2011 Matthias Bach <bach@compeng.uni-frankfurt.de>

add_executable(bandwidth bandwidth.cpp)
target_link_libraries( bandwidth
	optimal
)
