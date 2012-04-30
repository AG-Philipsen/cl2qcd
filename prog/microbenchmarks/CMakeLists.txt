# Build configuration for some microbenchmarks for OpenCL based HMC
#
# Copyright 2011 Matthias Bach <bach@compeng.uni-frankfurt.de>

find_package(Boost 1.33.0 REQUIRED COMPONENTS program_options)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

add_definitions( -D_PROFILING_ )

add_executable( bandwidth
	bandwidth.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../logger.cpp
	../exceptions.cpp
)
target_link_libraries( bandwidth
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard
	${PRNG_LIBRARIES}
)
add_dependencies(bandwidth git_commit_id_header)
