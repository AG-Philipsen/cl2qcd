# Build configuration for OpenCL based HMC
#
# Copyright 2011 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>

option( PERFORM_PROFILING "Perform Profiling" ON )
if( PERFORM_PROFILING )
	add_definitions( -D_PROFILING_ )
endif( PERFORM_PROFILING )

add_executable( heatbath_benchmarks
	heatbath_benchmarks.cpp
	../host_geometry.cpp
	../host_operations_complex.cpp
	../host_operations_gaugefield.cpp
	../host_operations_su3matrix.cpp
	../host_operations_matrix.cpp
	../gaugefield.cpp
	../gaugefield_heatbath.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl.cpp
	../opencl_heatbath.cpp
	../opencl_compiler.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../logger.cpp
)
target_link_libraries( heatbath_benchmarks
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	einhard
)

