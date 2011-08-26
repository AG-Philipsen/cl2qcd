# Build configuration for OpenCL based HMC
#
# Copyright 2011 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>

option( PERFORM_PROFILING "Perform Profiling" ON )
if( PERFORM_PROFILING )
	add_definitions( -D_PROFILING_ )
endif( PERFORM_PROFILING )

# the executables should include the same files as the "original" ones (see main CMakeLists.txt)

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
	../exceptions.cpp
)
target_link_libraries( heatbath_benchmarks
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	einhard
)

add_executable( inverter_benchmarks
	inverter_benchmarks.cpp
	../host_geometry.cpp
	../host_operations_complex.cpp
	../host_operations_gaugefield.cpp
	../host_operations_su3matrix.cpp
	../host_operations_matrix.cpp
	../gaugefield.cpp
	../gaugefield_inversion.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../opencl.cpp
	../opencl_compiler.cpp
	../opencl_fermions.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../logger.cpp
	../exceptions.cpp
)
target_link_libraries( inverter_benchmarks
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	einhard
)

add_executable( hmc_benchmarks
	hmc_benchmarks.cpp
	../host_geometry.cpp
	../host_operations_complex.cpp
	../host_operations_gaugefield.cpp
	../host_operations_su3matrix.cpp
	../host_operations_matrix.cpp
	../gaugefield.cpp
	../gaugefield_inversion.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl.cpp
	../opencl_compiler.cpp
	../opencl_fermions.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../logger.cpp
	../gaugefield_hmc.cpp
	../opencl_hmc.cpp
	../exceptions.cpp
)
target_link_libraries( hmc_benchmarks
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	einhard
)

