# Build configuration for OpenCL based HMC
#
# Copyright 2011 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>

# the folowing are pure profiling benchmarks
option(PERFORM_PROFILING "Perform Profiling" ON)
if(PERFORM_PROFILING)
	add_executable(heatbath_benchmarks heatbath_benchmarks.cpp)
	target_link_libraries(heatbath_benchmarks optimal)

	add_executable(inverter_benchmarks inverter_benchmarks.cpp)
	target_link_libraries(inverter_benchmarks optimal)

	add_executable(dslash_benchmark dslash_benchmark.cpp)
	target_link_libraries(dslash_benchmark optimal)

	add_executable(hmc_benchmarks hmc_benchmarks.cpp)
	target_link_libraries(hmc_benchmarks optimal)

endif(PERFORM_PROFILING)
