# Build configuration for OpenCL based HMC
#
# Copyright 2011 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>

# the folowing are pure profiling benchmarks
option(PERFORM_PROFILING "Perform Profiling" ON)
if(PERFORM_PROFILING)
	add_definitions( -D_PROFILING_ )

	# the executables should include the same files as the "original" ones (see main CMakeLists.txt)
	add_library(optimal_profile ${OPTIMAL_FILES})
	target_link_libraries(optimal_profile ${OPTIMAL_USED_LIBS})
	add_dependencies(optimal_profile git_commit_id_header)

	add_executable(heatbath_benchmarks heatbath_benchmarks.cpp)
	target_link_libraries(heatbath_benchmarks optimal_profile)

	add_executable(inverter_benchmarks inverter_benchmarks.cpp)
	target_link_libraries(inverter_benchmarks optimal_profile)

	add_executable(dslash_benchmark dslash_benchmark.cpp)
	target_link_libraries(dslash_benchmark optimal_profile)

	add_executable(hmc_benchmarks hmc_benchmarks.cpp)
	target_link_libraries(hmc_benchmarks optimal_profile)

endif(PERFORM_PROFILING)
