# Build configuration for OpenCL based HMC
#
# Copyright 2011 Matthias Bach <bach@compeng.uni-frankfurt.de>

cmake_minimum_required( VERSION 2.6.4 )
project( OpTiMaL )

# we have some custom CMake stuff
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# if a user doesn't specify anything he probably isn't interested in debug information
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# make the compiler bitchy
if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS           "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Wdisabled-optimization -std=c++0x")
	set(CMAKE_CXX_FLAGS_RELEASE   "-O3")
	set(CMAKE_CXX_FLAGS_DEBUG     "-g -O0 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
	set(CMAKE_CXX_FLAGS_DEBUGFULL "-g3 -fno-linline")
endif(CMAKE_COMPILER_IS_GNUCXX)
if(CMAKE_COMPILER_IS_GNUCC)
	set(CMAKE_C_FLAGS             "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic -Wc++-compat -Wdisabled-optimization -std=c99")
	set(CMAKE_C_FLAGS_RELEASE     "-O3")
	set(CMAKE_C_FLAGS_DEBUG       "-g -O0 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
	set(CMAKE_C_FLAGS_DEBUGFULL   "-g3 -fno-linline")
endif(CMAKE_COMPILER_IS_GNUCC)

# make sure we can compile agains opencl
find_package( OpenCL REQUIRED )
include_directories( SYSTEM ${OPENCL_INCLUDE_DIRS} )

find_package( LIME REQUIRED )
include_directories( SYSTEM ${LIME_INCLUDE_DIR} )

find_package( LibXml2 REQUIRED )
include_directories( SYSTEM ${LIBXML2_INCLUDE_DIR} )
if( LIBXML2_DEFINITIONS )
	add_definitions( ${LIBXML2_DEFINITIONS} )
endif( LIBXML2_DEFINITIONS )

find_package(Boost 1.33.0 REQUIRED COMPONENTS regex)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# Define where to find the OpenCL kernels
# TODO this should also include the installation place
# maybe even patch this during install step
add_definitions( -DSOURCEDIR="${CMAKE_SOURCE_DIR}" )

# Don't compile in debug output if compiling for release
if( CMAKE_BUILD_TYPE STREQUAL "Release" )
	add_definitions( -DNDEBUG )
endif( CMAKE_BUILD_TYPE STREQUAL "Release" )

# some options for various variants
option( USE_DOUBLE_PRECISION "Use double precision arithmetics" ON )
if( USE_DOUBLE_PRECISION )
	add_definitions( -D_USEDOUBLEPREC_ )
endif( USE_DOUBLE_PRECISION )
option( USE_RECONSTRUCT_TWELVE "Use compactified storage format (reconstruct twelve)" OFF )
if( USE_RECONSTRUCT_TWELVE )
	add_definitions( -D_RECONSTRUCT_TWELVE_ )
endif( USE_RECONSTRUCT_TWELVE )

# lattice dimensions (see feature #127)
set( NTIME 8 CACHE STRING "Number of lattice sites in time direction (t)" )
add_definitions( -DNTIME=${NTIME} )
set( NSPACE 4 CACHE STRING "Number of lattice sites in space direction (x,y,z each)" )
add_definitions( -DNSPACE=${NSPACE} )

# default random number seed
configure_file( rand_seeds rand_seeds COPYONLY )

add_executable( heatbath
	heatbath.cpp
	host_geometry.cpp
	host_operations_complex.cpp
	host_operations_gaugefield.cpp
	host_operations_su3matrix.cpp
	host_operations_matrix.cpp
	gaugefield.cpp
	gaugefield_heatbath.cpp
	inputparameters.cpp
	host_readgauge.cpp
	host_random.cpp
	opencl.cpp
	opencl_heatbath.cpp
	opencl_compiler.cpp
	host_use_timer.cpp
	host_writegaugefield.cpp
	logger.cpp
	exceptions.cpp
)
target_link_libraries( heatbath
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	einhard
)
add_executable( tk_kappa
	tk_kappa.cpp
	host_geometry.cpp
	host_operations_complex.cpp
	host_operations_gaugefield.cpp
	host_operations_su3matrix.cpp
	host_operations_matrix.cpp
	gaugefield.cpp
	gaugefield_k.cpp
	gaugefield_heatbath.cpp
	inputparameters.cpp
	host_readgauge.cpp
	host_random.cpp
	opencl.cpp
	opencl_compiler.cpp
	opencl_k.cpp
	opencl_heatbath.cpp
	host_use_timer.cpp
	host_writegaugefield.cpp
	logger.cpp
	exceptions.cpp
)
target_link_libraries( tk_kappa
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	einhard
)

add_executable( tk_kappa_hybrid
	tk_kappa_hybrid.cpp
	host_geometry.cpp
	host_operations_complex.cpp
	host_operations_gaugefield.cpp
	host_operations_su3matrix.cpp
	host_operations_matrix.cpp
	gaugefield_hybrid.cpp
	gaugefield_heatbath_kappa.cpp
	inputparameters.cpp
	host_readgauge.cpp
	host_random.cpp
	host_use_timer.cpp
	host_writegaugefield.cpp
	logger.cpp
	exceptions.cpp
	opencl_module.cpp
	opencl_compiler.cpp
)
target_link_libraries( tk_kappa_hybrid
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	einhard
)

add_executable( inverter
	inverter.cpp
	host_geometry.cpp
	host_operations_complex.cpp
	host_operations_gaugefield.cpp
	host_operations_su3matrix.cpp
	host_operations_matrix.cpp
	gaugefield.cpp
	gaugefield_inversion.cpp
	inputparameters.cpp
	host_readgauge.cpp
	opencl.cpp
	opencl_compiler.cpp
	opencl_fermions.cpp
	host_use_timer.cpp
	host_writegaugefield.cpp
	logger.cpp
	exceptions.cpp
)
target_link_libraries( inverter
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	einhard
)

add_executable( hmc
	hmc.cpp
	host_geometry.cpp
	host_operations_complex.cpp
	host_operations_gaugefield.cpp
	host_operations_su3matrix.cpp
	host_operations_matrix.cpp
	gaugefield.cpp
	gaugefield_inversion.cpp
	inputparameters.cpp
	host_readgauge.cpp
	host_random.cpp
	opencl.cpp
	opencl_compiler.cpp
	opencl_fermions.cpp
	host_use_timer.cpp
	host_writegaugefield.cpp
	logger.cpp
	gaugefield_hmc.cpp
	opencl_hmc.cpp
	exceptions.cpp
)
target_link_libraries( hmc
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	einhard
)

add_library( klepsydra
	klepsydra/timer.cpp
	klepsydra/monotonic.cpp )
if(NOT CMAKE_HOST_APPLE)
	target_link_libraries( klepsydra rt )
endif(NOT CMAKE_HOST_APPLE)

add_library( einhard
	einhard/einhard.cpp )

# Documentation
#include(FindDoxygen)
find_package( Doxygen QUIET )
if(DOXYGEN_FOUND)
    add_subdirectory(doc)
else(DOXYGEN_FOUND)
    message(STATUS "WARNING: Failed to find Doxygen - won't build documentation")
endif(DOXYGEN_FOUND)

option(BUILD_MICROBENCHMARKS "Build micro benchmarks" OFF)
if(BUILD_MICROBENCHMARKS)
	add_subdirectory(microbenchmarks)
endif(BUILD_MICROBENCHMARKS)

# Take care of unit tests
enable_testing()
add_subdirectory(tests)

# Benchmarking
option(BUILD_BENCHMARKS "Build benchmarks" OFF)
if(BUILD_BENCHMARKS)
	add_subdirectory(benchmarks)
endif(BUILD_BENCHMARKS)