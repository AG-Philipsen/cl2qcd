# Build configuration for OpenCL based HMC
#
# Copyright 2011 Matthias Bach <bach@compeng.uni-frankfurt.de>

cmake_minimum_required( VERSION 2.6.4 )
project( clhmc )

# we have some custom CMake stuff
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# make the compiler bitchy
if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS           "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Wdisabled-optimization -std=c++0x")
	set(CMAKE_CXX_FLAGS_RELEASE   "-O3")
	set(CMAKE_CXX_FLAGS_DEBUG     "-g -O0 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
	set(CMAKE_CXX_FLAGS_DEBUGFULL "-g3 -fno-linline")
endif(CMAKE_COMPILER_IS_GNUCXX)
if(CMAKE_COMPILER_IS_GNUC)
	set(CMAKE_C_FLAGS             "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic -Wc++-compat -Wdisabled-optimization -std=c99")
	set(CMAKE_C_FLAGS_RELEASE     "-O3")
	set(CMAKE_C_FLAGS_DEBUG       "-g -O0 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
	set(CMAKE_C_FLAGS_DEBUGFULL   "-g3 -fno-linline")
endif(CMAKE_COMPILER_IS_GNUC)

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

# Define where to find the OpenCL kernels
# TODO this should also include the installation place
# maybe even patch this during install step
add_definitions( -DSOURCEDIR="${CMAKE_SOURCE_DIR}" )

# some options for various variants
option( USE_DOUBLE_PRECISION "Use double precision arithmetics" ON )
if( USE_DOUBLE_PRECISION )
	add_definitions( -D_USEDOUBLEPREC_ )
endif( USE_DOUBLE_PRECISION )
option( USE_RECONSTRUCT_TWELVE "Use compactified storage format (reconstruct twelve)" OFF )
if( USE_RECONSTRUCT_TWELVE )
	add_definitions( -D_RECONSTRUCT_TWELVE_ )
endif( USE_RECONSTRUCT_TWELVE )
option( USE_FERMIONS "Use fermions" OFF )
if( USE_FERMIONS )
	add_definitions( -D_FERMIONS )
endif( USE_FERMIONS )
option( TESTING "Run optional test code" OFF )
if( TESTING )
	add_definitions( -D_TESTING )
endif( TESTING )

add_executable( hmc
	hmc.cpp
	host_operations_complex.cpp
	host_operations_gaugefield.cpp
	host_operations_spinor.cpp
	host_operations_spinorfield.cpp
	host_operations_fermionmatrix.cpp
	host_gaugeobservables.cpp
	host_testing.cpp
	host_geometry.cpp
	host_input.cpp
	host_readgauge.cpp
	host_random.cpp
	host_update_heatbath.cpp
	host_timer.cpp
	opencl.cpp
	host_use_timer.cpp
	host_gaugefieldoperations.cpp
	host_writegaugefield.cpp
	host_solver.cpp
	host_fermionobservables.cpp
)
target_link_libraries( hmc
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
)

