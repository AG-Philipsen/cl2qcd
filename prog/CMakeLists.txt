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
## USE_HMC is now moved to extra executable
#option( USE_HMC "Update with hybrid monte carlo algorithm" OFF)
#if (USE_HMC)
#	add_definitions( -D_USEHMC_ )
#endif( USE_HMC)
##LZ USE_FERMIONS is now a inputparameter option
#option( USE_FERMIONS "Use fermions" OFF )
#if( USE_FERMIONS )
#	add_definitions( -D_FERMIONS_ )
#endif( USE_FERMIONS )
#if( USE_FERMIONS )
#	option( USE_TWISTEDMASS "Use twisted mass fermions" ON)
#	if ( USE_TWISTEDMASS )
#		add_definitions( -D_TWISTEDMASS_ )
#	endif ( USE_TWISTEDMASS )
#endif( USE_FERMIONS )
#if( USE_FERMIONS AND !USE_TWISTEDMASS)
#	option( USE_CLOVER "Use clover fermions" OFF)
#		if ( USE_CLOVER )
#			add_definitions( -D_CLOVER_ )
#		endif (USE_CLOVER)
#endif (  )
## We do not want to use TESTING anymoe
#option( TESTING "Run optional test code" OFF )
#if( TESTING )
#	add_definitions( -D_TESTING_ )
#endif( TESTING )
option( USE_GPU "Use GPU" OFF )
if( USE_GPU )
	add_definitions( -D_USEGPU_ )
endif( USE_GPU )
## We do not want to use BENCHMARKS anymore
#option( PERFORM_BENCHMARKS "Perform Benchmarks" OFF )
#if( PERFORM_BENCHMARKS )
#	add_definitions( -D_PERFORM_BENCHMARKS_ )
#endif( PERFORM_BENCHMARKS )
option( USE_CHEM_POT_REAL "Apply real chemical potential" OFF )
if( USE_CHEM_POT_REAL )
	add_definitions( -D_CP_REAL_ )
endif( USE_CHEM_POT_REAL )
option( USE_CHEM_POT_IMAG "Apply imaginary chemical potential" OFF )
if( USE_CHEM_POT_IMAG )
	add_definitions( -D_CP_IMAG_ )
endif( USE_CHEM_POT_IMAG )
option( USE_NPBC_T "Apply non-periodic BC (temporal)" OFF )
if( USE_NPBC_T )
	add_definitions( -D_NPBC_T_ )
endif( USE_NPBC_T )
option( USE_NPBC_S "Apply non-periodic BC (spatial)" OFF )
if( USE_NPBC_S )
	add_definitions( -D_NPBC_S_ )
endif( USE_NPBC_S )
set( EXPIQ 2 CACHE STRING "Define order to which exp(iQ) is exact (2,3,ALL,MORNINGSTAR)")
if(EXPIQ MATCHES "2")
        add_definitions( -D_EXPONENTIATE_ALGEBRA_ORDER_2_ )
elseif(EXPIQ MATCHES "3")
        add_definitions( -D_EXPONENTIATE_ALGEBRA_ORDER_3_ )
elseif(EXPIQ MATCHES "ALL")
        add_definitions( -D_EXPONENTIATE_ALGEBRA_ORDER_ALL_ORDERS_ )
elseif(EXPIQ MATCHES "MORNINGSTAR")
        add_definitions( -D_USE_MORNINGSTAR_PEARDON_ )
endif (EXPIQ MATCHES "2")
option( USE_SMEARING "Use smearing of gauge-links" OFF)
if(USE_SMEARING)
	add_definitions( -D_USESMEARING_ )
endif(USE_SMEARING)
option( USE_EOPREC "Use even-odd-preconditioning" OFF)
if(USE_EOPREC)
	add_definitions( -D_USEEOPREC_ )
endif(USE_EOPREC)

# lattice dimensions (see feature #127)
set( NTIME 8 CACHE STRING "Number of lattice sites in time direction (t)" )
add_definitions( -DNTIME=${NTIME} )
set( NSPACE 4 CACHE STRING "Number of lattice sites in space direction (x,y,z each)" )
add_definitions( -DNSPACE=${NSPACE} )



# default random number seed
configure_file( rand_seeds rand_seeds COPYONLY )

add_executable( hmc
	hmc.cpp
	host_geometry.cpp
	host_operations_complex.cpp
	host_operations_matrix.cpp
	host_operations_su3matrix.cpp
	host_operations_gaugefield.cpp
	host_operations_gaugemomentum.cpp
	host_operations_spinor.cpp
	host_operations_spinorfield.cpp
	host_operations_fermionmatrix.cpp
#	host_gaugeobservables.cpp
# 	host_testing.cpp
	inputparameters.cpp
	host_readgauge.cpp
	host_random.cpp
	host_update_heatbath.cpp
	#opencl.cpp
	host_use_timer.cpp
#	host_gaugefieldoperations.cpp
	host_writegaugefield.cpp
	host_solver.cpp
	host_fermionobservables.cpp
	host_hmc.cpp
)
target_link_libraries( hmc
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	klepsydra
)
add_executable( heatbath
	heatbath.cpp
	host_geometry.cpp
	host_operations_complex.cpp
	host_operations_gaugefield.cpp
	host_operations_su3matrix.cpp
	host_operations_matrix.cpp	
	gaugefield.cpp
	inputparameters.cpp
	host_readgauge.cpp
	host_random.cpp
	opencl.cpp
	host_use_timer.cpp
	host_writegaugefield.cpp
	logger.cpp
)
target_link_libraries( heatbath
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
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
	inputparameters.cpp
	host_readgauge.cpp
	host_random.cpp
	opencl.cpp
	opencl_k.cpp
	host_use_timer.cpp
	host_writegaugefield.cpp
	logger.cpp
)
target_link_libraries( tk_kappa
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
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
	host_operations_spinor.cpp
	host_operations_spinorfield.cpp
	host_operations_fermionmatrix.cpp
	gaugefield.cpp
	gaugefield_inversion.cpp
	inputparameters.cpp
	host_readgauge.cpp
	host_random.cpp
	opencl.cpp
	opencl_fermions.cpp
	host_use_timer.cpp
	host_writegaugefield.cpp
	host_solver.cpp
	logger.cpp
)
target_link_libraries( inverter
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
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
