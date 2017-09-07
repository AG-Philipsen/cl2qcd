# Build configuration for OpenCL based HMC
#
# Copyright 2011-2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
# Copyright 2011-2013 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
# Copyright 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
#
# This file is part of CL2QCD.
#
# CL2QCD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CL2QCD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required( VERSION 2.8.3 )
project( OpTiMaL )

# we have some custom CMake stuff
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

include(Macros)

# if a user doesn't specify anything he probably isn't interested in debug information
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Have the user asked for tests
enable_testing()
# And enable Dashboard
include(CTest)

# make the compiler bitchy
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER MATCHES "/clang\\+\\+$")
	set(CMAKE_CXX_FLAGS           "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Wdisabled-optimization -std=c++11")
	set(CMAKE_CXX_FLAGS_RELEASE   "-O4") # LTO on clang, gcc will fall back to O3
	set(CMAKE_CXX_FLAGS_DEBUG     "-g -O0 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
	set(CMAKE_CXX_FLAGS_DEBUGFULL "-g3 -fno-linline")
endif(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER MATCHES "/clang\\+\\+$")
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_C_COMPILER MATCHES "/clang\\+\\+$")
	set(CMAKE_C_FLAGS             "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic -Wc++-compat -Wdisabled-optimization -std=c99")
	set(CMAKE_C_FLAGS_RELEASE     "-O4") # LTO on clang, gcc will fall back to O3
	set(CMAKE_C_FLAGS_DEBUG       "-g -O0 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
	set(CMAKE_C_FLAGS_DEBUGFULL   "-g3 -fno-linline")
endif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_C_COMPILER MATCHES "/clang\\+\\+$")

# add endianess information
include(TestBigEndian)
TEST_BIG_ENDIAN(BIG_ENDIAN)
message("Big Endian: ${BIG_ENDIAN}")
if(BIG_ENDIAN)
	message("Big endian arch..")
	add_definitions( "BIG_ENDIAN_ARCH=1" )
else()
	message("Little endian arch..")
endif()

# make sure we can compile agains opencl
find_package( OpenCL REQUIRED )
include_directories( SYSTEM ${OPENCL_INCLUDE_DIRS} )
# work around OpenCL implementations not shipping c++ bindings
if( NOT OPENCL_HAS_CPP_BINDINGS )
	message(STATUS "No OpenCL CPP bindings found, using hack" )
	include_directories( SYSTEM hacks/cl_hpp )
elseif( OPENCL_HAS_CPP_BINDINGS )
	message (STATUS "OpenCL CPP bindings found" )
endif( NOT OPENCL_HAS_CPP_BINDINGS )
# make sure we can uses OpenCL 1.1 features on AMD APP 2.7
add_definitions(-DCL_USE_DEPRECATED_OPENCL_1_1_APIS)

find_package( LIME REQUIRED )
include_directories( SYSTEM ${LIME_INCLUDE_DIR} )

find_package( GMP REQUIRED )
include_directories( SYSTEM ${GMP_INCLUDE_DIR} )

find_package( MPFR REQUIRED )
include_directories( SYSTEM ${MPFR_INCLUDE_DIR} )

find_package( LibXml2 REQUIRED )
include_directories( SYSTEM ${LIBXML2_INCLUDE_DIR} )
if( LIBXML2_DEFINITIONS )
	add_definitions( ${LIBXML2_DEFINITIONS} )
endif( LIBXML2_DEFINITIONS )

# require 1.48 as boost-regex of 1.46 can't be compiled with c++11-enabled g++-4.7
find_package(Boost 1.48.0 REQUIRED COMPONENTS regex filesystem system program_options unit_test_framework)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# Define where to find the OpenCL kernels
# TODO this should also include the installation place
# maybe even patch this during install step
add_definitions( -DSOURCEDIR="${CMAKE_SOURCE_DIR}" )

# Add code version to the code
# we need to include the binary dir as the custom target will create a header
# there
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_custom_target(git_commit_id_header ALL DEPENDS git_commit_id )
add_custom_command(OUTPUT git_commit_id
                   COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/GenerateGitVersion.cmake)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/gitcommitid.h
                            PROPERTIES GENERATED TRUE
                            HEADER_FILE_ONLY TRUE)

# Allow to choose the default log level
set(LOG_LEVEL "ALL" CACHE STRING "The minimum log level for output to be written. Options are: ALL TRACE DEBUG INFO WARN ERROR FATAL OFF. Note that some tests will break if INFO messages are not logged.")
string(REGEX MATCH "^(ALL)|(TRACE)|(DEBUG)|(INFO)|(WARN)|(ERROR)|(FATAL)|(OFF)$" TMP "${LOG_LEVEL}")
if(TMP STREQUAL "")
	message(SEND_ERROR "LOG_LEVEL must be one of ALL, TRACE, DEBUG, INFO, WARN, ERROR, FATAL or OFF. \"${LOG_LEVEL}\" is not a valid value.")
endif(TMP STREQUAL "")
add_definitions(-DLOG_LEVEL=${LOG_LEVEL})
# Don't compile in debug output if compiling for release
if( CMAKE_BUILD_TYPE STREQUAL "Release" )
	add_definitions( -DNDEBUG )
endif( CMAKE_BUILD_TYPE STREQUAL "Release" )

# some options for various variants
option( USE_DOUBLE_PRECISION "Use double precision arithmetics" ON )
if( USE_DOUBLE_PRECISION )
	add_definitions( -D_USEDOUBLEPREC_ )
endif( USE_DOUBLE_PRECISION )

# allow to choose a PRNG
set(PRNG "Ranlux" CACHE STRING "The pseudo random number generator to use. Options are: Ranlux (default)")
if(PRNG STREQUAL "Ranlux")
	add_definitions(-DUSE_PRNG_RANLUX)
	set(PRNG_LIBRARIES ranlux)
else()
	message(SEND_ERROR "PRNG must be Ranlux (default). \"${PRNG}\" is not a valid value.")
endif()

# whether to use lazy halo update
option(LAZY_HALO_UPDATES "Delay halo updates until the halo is required." ON)
if(LAZY_HALO_UPDATES)
	add_definitions(-DLAZY_HALO_UPDATES)
endif()

# whether to enable asynchroneous halo updates
option(ASYNC_HALO_UPDATES "Perform device-to-device transfers asynchroneous to kernel execution." ON)
if(ASYNC_HALO_UPDATES)
	add_definitions(-DASYNC_HALO_UPDATES)
endif()

# default random number seed
configure_file( rand_seeds rand_seeds COPYONLY )

# which directories contain documentation relevant stuff
# modules will add their files to this list
set_property(GLOBAL PROPERTY DOC_SOURCE_DIRS
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/klepsydra
	${CMAKE_SOURCE_DIR}/einhard
	${CMAKE_SOURCE_DIR}/crypto
)

# Definition of the core application
set_abs_paths(OPTIMAL_FILES
)
set(OPTIMAL_USED_LIBS
	${OPENCL_LIBRARIES}
#${GMP_LIBRARIES} automatically pulled by mpfr, relying on the implicit pull avoids libary version conflicts
	${MPFR_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	${PRNG_LIBRARIES}
	crypto
	executables
)

add_library(optimal ${OPTIMAL_FILES})
# add all modules
add_modules(optimal
	meta
	hardware
	geometry
	host_functionality
	ildg_io
	physics
	einhard
	executables
)
target_link_libraries(optimal ${OPTIMAL_USED_LIBS})
set_target_properties(optimal PROPERTIES LINK_INTERFACE_MULTIPLICITY 4)

# The actual executables (frontends)
add_executable(gaugeobservables executables/gaugeobservablesMain.cpp)
target_link_libraries(gaugeobservables optimal)

add_executable(su3heatbath executables/su3heatbathMain.cpp)
target_link_libraries(su3heatbath optimal)

add_executable(inverter executables/inverterMain.cpp)
target_link_libraries(inverter optimal)

add_executable(hmc executables/hmcMain.cpp)
target_link_libraries(hmc optimal)

add_executable(rhmc executables/rhmcMain.cpp)
target_link_libraries(rhmc optimal)

# Utilities
add_library(klepsydra
	klepsydra/timer.cpp
	klepsydra/monotonic.cpp)
if(NOT CMAKE_HOST_APPLE)
	target_link_libraries(klepsydra rt)
endif(NOT CMAKE_HOST_APPLE)

# Boost test
add_subdirectory(boost)

# Random number generation lib
add_subdirectory(ranlux)
add_subdirectory(ranluxcl)

# Crypto
add_subdirectory(crypto)

add_subdirectory(interfaceImplementations)

option(BUILD_MICROBENCHMARKS "Build micro benchmarks" OFF)
if(BUILD_MICROBENCHMARKS)
	add_subdirectory(microbenchmarks)
endif(BUILD_MICROBENCHMARKS)

# Benchmarking
option(BUILD_BENCHMARKS "Build benchmarks" OFF)
if(BUILD_BENCHMARKS)
	add_subdirectory(benchmarks)
endif(BUILD_BENCHMARKS)

# Documentation
#include(FindDoxygen)
find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
    add_subdirectory(doc)
else(DOXYGEN_FOUND)
    message(STATUS "WARNING: Failed to find Doxygen - won't build documentation")
endif(DOXYGEN_FOUND)

