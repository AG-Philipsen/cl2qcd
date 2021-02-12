# Cmake module
#
# Copyright (c) 2010-2011 Matthias Bach
# Copyright (c) 2021 Alessandro Sciarra
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
#

# Try to find OpenCL
#
# This module tries to find an OpenCL implementation on your system. It was written to
# support AMD/ATI, but for Apple and NVIDIA implementations might work, too.
#
# To set manually the paths, define these environment variables:
#   OpenCL_INCPATH    - Include path (e.g. OpenCL_INCPATH=/opt/cuda/4.0/cuda/include)
#   OpenCL_LIBPATH    - Library path (e.g. OpenCL_LIBPATH=/usr/lib64/nvidia)
#
# Once done this module will define
#   OpenCL_FOUND         - system has OpenCL
#   OpenCL_INCLUDE_DIRS  - the OpenCL include directories
#   OpenCL_LIBRARIES     - the OpenCL libraries
#
# WIN32 should work, but is untested

message(STATUS "Using customized FindOpenCL module")
find_package(PackageHandleStandardArgs)

message(WARNING "Not looking for OpenCL version, assuming version 1.0")
set(OpenCL_VERSION_STRING "1.0")
set(OpenCL_VERSION_MAJOR 1)
set(OpenCL_VERSION_MINOR 0)

if(APPLE)

    find_library(OpenCL_LIBRARIES OpenCL DOC "OpenCL lib for OSX")
    find_path(OpenCL_INCLUDE_DIRS OpenCL/cl.h DOC "Include for OpenCL on OSX")
    find_path(_OPENCL_CPP_INCLUDE_DIRS OpenCL/cl.hpp DOC "Include for OpenCL CPP bindings on OSX")

else(APPLE)

    if(WIN32)

        find_path(OpenCL_INCLUDE_DIRS CL/cl.h)
        find_path(_OPENCL_CPP_INCLUDE_DIRS CL/cl.hpp)

        # The AMD SDK currently installs both x86 and x86_64 libraries
        # This is only a hack to find out architecture
        if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "AMD64")
            set(_OPENCL_LIB_DIR "$ENV{ATISTREAMSDKROOT}/lib/x86_64")
        else(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "AMD64")
            set(_OPENCL_LIB_DIR "$ENV{ATISTREAMSDKROOT}/lib/x86")
        endif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "AMD64")
        find_library(OpenCL_LIBRARIES OpenCL.lib PATHS ${_OPENCL_LIB_DIR} ENV OpenCL_LIBPATH)

        get_filename_component(_OPENCL_INC_CAND ${_OPENCL_LIB_DIR}/../../include ABSOLUTE)

        # On Win32 search relative to the library
        find_path(OpenCL_INCLUDE_DIRS CL/cl.h PATHS "${_OPENCL_INC_CAND}" ENV OpenCL_INCPATH)
        find_path(_OPENCL_CPP_INCLUDE_DIRS CL/cl.hpp PATHS "${_OPENCL_INC_CAND}" ENV OpenCL_INCPATH)

    else(WIN32)

        # Unix style platforms
        find_library(OpenCL_LIBRARIES OpenCL PATHS ENV LD_LIBRARY_PATH ENV OpenCL_LIBPATH)

        get_filename_component(_OPENCL_LIB_DIR ${OpenCL_LIBRARIES} PATH)
        get_filename_component(_OPENCL_INC_CAND ${_OPENCL_LIB_DIR}/../../include ABSOLUTE)

        # The AMD SDK currently does not place its headers in /usr/include,
        # therefore also search relative to the library
        find_path(OpenCL_INCLUDE_DIRS CL/cl.h PATHS ${_OPENCL_INC_CAND} "/usr/local/cuda/include" "/opt/AMDAPP/include" ENV OpenCL_INCPATH)
        find_path(_OPENCL_CPP_INCLUDE_DIRS CL/cl.hpp PATHS ${_OPENCL_INC_CAND} "/usr/local/cuda/include" "/opt/AMDAPP/include" ENV OpenCL_INCPATH)

    endif(WIN32)

endif(APPLE)

find_package_handle_standard_args(OpenCL DEFAULT_MSG OpenCL_LIBRARIES OpenCL_INCLUDE_DIRS)

if(_OPENCL_CPP_INCLUDE_DIRS)
    set( OpenCL_HAS_CPP_BINDINGS TRUE )
    list( APPEND OpenCL_INCLUDE_DIRS ${_OPENCL_CPP_INCLUDE_DIRS} )
    # This is often the same, so clean up
    list( REMOVE_DUPLICATES OpenCL_INCLUDE_DIRS )
endif(_OPENCL_CPP_INCLUDE_DIRS)

unset(_OPENCL_CPP_INCLUDE_DIRS CACHE)

mark_as_advanced(
  OpenCL_LIBRARIES
  OpenCL_INCLUDE_DIRS
)
