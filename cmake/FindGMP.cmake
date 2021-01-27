# Script to locate the GMP library
#
# Copyright (c) 2013 Matthias Bach
# Copyright (c) 2018,2021 Alessandro Sciarra
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

find_package(PackageHandleStandardArgs)

if(APPLE)
    find_library(GMP_LIBRARIES gmp DOC "GMP lib for OSX" ENV DYLD_LIBRARY_PATH ENV LIBRARY_PATH)
else(APPLE)
    if(WIN32)
        find_library(GMP_LIBRARIES gmp ENV PATH)
    else(WIN32)
        # Unix style platforms
        find_library(GMP_LIBRARIES gmp ENV LD_LIBRARY_PATH ENV LIBRARY_PATH)
    endif(WIN32)
endif(APPLE)

find_path(GMP_INCLUDE_DIR gmp.h DOC "Include for GMP")

# Also search relative to lib (git build)
if(NOT GMP_INCLUDE_DIR)
    get_filename_component(_GMP_LIB_DIR ${GMP_LIBRARIES} PATH)
    get_filename_component(_GMP_INC_CAND ${_GMP_LIB_DIR}/../include ABSOLUTE)
    find_path(GMP_INCLUDE_DIR gmp.h PATHS ${_GMP_INC_CAND})
endif(NOT GMP_INCLUDE_DIR)

find_package_handle_standard_args(GMP DEFAULT_MSG GMP_LIBRARIES GMP_INCLUDE_DIR)

mark_as_advanced(
  GMP_LIBRARIES
  GMP_INCLUDE_DIR
)
