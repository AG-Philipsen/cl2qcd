# Script to locate the MPFR library
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
    find_library(MPFR_LIBRARIES mpfr DOC "MPFR lib for OSX" ENV DYLD_LIBRARY_PATH ENV LIBRARY_PATH)
else(APPLE)
    if(WIN32)
        find_library(MPFR_LIBRARIES mpfr ENV PATH)
    else(WIN32)
        # Unix style platforms
        find_library(MPFR_LIBRARIES mpfr ENV LD_LIBRARY_PATH ENV LIBRARY_PATH)
    endif(WIN32)
endif(APPLE)

find_path(MPFR_INCLUDE_DIR mpfr.h DOC "Include for MPFR")

# Also search relative to lib (git build)
if(NOT MPFR_INCLUDE_DIR)
    get_filename_component(_MPFR_LIB_DIR ${MPFR_LIBRARIES} PATH)
    get_filename_component(_MPFR_INC_CAND ${_MPFR_LIB_DIR}/../include ABSOLUTE)
    find_path(MPFR_INCLUDE_DIR mpfr.h PATHS ${_MPFR_INC_CAND})
endif(NOT MPFR_INCLUDE_DIR)

find_package_handle_standard_args(MPFR DEFAULT_MSG MPFR_LIBRARIES MPFR_INCLUDE_DIR)

mark_as_advanced(
  MPFR_LIBRARIES
  MPFR_INCLUDE_DIR
)
