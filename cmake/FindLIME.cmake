# Script to locate the LIME library
#
# Copyright (c) 2011 Matthias Bach
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
    find_library(LIME_LIBRARIES lime DOC "LIME lib for OSX"  ENV DYLD_LIBRARY_PATH ENV LIBRARY_PATH)
else(APPLE)
    if(WIN32)
        find_library(LIME_LIBRARIES lime ENV PATH)
    else(WIN32)
        # Unix style platforms
        find_library(LIME_LIBRARIES lime ENV LD_LIBRARY_PATH ENV LIBRARY_PATH)
    endif(WIN32)
endif(APPLE)

find_path(LIME_INCLUDE_DIR lime.h DOC "Include for LIME")

# Also search relative to lib (git build)
if(NOT LIME_INCLUDE_DIR)
    get_filename_component(_LIME_LIB_DIR ${LIME_LIBRARIES} PATH)
    get_filename_component(_LIME_INC_CAND ${_LIME_LIB_DIR}/../include ABSOLUTE)
    find_path(LIME_INCLUDE_DIR lime.h PATHS ${_LIME_INC_CAND})
endif(NOT LIME_INCLUDE_DIR)

find_package_handle_standard_args( LIME DEFAULT_MSG LIME_LIBRARIES LIME_INCLUDE_DIR )

mark_as_advanced(
  LIME_LIBRARIES
  LIME_INCLUDE_DIR
)
