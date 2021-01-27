# Script to locate the Nettle library
#
# Copyright (c) 2013 Matthias Bach
# Copyright (c) 2017-2018,2021 Alessandro Sciarra
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
    find_library(NETTLE_LIBRARIES nettle DOC "Nettle lib for OSX" ENV DYLD_LIBRARY_PATH ENV LIBRARY_PATH)
else(APPLE)
    if(WIN32)
        find_library(NETTLE_LIBRARIES nettle ENV PATH)
    else(WIN32)
        # Unix style platforms
        find_library(NETTLE_LIBRARIES nettle ENV LD_LIBRARY_PATH ENV LIBRARY_PATH)
    endif(WIN32)
endif(APPLE)

find_path(NETTLE_INCLUDE_DIR nettle/nettle-meta.h DOC "Include for Nettle")

# Also search relative to lib (git build)
if(NOT NETTLE_INCLUDE_DIR)
    get_filename_component(_Nettle_LIB_DIR ${NETTLE_LIBRARIES} PATH)
    get_filename_component(_Nettle_INC_CAND ${_Nettle_LIB_DIR}/../include ABSOLUTE)
    find_path(NETTLE_INCLUDE_DIR nettle/nettle-meta.h PATHS ${_Nettle_INC_CAND})
endif(NOT NETTLE_INCLUDE_DIR)

find_package_handle_standard_args(Nettle DEFAULT_MSG NETTLE_LIBRARIES NETTLE_INCLUDE_DIR)

mark_as_advanced(
  NETTLE_LIBRARIES
  NETTLE_INCLUDE_DIR
)
