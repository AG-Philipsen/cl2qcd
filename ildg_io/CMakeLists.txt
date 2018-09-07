# Build specifications for the ildg I/O functionality
#
# Copyright (c) 2014,2015 Christopher Pinke
# Copyright (c) 2016 Francesca Cuteri
# Copyright (c) 2018 Alessandro Sciarra
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

add_library(ildg_io
    ildgIo_gaugefield.cpp
    ildgIo.cpp
    ildgIoParameters.cpp
    matrixSu3_utilities.cpp
)

target_link_libraries(ildg_io
    ${LIME_LIBRARIES}
    ${LIBXML2_LIBRARIES}
    # this needs to be a stand-alone library...
    einhard
    logger
    exceptions
    meta
    host_functionality
    lime
    sourcefileParameters
    geometry
)

add_subdirectory(lime)
add_subdirectory(sourcefileParameters)

add_unit_test(NAME ildg_io/ildgIo_gaugefield   LIBRARIES ildg_io)
add_unit_test(NAME ildg_io/matrixSu3_utilities LIBRARIES ildg_io)
