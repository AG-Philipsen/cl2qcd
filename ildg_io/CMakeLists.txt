# Build specifications for the ildg I/O functionality
#
# Copyright (c) 2014 Christopher Pinke <pinke@th.physik.uni-frankfurt.de> 
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

#
# Collect all code in this part of the source tree
#

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

add_unit_test_withLibraries(ildgIo_gaugefield ildgIo_gaugefield_test ildg_io)
add_unit_test_withLibraries(matrixSu3_utilities matrixSu3_utilities_test ildg_io)

