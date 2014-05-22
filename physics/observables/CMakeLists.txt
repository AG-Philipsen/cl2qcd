# Build specifications for physics::observables
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

add_library(observables
	gaugeObservables.cpp
	wilsonTwoFlavourChiralCondensate.cpp
)

target_link_libraries(observables
	algorithms
	prng
	meta
	hardware
	host_functionality
	ildg_io
	lattices
)

#
# Definition of tests
#

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/conf.00200 conf.00200 COPYONLY )

add_unit_test_withLibraries(physics::gaugeObservables gaugeObservables_test observables)
add_unit_test_withLibraries(physics::observables::wilson::twoFlavourChiralCondensate wilsonTwoFlavourChiralCondensate_test observables)

