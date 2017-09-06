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

add_library(gaugeObservables
	gaugeObservables.cpp
)	

add_library(observables
	wilsonTwoFlavourChiralCondensate.cpp
	wilsonTwoFlavourCorrelators.cpp
	staggeredChiralCondensate.cpp
	staggeredTwoFlavourCorrelators.cpp
)

target_link_libraries(observables
	lattices
	gaugeObservables
	algorithms
	fermionmatrix
	prng
	meta
	hardware
	host_functionality
	ildg_io
)

#
# Definition of tests
#

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/conf.00200 conf.00200 COPYONLY )

add_unit_test_withLibraries(physics::gaugeObservables gaugeObservables_test observables)
add_unit_test_withLibraries(physics::observables::wilson::twoFlavourChiralCondensate wilsonTwoFlavourChiralCondensate_test lattices observables)
add_test(physics::observables::wilson::twoFlavourChiralCondensate_GPU wilsonTwoFlavourChiralCondensate_test --use_gpu=true)
add_unit_test_withLibraries(physics::observables::wilson::twoFlavourCorrelators wilsonTwoFlavourCorrelators_test lattices observables)
add_unit_test_withLibraries(physics::observables::staggered::ChiralCondensate staggeredChiralCondensate_test lattices observables)
add_unit_test_withLibraries(physics::observables::staggered::twoFlavourCorrelators staggeredTwoFlavourCorrelators_test lattices observables)

