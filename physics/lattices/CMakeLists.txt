# Build specifications for the contents of the hardware package
#
# Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

add_library(gaugefield
	gaugefield.cpp
# 	../../hardware/lattices/gaugefield.cpp
# 	../../hardware/lattices/gaugemomenta.cpp
# 	../../hardware/lattices/spinorfield.cpp
# 	../../hardware/lattices/prng.cpp
)

target_link_libraries(gaugefield
	hardware_lattices
)

add_library(lattices
	spinorfield.cpp
	spinorfield_eo.cpp
	scalar_complex.cpp
	algebra_real.cpp
	swappable_spinorfield.cpp
	gaugemomenta.cpp
	staggeredfield_eo.cpp
	rooted_staggeredfield_eo.cpp
)
target_link_libraries(lattices
	hardware
	gaugefield
	ildg_io
	host_functionality
	prng
	observables
)

#
# Definition of tests
#
add_unit_test_withLibraries(physics::lattices::Gaugefield gaugefield_test gaugefield meta hardware gaugeObservables prng ildg_io host_functionality)
add_unit_test_withLibraries(physics::lattices::Spinorfield spinorfield_test lattices)
add_unit_test_withLibraries(physics::lattices::Spinorfield_eo spinorfield_eo_test lattices)
add_unit_test_withLibraries(physics::lattices::Scalar<hmc_complex> scalar_complex_test lattices)
add_unit_test_withLibraries(physics::lattices::SwappableSpinorfield swappable_spinorfield_test lattices)
add_unit_test_withLibraries(physics::lattices::Gaugemomenta gaugemomenta_test lattices)
add_unit_test_withLibraries(physics::lattices::Staggeredfield_eo staggeredfield_eo_test lattices)
add_unit_test_withLibraries(physics::lattices::Algebra_Real algebra_real_test lattices)
add_unit_test_withLibraries(physics::lattices::Vector vector_test lattices)

add_unit_test(physics::lattices::Rooted_Staggeredfield_eo rooted_staggeredfield_eo_test)
