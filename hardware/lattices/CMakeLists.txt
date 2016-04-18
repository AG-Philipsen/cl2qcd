# Build specifications for the contents of the buffers package
#
# Copyright (c) 2016 Francesca Cuteri
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

add_library(hardware_lattices
	prng.cpp
	gaugefield.cpp
	gaugemomenta.cpp
	spinorfield.cpp
	spinorfield_eo.cpp
	staggeredfield_eo.cpp
)

target_link_libraries(hardware_lattices
	hardware
	code
	buffers
	host_functionality
)

#
# Definition of tests
#
# add_executable(Gaugefield gaugefield_test.cpp)
# target_link_libraries(Gaugefield 	${Boost_LIBRARIES} geometry crypto)
# add_test(gaugefield Gaugefield )

# add_unit_test_withLibraries(hardware::lattices::Gaugefield hardware_lattices_gaugefield_test lattices crypto)
# add_unit_test_withLibraries(hardware::lattices::Gaugemomenta hardware_lattices_gaugemomenta_test lattices crypto)
# add_unit_test_withLibraries(hardware::lattices::Spinorfield hardware_lattices_spinorfield_test lattices crypto)
# add_unit_test_withLibraries(hardware::lattices::Spinorfield_eo hardware_lattices_spinorfield_eo_test lattices crypto)
# add_unit_test_withLibraries(hardware::lattices::Prng hardware_lattices_prng_test lattices crypto)
