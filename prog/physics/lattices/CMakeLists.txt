# Build specifications for the contents of the hardware package
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

#
# Collect all code in this part of the source tree
#

add_library(lattices
	gaugefield.cpp
	spinorfield.cpp
	spinorfield_eo.cpp
	scalar_complex.cpp
	swappable_spinorfield.cpp
	gaugemomenta.cpp
)
target_link_libraries(lattices
	hardware
	physics
)

#
# Definition of tets
#
add_unit_test(physics::lattices::Gaugefield gaugefield_test)
add_unit_test(physics::lattices::Spinorfield spinorfield_test)
add_unit_test(physics::lattices::Spinorfield_eo spinorfield_eo_test)
add_unit_test(physics::lattices::Scalar<hmc_complex> scalar_complex_test)
add_unit_test(physics::lattices::SwappableSpinorfield swappable_spinorfield_test)
add_unit_test(physics::lattices::Gaugemomenta gaugemomenta_test)
