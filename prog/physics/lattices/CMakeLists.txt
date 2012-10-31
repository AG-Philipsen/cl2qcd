# Build specifications for the contents of the hardware package
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

#
# Collect all code in this part of the source tree
#

add_library(lattices
	gaugefield.cpp
)

#
# Definition of tets
#
add_unit_test(physics::lattices::Gaugefield gaugefield_test)
