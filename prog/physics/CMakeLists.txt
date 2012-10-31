# Build specifications for the contents of the physics package
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

#
# Collect all code in this part of the source tree
#

add_library(physics
	prng.cpp
)
add_modules(physics
	lattices
#	algorithms
)

#
# Definition of tets
#

add_unit_test(physics::PRNG prng_test)
