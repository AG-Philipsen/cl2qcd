# Build specifications for the contents of the physics package
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

#
# Collect all code in this part of the source tree
#

add_library(physics
	prng.cpp
	sources.cpp
	../host_random.cpp
)
add_modules(physics
	lattices
	fermionmatrix
	algorithms
)
target_link_libraries(physics
	meta
	hardware
)

#
# Definition of tets
#

add_unit_test(physics::PRNG prng_test)
add_unit_test(physics::<sources> sources_test)
