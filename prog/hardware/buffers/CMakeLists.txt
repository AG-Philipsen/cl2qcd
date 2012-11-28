# Build specifications for the contents of the buffers package
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

#
# Collect all code in this part of the source tree
#

add_library(buffers
	buffer.cpp
	prng_buffer.cpp
	su3.cpp
	spinor.cpp
	gaugemomentum.cpp
	3x3.cpp
)
target_link_libraries(buffers
	hardware
	code
	optimal
)

#
# Definition of tets
#

add_unit_test(hardware::buffers::Buffer buffer_test)
add_unit_test(hardware::buffers::PRNGBuffer prng_buffer_test)
add_unit_test(hardware::buffers::SU3 su3_test)
add_unit_test(hardware::buffers::Plain plain_test)
add_unit_test(hardware::buffers::Spinor spinor_test)
add_unit_test(hardware::buffers::Gaugemomentum gaugemomentum_test)
add_unit_test(hardware::buffers::Matrix3x3 3x3_test)
