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
)

#
# Definition of tets
#

add_unit_test(hardware::buffers::Buffer buffer_test)
add_unit_test(hardware::buffers::PRNGBuffer prng_buffer_test)
add_unit_test(hardware::buffers::SU3 su3_test)
add_unit_test(hardware::buffers::ScalarBuffer scalar_buffer_test)
