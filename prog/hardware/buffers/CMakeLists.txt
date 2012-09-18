# Build specifications for the contents of the buffers package
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

#
# Collect all code in this part of the source tree
#

add_library(buffers
	buffer.cpp
)

#
# Definition of tets
#

add_unit_test(hardware::buffers::Buffer buffer_test)
