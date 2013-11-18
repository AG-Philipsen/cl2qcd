# Build specifications for the contents of the buffers package
#
# Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

#
# Collect all code in this part of the source tree
#

add_library(transfer
	transfer.cpp
	ocl_copy.cpp
	async_ocl_copy.cpp
	dgma.cpp
)
target_link_libraries(transfer
	hardware
	buffers
	optimal
)

#
# Definition of tets
#

