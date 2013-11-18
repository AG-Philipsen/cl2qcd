# Build specifications for the contents of the fermionmatrix package
#
# Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>

#
# Collect all code in this part of the source tree
#

add_library(fermionmatrix
	fermionmatrix.cpp
	explicit.cpp
	fermionmatrix_stagg.cpp
	explicit_stagg.cpp
)
target_link_libraries(fermionmatrix
	meta
	hardware
	lattices
)

add_unit_test(physics::fermionmatrix::<explicit> explicit_test)
add_unit_test(physics::fermionmatrix::<fermionmatrix> fermionmatrix_test)
add_unit_test(physics::fermionmatrix::<explicit_stagg> explicit_stagg_test)
add_unit_test(physics::fermionmatrix::<fermionmatrix_stagg> fermionmatrix_stagg_test)
