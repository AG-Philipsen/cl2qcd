# Build specifications for the contents of the fermionmatrix package
#
# Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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
