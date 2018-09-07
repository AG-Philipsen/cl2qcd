# Build specifications for the contents of the hardware package
#
# Copyright (c) 2016 Christopher Pinke
# Copyright (c) 2016 Francesca Cuteri
# Copyright (c) 2018 Alessandro Sciarra
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.

add_library(geometry
    latticeExtents.cpp
    index.cpp
    latticeGrid.cpp
    parallelization.cpp
)

target_link_libraries(geometry
    einhard
    logger
)

#
# Definition of tests
#
add_unit_test(NAME geometry/index           LIBRARIES geometry)
add_unit_test(NAME geometry/latticeExtents  LIBRARIES geometry)
add_unit_test(NAME geometry/latticeGrid     LIBRARIES geometry)
add_unit_test(NAME geometry/parallelization LIBRARIES geometry)
