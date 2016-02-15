# Build specifications for the contents of the hardware package
#
# Copyright (c) 2016 Francesca Cuteri <cuteri@th.physik.uni-frankfurt.de>
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

add_library(geometry
	latticeExtents.cpp
	index.cpp
	latticeGrid.cpp
	parallelization.cpp
)

target_link_libraries(geometry einhard logger)

#
# Definition of tests
#
add_executable(index index_test.cpp)
target_link_libraries(index 	${Boost_LIBRARIES} geometry )
add_test(Index index )

add_executable(latticeExtents latticeExtents_test.cpp)
target_link_libraries(latticeExtents 	${Boost_LIBRARIES} geometry )
add_test(LatticeExtents latticeExtents )

add_executable(latticeGrid latticeGrid_test.cpp)
target_link_libraries(latticeGrid 	${Boost_LIBRARIES} geometry )
add_test(Grid latticeGrid )

add_executable(parallelizationTest parallelization_test.cpp)
target_link_libraries(parallelizationTest 	${Boost_LIBRARIES} geometry )
add_test(parallelization parallelizationTest )
