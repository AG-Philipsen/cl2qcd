# Boost tests
#
# Copyright 2014 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
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

add_executable(boost_test boost_test.cpp)
target_link_libraries(boost_test
	${Boost_LIBRARIES}
)
add_test(boost_ARGC_1 boost_test --run_test=BOOST_ARGUMENTS/BOOST_ARGC_1)
add_test(boost_ARGC_2 boost_test ${CMAKE_CURRENT_SOURCE_DIR} --run_test=BOOST_ARGUMENTS/BOOST_ARGC_2)
add_test(boost_ARGV boost_test "firstArgument" "secondArgument" "thirdArgument" --run_test=BOOST_ARGUMENTS/BOOST_ARGV )
