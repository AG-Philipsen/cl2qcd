# Copyright 2012, 2013 Matthias Bach
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

add_executable(ranluxcltest ranluxcltest.cpp)
target_link_libraries( ranluxcltest
	${OPENCL_LIBRARIES}
)
add_test(NAME Ranluxcl_CPU COMMAND ranluxcltest 1 0 1 1)
add_test(NAME Ranluxcl_GPU COMMAND ranluxcltest 1 1 1 1)
set_tests_properties(Ranluxcl_CPU Ranluxcl_GPU PROPERTIES PASS_REGULAR_EXPRESSION "Correctness check PASSED!")
