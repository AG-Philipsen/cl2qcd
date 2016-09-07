# Build specifications for the contents of the hardware package
#
# Copyright (c) 2015 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
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

add_unit_test_withLibraries(hardware::Parameters hardwareParameters_test  meta)
add_unit_test_withLibraries(hardware::code::Parameters openClKernelParameters_test meta)
add_unit_test_withLibraries(physics::Parameters physicsParameters_test meta exceptions)
add_unit_test_withLibraries(physics::observables::Parameters observablesParameters_test meta)
add_unit_test_withLibraries(physics::lattices::Parameters latticesParameters_test meta)
add_unit_test_withLibraries(physics::fermionmatrix::Parameters fermionmatrixParameters_test meta)
add_unit_test_withLibraries(physics::algorithms::Parameters algorithmsParameters_test meta)
add_unit_test_withLibraries(physics::interfacesHandler interfacesHandler_test meta exceptions)