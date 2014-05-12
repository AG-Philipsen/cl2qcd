# Build specifications for the contents of the code package
#
# Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>,
# 	        2014 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
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

add_library(code
		opencl_module.cpp
		gaugefield.cpp
		prng.cpp
		heatbath.cpp
		kappa.cpp
		correlator.cpp
		correlator_staggered.cpp
		real.cpp
		complex.cpp
		spinors.cpp
		spinors_staggered.cpp
		fermions.cpp
		fermions_staggered.cpp
		gaugemomentum.cpp
		molecular_dynamics.cpp
		buffer.cpp
)

target_link_libraries(code
		meta
		buffers
		klepsydra
		einhard
)

#
# Definition of tests
#

add_library(kernelTestUtilities
	testUtilities.cpp
	kernelTester.cpp
)

target_link_libraries(kernelTestUtilities
	${Boost_LIBRARIES}
	host_functionality
	meta
	hardware
	# todo: remove this dependency
	lattices
)

# todo: perhaps one can remove the quite unspecific dependency on physics library
add_unit_test_withLibraries(kernelTester_test KernelTester_test kernelTestUtilities ildg_io physics)

