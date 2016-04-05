# Build specifications for the contents of the hardware package
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

add_library(hardware
	system.cpp
	device.cpp
	device_info.cpp
	profiling_data.cpp
	synchronization_event.cpp
	opencl_compiler.cpp
)

add_library(hardwareTestUtilities
	hardwareTestUtilities.cpp
)


add_modules(hardware
	buffers
	code
	transfer
)

#TODO: the following two lines are a workaround for the linking against the hardware_lattices library whose name cannot coincide with the name of the corresponding directory because of homonymy with (physics/)lattices
add_subdirectory(lattices)
set_property(GLOBAL APPEND PROPERTY DOC_SOURCE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lattices")
target_link_libraries(hardware 	${OPENCL_LIBRARIES} exceptions crypto meta geometry hardware_lattices)

#
# Definition of tests
#

add_unit_test_withLibraries(hardware::System system_test hardware crypto)
add_unit_test_withLibraries(hardware::Device device_test hardware crypto)
add_unit_test_withLibraries(hardware::ProfilingData profiling_data_test hardware crypto)
add_unit_test_withLibraries(hardware::SynchronizationEvent synchronization_event_test hardware crypto)

add_executable(openclCompiler_test opencl_compiler_test.cpp opencl_compiler.cpp ../host_functionality/logger.cpp ../executables/exceptions.cpp)
target_link_libraries(openclCompiler_test ${Boost_LIBRARIES} ${OPENCL_LIBRARIES} crypto einhard)
add_test(openclCompiler openclCompiler_test)

