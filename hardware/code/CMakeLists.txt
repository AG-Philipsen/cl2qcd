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

# todo: remove dependencies on optimal library!

# default random number seed
configure_file( ${CMAKE_SOURCE_DIR}/rand_seeds rand_seeds COPYONLY )

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/conf.00200 conf.00200 COPYONLY )

add_library(kernelTestUtilities
	testUtilities.cpp
	kernelTester.cpp
	SpinorTester.cpp
	GaugemomentumTester.cpp
	SpinorStaggeredTester.cpp
)

target_link_libraries(kernelTestUtilities
	${Boost_LIBRARIES}
	host_functionality
	meta
	hardware
)

add_unit_test_withLibraries(kernelTester_test KernelTester_test kernelTestUtilities)

add_executable(hardware_code_gaugefield_test gaugefield_test.cpp)
target_link_libraries(hardware_code_gaugefield_test
	kernelTestUtilities
	lattices
)
add_test(hardware::code::Gaugefield_CPU       hardware_code_gaugefield_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=false)
add_test(hardware::code::Gaugefield_CPU_REC12 hardware_code_gaugefield_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=true)
add_test(hardware::code::Gaugefield_GPU       hardware_code_gaugefield_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=false)
add_test(hardware::code::Gaugefield_GPU_REC12 hardware_code_gaugefield_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=true)

add_executable(hardware_code_spinors_test spinors_test.cpp)
target_link_libraries(hardware_code_spinors_test
	prng
	lattices
	kernelTestUtilities
)
add_test(hardware_code_spinors_CPU       hardware_code_spinors_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false)
add_test(hardware_code_spinors_GPU       hardware_code_spinors_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true)

add_executable(hardware_code_fermions_test fermions_test.cpp)
target_link_libraries(hardware_code_fermions_test
	prng
	lattices
	kernelTestUtilities
)
add_test(hardware_code_fermions_CPU       hardware_code_fermions_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=false)
add_test(hardware_code_fermions_CPU_REC12 hardware_code_fermions_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=true)
add_test(hardware_code_fermions_GPU       hardware_code_fermions_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=false)
add_test(hardware_code_fermions_GPU_REC12 hardware_code_fermions_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=true)

add_executable(hardware_code_gaugemomentum_test gaugemomentum_test.cpp)
target_link_libraries(hardware_code_gaugemomentum_test
	prng
	kernelTestUtilities
)
add_test(hardware_code_gaugemomentum_CPU       hardware_code_gaugemomentum_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=false)
add_test(hardware_code_gaugemomentum_CPU_REC12 hardware_code_gaugemomentum_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=true)
add_test(hardware_code_gaugemomentum_GPU       hardware_code_gaugemomentum_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=false)
add_test(hardware_code_gaugemomentum_GPU_REC12 hardware_code_gaugemomentum_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=true)

add_executable(hardware_code_molecular_dynamics_test molecular_dynamics_test.cpp)
target_link_libraries(hardware_code_molecular_dynamics_test
	optimal
	kernelTestUtilities
)
add_test(hardware_code_molecular_dynamics_CPU       hardware_code_molecular_dynamics_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=false)
add_test(hardware_code_molecular_dynamics_CPU_REC12 hardware_code_molecular_dynamics_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=true)
add_test(hardware_code_molecular_dynamics_GPU       hardware_code_molecular_dynamics_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=false)
add_test(hardware_code_molecular_dynamics_GPU_REC12 hardware_code_molecular_dynamics_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=true)

add_executable(hardware_code_spinors_staggered_test spinors_staggered_test.cpp)
target_link_libraries(hardware_code_spinors_staggered_test
	optimal
	kernelTestUtilities
)
add_test(hardware_code_spinors_staggered_CPU       hardware_code_spinors_staggered_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=false)
add_test(hardware_code_spinors_staggered_GPU       hardware_code_spinors_staggered_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=false)

add_executable(hardware_code_fermions_staggered_test fermions_staggered_test.cpp)
target_link_libraries(hardware_code_fermions_staggered_test
	optimal
	kernelTestUtilities
)
add_test(hardware_code_fermions_staggered_CPU       hardware_code_fermions_staggered_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=false)
add_test(hardware_code_fermions_staggered_CPU_REC12 hardware_code_fermions_staggered_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=true)
add_test(hardware_code_fermions_staggered_GPU       hardware_code_fermions_staggered_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=false)
add_test(hardware_code_fermions_staggered_GPU_REC12 hardware_code_fermions_staggered_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=true)

add_executable(hardware_code_complex_test complex_test.cpp)
target_link_libraries(hardware_code_complex_test
	optimal
	kernelTestUtilities
)
add_test(hardware_code_complex_CPU       hardware_code_complex_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=false)
add_test(hardware_code_complex_GPU       hardware_code_complex_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=false)

add_executable(hardware_code_real_test real_test.cpp)
target_link_libraries(hardware_code_real_test
	optimal
	kernelTestUtilities
)
add_test(hardware_code_real_CPU       hardware_code_real_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=false)
add_test(hardware_code_real_GPU       hardware_code_real_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=false)