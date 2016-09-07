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
		buffers
		klepsydra
		einhard
		hardwareTestUtilities
)

#
# Definition of tests
#

add_subdirectory(miscellaneousTests)

# default random number seed
configure_file( ${CMAKE_SOURCE_DIR}/rand_seeds rand_seeds COPYONLY )

add_library(kernelTestUtilities
	testUtilities.cpp
	kernelTester.cpp
	GaugefieldTester.cpp
	SpinorTester.cpp
	GaugemomentumTester.cpp
	SpinorStaggeredTester.cpp
)

target_link_libraries(kernelTestUtilities
	${Boost_LIBRARIES}
	host_functionality
	hardware
	geometry
)

add_unit_test_withLibraries(kernelTester_test KernelTester_test kernelTestUtilities)

add_executable(hardware_code_gaugefield_test gaugefield_test.cpp)
target_link_libraries(hardware_code_gaugefield_test kernelTestUtilities)
add_test(hardware::code::Gaugefield_CPU       hardware_code_gaugefield_test --use_gpu=false --use_rec12=false)
add_test(hardware::code::Gaugefield_CPU_REC12 hardware_code_gaugefield_test --use_gpu=false --use_rec12=true)
add_test(hardware::code::Gaugefield_GPU       hardware_code_gaugefield_test --use_gpu=true --use_rec12=false)
add_test(hardware::code::Gaugefield_GPU_REC12 hardware_code_gaugefield_test --use_gpu=true --use_rec12=true)

add_executable(hardware_code_spinors_test spinors_test.cpp)
target_link_libraries(hardware_code_spinors_test kernelTestUtilities)
add_test(hardware_code_spinors_CPU       hardware_code_spinors_test --use_gpu=false)
add_test(hardware_code_spinors_GPU       hardware_code_spinors_test --use_gpu=true)

add_executable(hardware_code_spinors_merged_kernels_test spinors_merged_kernels_test.cpp)
target_link_libraries(hardware_code_spinors_merged_kernels_test kernelTestUtilities )
add_test(hardware_code_spinors_merged_kernels_test_CPU       hardware_code_spinors_merged_kernels_test --use_gpu=false --use_rec12=false)
add_test(hardware_code_spinors_merged_kernels_test_GPU       hardware_code_spinors_merged_kernels_test --use_gpu=true --use_rec12=false)

add_executable(hardware_code_fermions_test fermions_test.cpp)
target_link_libraries(hardware_code_fermions_test kernelTestUtilities )
add_test(hardware_code_fermions_CPU       hardware_code_fermions_test --use_gpu=false --use_rec12=false)
add_test(hardware_code_fermions_CPU_REC12 hardware_code_fermions_test --use_gpu=false --use_rec12=true)
add_test(hardware_code_fermions_GPU       hardware_code_fermions_test --use_gpu=true --use_rec12=false)
add_test(hardware_code_fermions_GPU_REC12 hardware_code_fermions_test --use_gpu=true --use_rec12=true)

add_executable(hardware_code_fermions_merged_kernels_test fermions_merged_kernels_test.cpp)
target_link_libraries(hardware_code_fermions_merged_kernels_test kernelTestUtilities)
add_test(hardware_code_fermions_merged_kernels_test_CPU       hardware_code_fermions_merged_kernels_test --use_gpu=false --use_rec12=false)
add_test(hardware_code_fermions_merged_kernels_test_CPU_REC12 hardware_code_fermions_merged_kernels_test --use_gpu=false --use_rec12=true)
add_test(hardware_code_fermions_merged_kernels_test_GPU       hardware_code_fermions_merged_kernels_test --use_gpu=true --use_rec12=false)
add_test(hardware_code_fermions_merged_kernels_test_GPU_REC12 hardware_code_fermions_merged_kernels_test --use_gpu=true --use_rec12=true)

add_executable(hardware_code_spinors_staggered_test spinors_staggered_test.cpp)
target_link_libraries(hardware_code_spinors_staggered_test kernelTestUtilities)
add_test(hardware_code_spinors_staggered_CPU       hardware_code_spinors_staggered_test --use_gpu=false --use_rec12=false)
add_test(hardware_code_spinors_staggered_GPU       hardware_code_spinors_staggered_test --use_gpu=true --use_rec12=false)

add_executable(hardware_code_correlator_test correlator_test.cpp)
target_link_libraries(hardware_code_correlator_test kernelTestUtilities )
add_test(hardware_code_correlator_CPU       hardware_code_correlator_test --use_gpu=false --use_rec12=false)
add_test(hardware_code_correlator_GPU       hardware_code_correlator_test --use_gpu=true --use_rec12=false)

add_executable(hardware_code_real_test real_test.cpp)
target_link_libraries(hardware_code_real_test kernelTestUtilities)
add_test(hardware_code_real_CPU       hardware_code_real_test --use_gpu=false --use_rec12=false)
add_test(hardware_code_real_GPU       hardware_code_real_test --use_gpu=true --use_rec12=false)

# todo: these contain only dummy tests
add_executable(hardware_code_gaugemomentum_test gaugemomentum_test.cpp)
target_link_libraries(hardware_code_gaugemomentum_test kernelTestUtilities )
add_test(hardware_code_gaugemomentum_CPU       hardware_code_gaugemomentum_test --use_gpu=false --use_rec12=false)
add_test(hardware_code_gaugemomentum_CPU_REC12 hardware_code_gaugemomentum_test --use_gpu=false --use_rec12=true)
add_test(hardware_code_gaugemomentum_GPU       hardware_code_gaugemomentum_test --use_gpu=true --use_rec12=false)
add_test(hardware_code_gaugemomentum_GPU_REC12 hardware_code_gaugemomentum_test --use_gpu=true --use_rec12=true)

add_executable(hardware_code_molecular_dynamics_test molecular_dynamics_test.cpp)
target_link_libraries(hardware_code_molecular_dynamics_test kernelTestUtilities )
add_test(hardware_code_molecular_dynamics_CPU       hardware_code_molecular_dynamics_test --use_gpu=false --use_rec12=false)
add_test(hardware_code_molecular_dynamics_CPU_REC12 hardware_code_molecular_dynamics_test --use_gpu=false --use_rec12=true)
add_test(hardware_code_molecular_dynamics_GPU       hardware_code_molecular_dynamics_test --use_gpu=true --use_rec12=false)
add_test(hardware_code_molecular_dynamics_GPU_REC12 hardware_code_molecular_dynamics_test --use_gpu=true --use_rec12=true)

add_executable(hardware_code_fermions_staggered_test fermions_staggered_test.cpp)
target_link_libraries(hardware_code_fermions_staggered_test kernelTestUtilities)
add_test(hardware_code_fermions_staggered_CPU       hardware_code_fermions_staggered_test --use_gpu=false --use_rec12=false)
add_test(hardware_code_fermions_staggered_CPU_REC12 hardware_code_fermions_staggered_test --use_gpu=false --use_rec12=true)
add_test(hardware_code_fermions_staggered_GPU       hardware_code_fermions_staggered_test --use_gpu=true --use_rec12=false)
add_test(hardware_code_fermions_staggered_GPU_REC12 hardware_code_fermions_staggered_test --use_gpu=true --use_rec12=true)

add_executable(hardware_code_correlator_staggered_test correlator_staggered_test.cpp)
target_link_libraries(hardware_code_correlator_staggered_test kernelTestUtilities)
add_test(hardware_code_correlator_staggered_CPU       hardware_code_correlator_staggered_test --use_gpu=false --use_rec12=false)
add_test(hardware_code_correlator_staggered_GPU       hardware_code_correlator_staggered_test --use_gpu=true --use_rec12=false)

add_executable(hardware_code_complex_test complex_test.cpp)
target_link_libraries(hardware_code_complex_test kernelTestUtilities )
add_test(hardware_code_complex_CPU       hardware_code_complex_test --use_gpu=false --use_rec12=false)
add_test(hardware_code_complex_GPU       hardware_code_complex_test --use_gpu=true --use_rec12=false)







