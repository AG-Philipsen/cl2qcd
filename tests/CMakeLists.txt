# Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
# Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

message(STATUS "Building Tests")

find_package(PythonInterp REQUIRED)

# default random number seed
configure_file( ${CMAKE_SOURCE_DIR}/rand_seeds rand_seeds COPYONLY )

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/conf.00200 conf.00200 COPYONLY )

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

add_executable(boost_test boost_test.cpp)
target_link_libraries(boost_test
	${Boost_LIBRARIES}
)
add_test(boost_ARGC_1 boost_test --run_test=BOOST_ARGUMENTS/BOOST_ARGC_1)
add_test(boost_ARGC_2 boost_test ${CMAKE_CURRENT_SOURCE_DIR} --run_test=BOOST_ARGUMENTS/BOOST_ARGC_2)
add_test(boost_ARGV boost_test "firstArgument" "secondArgument" "thirdArgument" --run_test=BOOST_ARGUMENTS/BOOST_ARGV )

add_executable(su2su3extend su2su3extend.cpp)
target_link_libraries(su2su3extend
	optimal
	${Boost_LIBRARIES}
)
add_test(Su2su3extend su2su3extend)

add_executable(geometry geometry.cpp)
target_link_libraries(geometry
	optimal
	${Boost_LIBRARIES}
)
add_test(geometry geometry)

add_executable(gaugefield_convert gaugefield_convert.cpp)
target_link_libraries(gaugefield_convert
	optimal
	${Boost_LIBRARIES}
)
add_test(Gaugefield_Conversion gaugefield_convert)

add_executable(calc_staple Calc_staple_test.cpp)
target_link_libraries(calc_staple
	optimal
	${Boost_LIBRARIES}
)

add_test(Calc_staple_1_CPU       calc_staple ${CMAKE_CURRENT_SOURCE_DIR}/staple_input_1 --use_gpu=false --use_rec12=false)
add_test(Calc_staple_1_CPU_REC12 calc_staple ${CMAKE_CURRENT_SOURCE_DIR}/staple_input_1 --use_gpu=false --use_rec12=true)
add_test(Calc_staple_1_GPU       calc_staple ${CMAKE_CURRENT_SOURCE_DIR}/staple_input_1 --use_gpu=true  --use_rec12=false)
add_test(Calc_staple_1_GPU_REC12 calc_staple ${CMAKE_CURRENT_SOURCE_DIR}/staple_input_1 --use_gpu=true  --use_rec12=true)

add_executable(localQ local_Q_test.cpp)
target_link_libraries(localQ
	optimal
	${Boost_LIBRARIES}
)
add_test(LocalQ_1_CPU       localQ ${CMAKE_CURRENT_SOURCE_DIR}/local_Q_input_1 --use_gpu=false --use_rec12=false)
add_test(LocalQ_1_CPU_REC12 localQ ${CMAKE_CURRENT_SOURCE_DIR}/local_Q_input_1 --use_gpu=false --use_rec12=true)
add_test(LocalQ_1_GPU       localQ ${CMAKE_CURRENT_SOURCE_DIR}/local_Q_input_1 --use_gpu=true  --use_rec12=false)
add_test(LocalQ_1_GPU_REC12 localQ ${CMAKE_CURRENT_SOURCE_DIR}/local_Q_input_1 --use_gpu=true  --use_rec12=true)

# ignore some tests
configure_file(CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake @ONLY)

# Test the inverter for correct inversion
## even testnumbers correspond to tests using even-odd-precond., 
## odd numbers to those who don`t
## tests on a cold config
add_test(Inverter1_CPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_rec12=false --use_eo=false)
add_test(Inverter2_CPU_TM_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_rec12=false --use_eo=true)
add_test(Inverter1_CPU_TM_CG ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_rec12=false --use_eo=false --solver=cg)
add_test(Inverter2_CPU_TM_CG_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_rec12=false --use_eo=true --solver=cg)
add_test(Inverter1_CPU_TM_SAVE ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_rec12=false --use_eo=false --solver=bicgstab_save)
add_test(Inverter2_CPU_TM_SAVE_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_rec12=false --use_eo=true --solver=bicgstab_save)
add_test(Inverter3_CPU_TM_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_rec12=true --use_eo=false)
add_test(Inverter4_CPU_TM_EO_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_rec12=true --use_eo=true)
add_test(Inverter5_GPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_rec12=false --use_eo=false --use_cpu=false)
add_test(Inverter6_GPU_TM_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_rec12=false --use_eo=true --use_cpu=false)
add_test(Inverter5_GPU_TM_CG ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_rec12=false --use_eo=false --solver=cg --use_cpu=false)
add_test(Inverter6_GPU_TM_CG_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_rec12=false --use_eo=true --solver=cg --use_cpu=false)
add_test(Inverter5_GPU_TM_SAVE ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_rec12=false --use_eo=false --solver=bicgstab_save --use_cpu=false)
add_test(Inverter6_GPU_TM_SAVE_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_rec12=false --use_eo=true --solver=bicgstab_save --use_cpu=false)
add_test(Inverter7_GPU_TM_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_rec12=true --use_eo=false --use_cpu=false)
add_test(Inverter8_GPU_TM_EO_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_rec12=true --use_eo=true --use_cpu=false)
## tests on a given config
## the config does not contain a valid checksum, the resulting error will be ignored
add_test(Inverter9_CPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0 --use_rec12=false --use_eo=false --ignore_checksum_errors=true)
add_test(Inverter10_CPU_TM_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0 --use_rec12=false --use_eo=true --ignore_checksum_errors=true)
add_test(Inverter9_CPU_TM_CG ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0 --use_rec12=false --use_eo=false --ignore_checksum_errors=true --solver=cg)
add_test(Inverter10_CPU_TM_CG_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0 --use_rec12=false --use_eo=true --ignore_checksum_errors=true --solver=cg)
add_test(Inverter9_CPU_TM_SAVE ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0 --use_rec12=false --use_eo=false --ignore_checksum_errors=true --solver=bicgstab_save)
add_test(Inverter10_CPU_TM_SAVE_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0 --use_rec12=false --use_eo=true --ignore_checksum_errors=true --solver=bicgstab_save)
add_test(Inverter11_CPU_TM_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 
--use_gpu=0 --use_rec12=true --use_eo=false  --ignore_checksum_errors=true)
add_test(Inverter12_CPU_TM_EO_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 
--use_gpu=0 --use_rec12=true --use_eo=true  --ignore_checksum_errors=true)
add_test(Inverter13_GPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=1 --use_cpu=false 
--use_rec12=false --use_eo=false  --ignore_checksum_errors=true)
add_test(Inverter14_GPU_TM_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_cpu=false 
--use_gpu=1 --use_rec12=false --use_eo=true  --ignore_checksum_errors=true)
add_test(Inverter13_GPU_TM_CG ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=1 --use_rec12=false --use_eo=false  --ignore_checksum_errors=true --solver=cg --use_cpu=false)
add_test(Inverter14_GPU_TM_CG_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=1 --use_rec12=false --use_eo=true  --ignore_checksum_errors=true --solver=cg --use_cpu=false)
add_test(Inverter13_GPU_TM_SAVE ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=1 --use_rec12=false --use_eo=false  --ignore_checksum_errors=true --solver=bicgstab_save --use_cpu=false)
add_test(Inverter14_GPU_TM_SAVE_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=1 --use_rec12=false --use_eo=true  --ignore_checksum_errors=true --solver=bicgstab_save --use_cpu=false)
add_test(Inverter15_GPU_TM_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_cpu=false 
--use_gpu=1 --use_rec12=true --use_eo=false  --ignore_checksum_errors=true)  
add_test(Inverter16_GPU_TM_EO_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_cpu=false 
--use_gpu=1 --use_rec12=true --use_eo=true  --ignore_checksum_errors=true)

# Test the correlators
## the inverter test above have tested also the z-dir correlators, the t-dir correlators are pending
add_test(Inverter17_CPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref 
${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_eo=false --corr_dir=0)
add_test(Inverter18_CPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_eo=false --corr_dir=0)
add_test(Inverter19_CPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref3 
${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0 --use_eo=false --corr_dir=0 --ignore_checksum_errors=true )
add_test(Inverter20_CPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref3 
${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=1 --use_eo=false --corr_dir=0 --ignore_checksum_errors=true )

# Tests for pure Wilson fermions
## compared to the inverter tests above, only the fermionmatrix is different
#in that it needs less matrix kernels.
#thus only redo tests 1,2,5,6,9,10,13,14 from above
add_test(Inverter21_CPU_W ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref4 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0  
--use_eo=false --fermact=WILSON)
add_test(Inverter22_CPU_W_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref4 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 
--use_eo=true --fermact=WILSON)
add_test(Inverter23_GPU_W ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref4 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1  --use_cpu=false 
--use_eo=false --fermact=WILSON)
add_test(Inverter24_GPU_W_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref4 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_cpu=false 
--use_eo=true --fermact=WILSON)
add_test(Inverter21_CPU_W_CG ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref4 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_eo=false --fermact=WILSON --solver=cg)
add_test(Inverter22_CPU_W_CG_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref4 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_eo=true --fermact=WILSON --solver=cg)
add_test(Inverter23_GPU_W_CG ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref4 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_eo=false --fermact=WILSON --solver=cg --use_cpu=false)
add_test(Inverter24_GPU_W_CG_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref4 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_eo=true --fermact=WILSON --solver=cg --use_cpu=false)
add_test(Inverter21_CPU_W_SAVE ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref4 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_eo=false --fermact=WILSON --solver=bicgstab_save)
add_test(Inverter22_CPU_W_SAVE_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref4 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_eo=true --fermact=WILSON --solver=bicgstab_save)
add_test(Inverter23_GPU_W_SAVE ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref4 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_eo=false --fermact=WILSON --solver=bicgstab_save --use_cpu=false)
add_test(Inverter24_GPU_W_SAVE_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref4 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_eo=true --fermact=WILSON --solver=bicgstab_save --use_cpu=false)
## tests on a given config
## the config does not contain a valid checksum, the resulting error will be ignored
add_test(Inverter25_CPU_W ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref5 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0  
--use_eo=false --ignore_checksum_errors=true --fermact=WILSON)
add_test(Inverter26_CPU_W_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref5 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0
  --use_eo=true --ignore_checksum_errors=true --fermact=WILSON)
add_test(Inverter27_GPU_W ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref5 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=1 --use_cpu=false 
 --use_eo=false  --ignore_checksum_errors=true --fermact=WILSON)
add_test(Inverter28_GPU_W_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref5 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_cpu=false 
--use_gpu=1 --use_eo=true  --ignore_checksum_errors=true --fermact=WILSON)
add_test(Inverter25_CPU_W_CG ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref5 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0 --use_eo=false --ignore_checksum_errors=true --fermact=WILSON --solver=cg)
add_test(Inverter26_CPU_W_CG_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref5 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0 --use_eo=true --ignore_checksum_errors=true --fermact=WILSON --solver=cg)
add_test(Inverter27_GPU_W_CG ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref5 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=1 --use_eo=false  --ignore_checksum_errors=true --fermact=WILSON --solver=cg --use_cpu=false)
add_test(Inverter28_GPU_W_CG_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref5 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=1 --use_eo=true  --ignore_checksum_errors=true --fermact=WILSON --solver=cg --use_cpu=false)
add_test(Inverter25_CPU_W_SAVE ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref5 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0 --use_eo=false --ignore_checksum_errors=true --fermact=WILSON --solver=bicgstab_save)
add_test(Inverter26_CPU_W_SAVE_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref5 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0 --use_eo=true --ignore_checksum_errors=true --fermact=WILSON --solver=bicgstab_save)
add_test(Inverter27_GPU_W_SAVE ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref5 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=1 --use_eo=false  --ignore_checksum_errors=true --fermact=WILSON --solver=bicgstab_save --use_cpu=false)
add_test(Inverter28_GPU_W_SAVE_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref5 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=1 --use_eo=true  --ignore_checksum_errors=true --fermact=WILSON --solver=bicgstab_save --use_cpu=false)

# Test the chiral condensate
## Suppose that the inverter works fine (tests above)
add_test(Inverter29_CPU_PBP_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref6 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_cpu=1 --measure_pbp=true --measure_corr=false --sourcetype=volume --sourcecontent=one )
add_test(Inverter30_GPU_PBP_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref6 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_cpu=false --measure_pbp=true --measure_corr=false --sourcetype=volume --sourcecontent=one)
add_test(Inverter31_CPU_PBP_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref7_cpu ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_cpu=1 --measure_pbp=true --measure_corr=false --sourcetype=volume --sourcecontent=z4  )
add_test(Inverter32_GPU_PBP_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref7_gpu ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_cpu=false --measure_pbp=true --measure_corr=false --sourcetype=volume --sourcecontent=z4)
add_test(Inverter33_CPU_PBP_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref8_cpu ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_cpu=1 --measure_pbp=true --measure_corr=false --sourcetype=volume --sourcecontent=gaussian)
add_test(Inverter34_GPU_PBP_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref8_gpu ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_cpu=false --measure_pbp=true --measure_corr=false --sourcetype=volume --sourcecontent=gaussian)
add_test(Inverter35_CPU_PBP_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref9_cpu ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_cpu=1 --measure_pbp=true --measure_corr=false --sourcetype=timeslice --sourcecontent=gaussian)
add_test(Inverter36_GPU_PBP_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref9_gpu ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_cpu=false --measure_pbp=true --measure_corr=false --sourcetype=timeslice --sourcecontent=gaussian)
add_test(Inverter37_CPU_PBP_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref10_cpu ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_cpu=1 --measure_pbp=true --measure_corr=false --sourcetype=zslice --sourcecontent=gaussian)
add_test(Inverter38_GPU_PBP_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref10_gpu ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_cpu=false --measure_pbp=true --measure_corr=false --sourcetype=zslice --sourcecontent=gaussian)

# Test the heatbath
add_test(Heatbath1_CPU_DP ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHeatbath.py ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_input_1)
add_test(Heatbath3_GPU_DP ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHeatbath.py ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_input_3 --use_cpu=false)
add_test(Heatbath5_CPU_SP ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHeatbath.py ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_input_5)
add_test(Heatbath7_GPU_SP ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHeatbath.py ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_input_7 --use_cpu=false)

# Test the HMC
## odd testnumbers correspond to tests using even-odd-precond., 
## even numbers to those who don`t
## LF: Leapfrog
## 2MN: 2mn
## XTS: X timescales
## CONF: starting from config
## ANTI: antiperiodic boundary conditions
add_test(HMC1_CPU_TM_EO_LF_1TS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_1)
add_test(HMC2_CPU_TM_NONEO_LF_1TS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_2)
add_test(HMC3_CPU_TM_EO_2MN_1TS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_3)
add_test(HMC5_CPU_TM_EO_LF_CG_1TS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_5)
add_test(HMC7_CPU_TM_EO_2MN_2TS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_7)
add_test(HMC9_CPU_TM_EO_2MN_3TS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_9)
add_test(HMC11_CPU_TM_EO_LF_1TS_CONF ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref2  ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_11)
add_test(HMC12_CPU_TM_NONEO_LF_1TS_CONF ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_12)
add_test(HMC13_CPU_PGT_NONEO_LF_1TS_COLD ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref3 ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_13)
add_test(HMC14_CPU_PGT_RECT_NONEO_LF_1TS_COLD ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref4 ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_14)
add_test(HMC15_CPU_TM_EO_2MN_2TS_ANTI ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_15)
add_test(HMC17_CPU_TM_EO_2MN_3TS_ANTI ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_17)

# Test the RHMC
## LF: Leapfrog
## 2MN: 2min-norm
## XTS: X timescales
## CONF: starting from config COLD: starting from cold config
## ANTI: antiperiodic boundary conditions
### All these tests are with: EO preconditioning, fermact=rooted_stagg
## The options for the Rational Approximation files are given here to have a global path!!
set(APPROX_FILES "--approx_heatbath_file=${CMAKE_CURRENT_SOURCE_DIR}/Approx_Heatbath" "--approx_md_file=${CMAKE_CURRENT_SOURCE_DIR}/Approx_MD" "--approx_metropolis_file=${CMAKE_CURRENT_SOURCE_DIR}/Approx_Metropolis")
add_test(RHMC_01_CPU_LF_1TS_COLD        ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_1 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_01 ${APPROX_FILES})
add_test(RHMC_02_CPU_LF_1TS_COLD_ANTI   ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_2 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_02 ${APPROX_FILES})
add_test(RHMC_03_CPU_LF_1TS_CONF        ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_3 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_03 ${APPROX_FILES})
add_test(RHMC_04_CPU_LF_1TS_CONF_ANTI   ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_4 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_04 ${APPROX_FILES})
add_test(RHMC_05_CPU_LF_2TS_COLD        ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_1 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_05 ${APPROX_FILES})
add_test(RHMC_06_CPU_LF_2TS_COLD_ANTI   ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_2 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_06 ${APPROX_FILES})
add_test(RHMC_07_CPU_LF_2TS_CONF        ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_3 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_07 ${APPROX_FILES})
add_test(RHMC_08_CPU_LF_2TS_CONF_ANTI   ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_4 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_08 ${APPROX_FILES})
add_test(RHMC_09_CPU_2MN_1TS_COLD       ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_1 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_09 ${APPROX_FILES})
add_test(RHMC_10_CPU_2MN_1TS_COLD_ANTI  ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_2 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_10 ${APPROX_FILES})
add_test(RHMC_11_CPU_2MN_1TS_CONF       ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_3 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_11 ${APPROX_FILES})
add_test(RHMC_12_CPU_2MN_1TS_CONF_ANTI  ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_4 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_12 ${APPROX_FILES})
add_test(RHMC_13_CPU_2MN_2TS_COLD       ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_1 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_13 ${APPROX_FILES})
add_test(RHMC_14_CPU_2MN_2TS_COLD_ANTI  ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_2 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_14 ${APPROX_FILES})
add_test(RHMC_15_CPU_2MN_2TS_CONF       ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_3 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_15 ${APPROX_FILES})
add_test(RHMC_16_CPU_2MN_2TS_CONF_ANTI  ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testRhmc.py ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_ref_4 ${CMAKE_CURRENT_SOURCE_DIR}/rhmc_test_input_16 ${APPROX_FILES})


