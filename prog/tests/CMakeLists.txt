message(STATUS "Building Tests")

find_package(PythonInterp REQUIRED)

# default random number seed
configure_file( ${CMAKE_SOURCE_DIR}/rand_seeds rand_seeds COPYONLY )

# default random number seed
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/conf.00200 conf.00200 COPYONLY )

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

add_executable(Boost_test boost_test.cpp)
target_link_libraries(Boost_test
	${Boost_LIBRARIES}
)
add_test(Boost_ARGC_1 Boost_test --run_test=BOOST_ARGC_1)
add_test(Boost_ARGC_2 Boost_test ${CMAKE_CURRENT_SOURCE_DIR} --run_test=BOOST_ARGC_2)
add_test(Boost_ARGV Boost_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=false   --run_test=BOOST_ARGV )

add_executable(openclCompiler opencl_compiler.cpp ../opencl_compiler.cpp ../logger.cpp ../exceptions.cpp)
target_link_libraries(openclCompiler ${Boost_LIBRARIES} ${OPENCL_LIBRARIES} crypto einhard)
add_test(OpenclCompiler openclCompiler)

add_executable(cryptoLibTest crypto.cpp)
target_link_libraries(cryptoLibTest ${Boost_LIBRARIES} crypto)
add_test(Crypto cryptoLibTest)

add_executable(complex_mem_access complex_mem_access.cpp)
target_link_libraries(complex_mem_access
	optimal
	${Boost_LIBRARIES}
)
add_test(MemoryAccessViaComplex complex_mem_access)

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

add_executable(saxsbypz saxsbypz.cpp)
target_link_libraries(saxsbypz
	optimal
	${Boost_LIBRARIES}
)
add_test(Saxsbypz saxsbypz)

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

add_executable(opencl_module_test opencl_module_test.cpp)
target_link_libraries(opencl_module_test
	optimal
	${Boost_LIBRARIES}
)
add_test(Opencl_module_CPU       opencl_module_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=false)
add_test(Opencl_module_CPU_REC12 opencl_module_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=true)
add_test(Opencl_module_GPU       opencl_module_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=false)
add_test(Opencl_module_GPU_REC12 opencl_module_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=true)

add_executable(opencl_module_fermions_test opencl_module_fermions_test.cpp)
target_link_libraries(opencl_module_fermions_test
	optimal
	${Boost_LIBRARIES}
)
add_test(Opencl_module_fermions_CPU       opencl_module_fermions_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=false)
add_test(Opencl_module_fermions_CPU_REC12 opencl_module_fermions_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=true)
add_test(Opencl_module_fermions_GPU       opencl_module_fermions_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=false)
add_test(Opencl_module_fermions_GPU_REC12 opencl_module_fermions_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=true)

add_executable(opencl_module_gaugemomentum_test opencl_module_gaugemomentum_test.cpp)
target_link_libraries(opencl_module_gaugemomentum_test
	optimal
	${Boost_LIBRARIES}
)
add_test(Opencl_module_gaugemomentum_CPU       opencl_module_gaugemomentum_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=false)
add_test(Opencl_module_gaugemomentum_CPU_REC12 opencl_module_gaugemomentum_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=true)
add_test(Opencl_module_gaugemomentum_GPU       opencl_module_gaugemomentum_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=false)
add_test(Opencl_module_gaugemomentum_GPU_REC12 opencl_module_gaugemomentum_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=true)

add_executable(opencl_module_molecular_dynamics_test opencl_module_molecular_dynamics_test.cpp)
target_link_libraries(opencl_module_molecular_dynamics_test
	optimal
	${Boost_LIBRARIES}
)
add_test(Opencl_module_molecular_dynamics_CPU       opencl_module_molecular_dynamics_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=false)
add_test(Opencl_module_molecular_dynamics_CPU_REC12 opencl_module_molecular_dynamics_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=false --use_rec12=true)
add_test(Opencl_module_molecular_dynamics_GPU       opencl_module_molecular_dynamics_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=false)
add_test(Opencl_module_molecular_dynamics_GPU_REC12 opencl_module_molecular_dynamics_test ${CMAKE_CURRENT_SOURCE_DIR} --use_gpu=true --use_rec12=true)

# ignore some tests
configure_file(CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake @ONLY)

# Test the inverter for correct inversion
## even testnumbers correspond to tests using even-odd-precond., 
## odd numbers to those who don`t
## tests on a cold config
add_test(Inverter1_CPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_rec12=false --use_eo=false)
add_test(Inverter2_CPU_TM_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_rec12=false --use_eo=true)
add_test(Inverter3_CPU_TM_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_rec12=true --use_eo=false)
add_test(Inverter4_CPU_TM_EO_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=0 --use_rec12=true --use_eo=true)
add_test(Inverter5_GPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_rec12=false --use_eo=false)
add_test(Inverter6_GPU_TM_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_rec12=false --use_eo=true)
add_test(Inverter7_GPU_TM_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_rec12=true --use_eo=false)
add_test(Inverter8_GPU_TM_EO_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1 --use_gpu=1 --use_rec12=true --use_eo=true)
## tests on a given config
## the config does not contain a valid checksum, the resulting error will be ignored
add_test(Inverter9_CPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0 --use_rec12=false --use_eo=false --ignore_checksum_errors=true )
add_test(Inverter10_CPU_TM_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=0
 --use_rec12=false --use_eo=true --ignore_checksum_errors=true)
add_test(Inverter11_CPU_TM_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 
--use_gpu=0 --use_rec12=true --use_eo=false  --ignore_checksum_errors=true)
add_test(Inverter12_CPU_TM_EO_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 
--use_gpu=0 --use_rec12=true --use_eo=true  --ignore_checksum_errors=true)
add_test(Inverter13_GPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 --use_gpu=1 
--use_rec12=false --use_eo=false  --ignore_checksum_errors=true)
add_test(Inverter14_GPU_TM_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 
--use_gpu=1 --use_rec12=false --use_eo=true  --ignore_checksum_errors=true)
add_test(Inverter15_GPU_TM_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 
--use_gpu=1 --use_rec12=true --use_eo=false  --ignore_checksum_errors=true)  
add_test(Inverter16_GPU_TM_EO_REC12 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref2 ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2 
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


# Test the heatbath
add_test(Heatbath1_CPU_DP ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHeatbath.py ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_input_1)
add_test(Heatbath3_GPU_DP ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHeatbath.py ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_input_3)
add_test(Heatbath5_CPU_SP ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHeatbath.py ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_input_5)
add_test(Heatbath7_GPU_SP ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHeatbath.py ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_input_7)

# Test the HMC
## odd testnumbers correspond to tests using even-odd-precond., 
## even numbers to those who don`t
## LF: Leapfrog
## 2MN: 2mn
## XTS: X timescales
## CONF: starting from config
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


