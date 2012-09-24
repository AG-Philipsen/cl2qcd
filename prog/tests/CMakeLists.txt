message(STATUS "Building Tests")

find_package(PythonInterp REQUIRED)

# default random number seed
configure_file( ${CMAKE_SOURCE_DIR}/rand_seeds rand_seeds COPYONLY )

# default random number seed
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/conf.00200 conf.00200 COPYONLY )

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

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

add_executable(gaugefield_observables gaugefield_observables.cpp)
target_link_libraries(gaugefield_observables
	optimal
	${Boost_LIBRARIES}
)
add_test(Gaugeobservables gaugefield_observables)

add_executable(calc_staple Calc_staple_test.cpp)
target_link_libraries(calc_staple
	optimal
	${Boost_LIBRARIES}
)
add_test(Calc_staple_CPU calc_staple ${CMAKE_CURRENT_SOURCE_DIR}/staple_input_1 --use_gpu=false)
add_test(Calc_staple_CPU_REC12 calc_staple ${CMAKE_CURRENT_SOURCE_DIR}/staple_input_rec12_1 --use_gpu=false)
add_test(Calc_staple_GPU calc_staple ${CMAKE_CURRENT_SOURCE_DIR}/staple_input_1 --use_gpu=true)
add_test(Calc_staple_GPU_REC12 calc_staple ${CMAKE_CURRENT_SOURCE_DIR}/staple_input_rec12_1 --use_gpu=true)

add_executable(f_gauge F_gauge_test.cpp)
target_link_libraries(f_gauge
	optimal
	${Boost_LIBRARIES}
)
add_test(F_gauge_CPU f_gauge ${CMAKE_CURRENT_SOURCE_DIR}/f_gauge_input_1 --use_gpu=false)
add_test(F_gauge_CPU_REC12 f_gauge ${CMAKE_CURRENT_SOURCE_DIR}/f_gauge_input_rec12 --use_gpu=false)
add_test(F_gauge_GPU f_gauge ${CMAKE_CURRENT_SOURCE_DIR}/f_gauge_input_1 --use_gpu=true)
add_test(F_gauge_GPU_REC12 f_gauge ${CMAKE_CURRENT_SOURCE_DIR}/f_gauge_input_rec12 --use_gpu=true)

add_executable(rectangles Rectangles_test.cpp)
target_link_libraries(rectangles
	optimal
	${Boost_LIBRARIES}
)
add_test(Rectangles_CPU rectangles  ${CMAKE_CURRENT_SOURCE_DIR}/rectangles_input_1 --use_gpu=false)
add_test(Rectangles_CPU_REC12 rectangles  ${CMAKE_CURRENT_SOURCE_DIR}/rectangles_input_rec12_1 --use_gpu=false)
add_test(Rectangles_GPU rectangles  ${CMAKE_CURRENT_SOURCE_DIR}/rectangles_input_1 --use_gpu=true)
add_test(Rectangles_GPU_REC12 rectangles  ${CMAKE_CURRENT_SOURCE_DIR}/rectangles_input_rec12_1 --use_gpu=true)

add_executable(f_gauge_tlsym F_gauge_tlsym_test.cpp)
target_link_libraries(f_gauge_tlsym
	optimal
	${Boost_LIBRARIES}
)
add_test(F_gauge_tlsym_CPU f_gauge_tlsym ${CMAKE_CURRENT_SOURCE_DIR}/f_gauge_tlsym_input_1 --use_gpu=false)
add_test(F_gauge_tlsym_CPU_REC12 f_gauge_tlsym ${CMAKE_CURRENT_SOURCE_DIR}/f_gauge_tlsym_input_rec12_1 --use_gpu=false)
add_test(F_gauge_tlsym_GPU f_gauge_tlsym ${CMAKE_CURRENT_SOURCE_DIR}/f_gauge_tlsym_input_1 --use_gpu=true)
add_test(F_gauge_tlsym_GPU_REC12 f_gauge_tlsym ${CMAKE_CURRENT_SOURCE_DIR}/f_gauge_tlsym_input_rec12_1 --use_gpu=true)

add_executable(localQ local_Q_test.cpp)
target_link_libraries(localQ
	optimal
	${Boost_LIBRARIES}
)
add_test(LocalQ_CPU localQ ${CMAKE_CURRENT_SOURCE_DIR}/local_Q_input_1 --use_gpu=false)
add_test(LocalQ_CPU_REC12 localQ ${CMAKE_CURRENT_SOURCE_DIR}/local_Q_input_rec12_1 --use_gpu=false)
add_test(LocalQ_GPU localQ ${CMAKE_CURRENT_SOURCE_DIR}/local_Q_input_1 --use_gpu=true)
add_test(LocalQ_GPU_REC12 localQ ${CMAKE_CURRENT_SOURCE_DIR}/local_Q_input_rec12_1 --use_gpu=true)

add_executable(gf_update GF_update_test.cpp)
target_link_libraries(gf_update
	optimal
	${Boost_LIBRARIES}
)
add_test(Gf_update_CPU gf_update ${CMAKE_CURRENT_SOURCE_DIR}/gf_update_input_1 --use_gpu=false)
add_test(Gf_update_CPU_REC12 gf_update ${CMAKE_CURRENT_SOURCE_DIR}/gf_update_input_rec12_1 --use_gpu=false)
add_test(Gf_update_GPU gf_update ${CMAKE_CURRENT_SOURCE_DIR}/gf_update_input_1 --use_gpu=true)
add_test(Gf_update_GPU_REC12 gf_update ${CMAKE_CURRENT_SOURCE_DIR}/gf_update_input_rec12_1 --use_gpu=true)

add_executable(stout_smear Stout_Smear_test.cpp)
target_link_libraries(stout_smear
	optimal
	${Boost_LIBRARIES}
)
add_test(Stout_smear_CPU stout_smear ${CMAKE_CURRENT_SOURCE_DIR}/stout_smear_input_1 --use_gpu=false)

add_executable(f_update F_update_test.cpp)
target_link_libraries(f_update
	optimal
	${Boost_LIBRARIES}
)
add_test(F_update_CPU f_update ${CMAKE_CURRENT_SOURCE_DIR}/f_update_input_1  --use_gpu=false)
add_test(F_update_GPU f_update ${CMAKE_CURRENT_SOURCE_DIR}/f_update_input_1  --use_gpu=true)

add_executable(m_tm_plus M_tm_plus_test.cpp)
target_link_libraries(m_tm_plus
	optimal
	${Boost_LIBRARIES}
)
add_test(M_tm_plus_CPU m_tm_plus  ${CMAKE_CURRENT_SOURCE_DIR}/m_tm_plus_input_1 --use_gpu=false)
add_test(M_tm_plus_CPU_REC12 m_tm_plus  ${CMAKE_CURRENT_SOURCE_DIR}/m_tm_plus_input_rec12_1 --use_gpu=false)
add_test(M_tm_plus_GPU m_tm_plus  ${CMAKE_CURRENT_SOURCE_DIR}/m_tm_plus_input_1 --use_gpu=true)
add_test(M_tm_plus_GPU_REC12 m_tm_plus  ${CMAKE_CURRENT_SOURCE_DIR}/m_tm_plus_input_rec12_1 --use_gpu=true)

add_executable(dslash_eo dslash_eo_test.cpp)
target_link_libraries(dslash_eo
	optimal
	${Boost_LIBRARIES}
)
add_test(Dslash_eo_CPU_1 dslash_eo ${CMAKE_CURRENT_SOURCE_DIR}/dslash_eo_input_1 --use_gpu=false)
add_test(Dslash_eo_CPU_REC12_1 dslash_eo ${CMAKE_CURRENT_SOURCE_DIR}/dslash_eo_input_rec12_1 --use_gpu=false)
add_test(Dslash_eo_GPU_1 dslash_eo ${CMAKE_CURRENT_SOURCE_DIR}/dslash_eo_input_1 --use_gpu=true)
add_test(Dslash_eo_GPU_REC12_1 dslash_eo ${CMAKE_CURRENT_SOURCE_DIR}/dslash_eo_input_rec12_1 --use_gpu=true)
add_test(Dslash_eo_CPU_2 dslash_eo ${CMAKE_CURRENT_SOURCE_DIR}/dslash_eo_input_2 --use_gpu=false)
add_test(Dslash_eo_CPU_REC12_2 dslash_eo ${CMAKE_CURRENT_SOURCE_DIR}/dslash_eo_input_rec12_2 --use_gpu=false)
add_test(Dslash_eo_GPU_2 dslash_eo ${CMAKE_CURRENT_SOURCE_DIR}/dslash_eo_input_2 --use_gpu=true)
add_test(Dslash_eo_GPU_REC12_2 dslash_eo ${CMAKE_CURRENT_SOURCE_DIR}/dslash_eo_input_rec12_2 --use_gpu=true)

add_executable(m_tm_plus_compare_noneo_eo M_tm_plus_compare_noneo_eo.cpp)
target_link_libraries(m_tm_plus_compare_noneo_eo
	optimal
	${Boost_LIBRARIES}
)
add_test(M_tm_plus_compare_noneo_eo_CPU m_tm_plus_compare_noneo_eo ${CMAKE_CURRENT_SOURCE_DIR}/m_tm_plus_compare_noneo_eo_input_1 --use_gpu=false)
add_test(M_tm_plus_compare_noneo_eo_CPU_REC12 m_tm_plus_compare_noneo_eo ${CMAKE_CURRENT_SOURCE_DIR}/m_tm_plus_compare_noneo_eo_input_rec12_1 --use_gpu=false)
add_test(M_tm_plus_compare_noneo_eo_GPU m_tm_plus_compare_noneo_eo ${CMAKE_CURRENT_SOURCE_DIR}/m_tm_plus_compare_noneo_eo_input_1 --use_gpu=true)
add_test(M_tm_plus_compare_noneo_eo_GPU_REC12 m_tm_plus_compare_noneo_eo ${CMAKE_CURRENT_SOURCE_DIR}/m_tm_plus_compare_noneo_eo_input_rec12_1 --use_gpu=true)

add_executable(f_fermion Force_fermion_GPU.cpp)
target_link_libraries(f_fermion
	optimal
	${Boost_LIBRARIES}
)
add_test(F_fermion_CPU f_fermion ${CMAKE_CURRENT_SOURCE_DIR}/f_fermion_input_1 --use_gpu=false)
add_test(F_fermion_CPU_REC12 f_fermion ${CMAKE_CURRENT_SOURCE_DIR}/f_fermion_input_rec12_1 --use_gpu=false)
add_test(F_fermion_GPU f_fermion ${CMAKE_CURRENT_SOURCE_DIR}/f_fermion_input_1 --use_gpu=true)
add_test(F_fermion_GPU_REC12 f_fermion ${CMAKE_CURRENT_SOURCE_DIR}/f_fermion_input_rec12_1 --use_gpu=true)

add_executable(f_fermion_eo_gpu Force_fermion_eo_GPU.cpp)
target_link_libraries(f_fermion_eo_gpu
	optimal
	${Boost_LIBRARIES}
)
add_test(F_fermion_eo_gpu f_fermion_eo_gpu)

add_executable(f_fermion_compare_non_eo_eo_gpu Force_fermion_compare_non_eo_eo_GPU.cpp)
target_link_libraries(f_fermion_compare_non_eo_eo_gpu
	optimal
	${Boost_LIBRARIES}
)
add_test(F_fermion_compare_non_eo_eo_gpu f_fermion_compare_non_eo_eo_gpu)

# ignore some tests
configure_file(CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake @ONLY)

# Test the inverter for correct inversion
add_test(Inverter1_CPU_TM_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1)
add_test(Inverter3_GPU_TM_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_3)
add_test(Inverter5_CPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_5)
add_test(Inverter7_GPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_7)
add_test(Inverter9_CPU_TM_EO_CG ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref_cg ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_9)
add_test(Inverter11_GPU_TM_EO_CG ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref_cg ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_11)
add_test(Inverter13_CPU_TM_CG ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref_cg ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_13)
add_test(Inverter15_GPU_TM_CG ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref_cg ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_15)

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


