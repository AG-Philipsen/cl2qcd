message(STATUS "Building Tests")

set(Boost_USE_STATIC_LIBS   OFF)
find_package(Boost 1.34.0 REQUIRED COMPONENTS unit_test_framework)

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

add_executable(complex_mem_access
	complex_mem_access.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(complex_mem_access
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(MemoryAccessViaComplex complex_mem_access)
add_dependencies(complex_mem_access git_commit_id_header)

add_executable(su2su3extend
	su2su3extend.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(su2su3extend
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(Su2su3extend su2su3extend)
add_dependencies(su2su3extend git_commit_id_header)

add_executable(geometry
	geometry.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(geometry
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(geometry geometry)
add_dependencies(geometry git_commit_id_header)

add_executable(saxsbypz
	saxsbypz.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(saxsbypz
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(Saxsbypz saxsbypz)
add_dependencies(saxsbypz git_commit_id_header)

add_executable(m_gpu
	M_GPU.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../opencl_module_ran.cpp
	../opencl_module_spinors.cpp
	../opencl_module_fermions.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(m_gpu
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(M_gpu m_gpu)
add_dependencies(m_gpu git_commit_id_header)

add_executable(m_eo_gpu
	M_eo_GPU.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../opencl_module_ran.cpp
	../opencl_module_spinors.cpp
	../opencl_module_fermions.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(m_eo_gpu
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(M_eo_gpu m_eo_gpu)
add_dependencies(m_eo_gpu git_commit_id_header)

add_executable(f_gauge_gpu
	F_gauge_GPU.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../opencl_module_ran.cpp
	../opencl_module_spinors.cpp
	../opencl_module_fermions.cpp
	../opencl_module_hmc.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(f_gauge_gpu
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(F_gauge_gpu f_gauge_gpu)
add_dependencies(f_gauge_gpu git_commit_id_header)

add_executable(f_fermion_gpu
	Force_fermion_GPU.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../opencl_module_ran.cpp
	../opencl_module_spinors.cpp
	../opencl_module_fermions.cpp
	../opencl_module_hmc.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(f_fermion_gpu
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(F_fermion_gpu f_fermion_gpu)
add_dependencies(f_fermion_gpu git_commit_id_header)

add_executable(f_fermion_eo_gpu
	Force_fermion_eo_GPU.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../opencl_module_ran.cpp
	../opencl_module_spinors.cpp
	../opencl_module_fermions.cpp
	../opencl_module_hmc.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(f_fermion_eo_gpu
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(F_fermion_eo_gpu f_fermion_eo_gpu)
add_dependencies(f_fermion_eo_gpu git_commit_id_header)

add_executable(f_fermion_compare_non_eo_eo_gpu
	Force_fermion_compare_non_eo_eo_GPU.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../opencl_module_ran.cpp
	../opencl_module_spinors.cpp
	../opencl_module_fermions.cpp
	../opencl_module_hmc.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(f_fermion_compare_non_eo_eo_gpu
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(F_fermion_compare_non_eo_eo_gpu f_fermion_compare_non_eo_eo_gpu)
add_dependencies(f_fermion_compare_non_eo_eo_gpu git_commit_id_header)

add_executable(f_update_gpu
	Force_update_GPU.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../opencl_module_ran.cpp
	../opencl_module_spinors.cpp
	../opencl_module_fermions.cpp
	../opencl_module_hmc.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(f_update_gpu
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(F_update_gpu f_update_gpu)
add_dependencies(f_update_gpu git_commit_id_header)

add_executable(gf_update_gpu
	GF_update_GPU.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../opencl_module_ran.cpp
	../opencl_module_spinors.cpp
	../opencl_module_fermions.cpp
	../opencl_module_hmc.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(gf_update_gpu
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(Gf_update_gpu gf_update_gpu)
add_dependencies(gf_update_gpu git_commit_id_header)

add_executable(calc_staple
	Calc_staple_test.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(calc_staple
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(Calc_staple calc_staple)
add_dependencies(calc_staple git_commit_id_header)

add_executable(stout_smear
	Stout_Smear_test.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(stout_smear
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(Stout_smear stout_smear)
add_dependencies(stout_smear git_commit_id_header)

add_executable(localQ
	local_Q_test.cpp
	../host_geometry.cpp
	../host_operations_gaugefield.cpp
	../gaugefield_hybrid.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl_module.cpp
	../opencl_compiler.cpp
	../host_use_timer.cpp
	../host_writegaugefield.cpp
	../exceptions.cpp
	../logger.cpp)
target_link_libraries(localQ
	${OPENCL_LIBRARIES}
	${LIME_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_LIBRARIES}
	klepsydra
	crypto
	einhard)
add_test(LocalQ localQ)
add_dependencies(localQ git_commit_id_header)


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
add_test(HMC1_CPU_TM_EO_LF_1TS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_1)
add_test(HMC2_CPU_TM_NONEO_LF_1TS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_2)
add_test(HMC3_CPU_TM_EO_2MN_1TS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_3)
add_test(HMC5_CPU_TM_EO_LF_CG_1TS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_5)
add_test(HMC7_CPU_TM_EO_2MN_2TS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_7)
add_test(HMC9_CPU_TM_EO_2MN_3TS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHmc.py ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/hmc_test_input_9)
