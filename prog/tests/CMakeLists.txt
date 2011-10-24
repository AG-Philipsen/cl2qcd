message(STATUS "Building Tests")

set(Boost_USE_STATIC_LIBS   OFF)
find_package(Boost 1.34.0 REQUIRED COMPONENTS unit_test_framework)

find_package(PythonInterp REQUIRED)

# default random number seed
configure_file( ${CMAKE_SOURCE_DIR}/rand_seeds rand_seeds COPYONLY )

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

add_executable(openclCompiler opencl_compiler.cpp ../opencl_compiler.cpp ../logger.cpp ../exceptions.cpp)
target_link_libraries(openclCompiler ${Boost_LIBRARIES} ${OPENCL_LIBRARIES} einhard)
add_test(OpenclCompiler openclCompiler)

add_executable(su2su3extend
	su2su3extend.cpp
	../host_geometry.cpp
	../host_operations_complex.cpp
	../host_operations_gaugefield.cpp
	../host_operations_su3matrix.cpp
	../host_operations_matrix.cpp
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
	einhard)
add_test(Su2su3extend su2su3extend)

add_executable(geometry
	geometry.cpp
	../host_geometry.cpp
	../host_operations_complex.cpp
	../host_operations_gaugefield.cpp
	../host_operations_su3matrix.cpp
	../host_operations_matrix.cpp
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
	einhard)
add_test(geometry geometry)

add_executable(saxsbypz
	saxsbypz.cpp
	../host_geometry.cpp
	../host_operations_complex.cpp
	../host_operations_gaugefield.cpp
	../host_operations_su3matrix.cpp
	../host_operations_matrix.cpp
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
	einhard)
add_test(Saxsbypz saxsbypz)

add_executable(m_gpu
	M_GPU.cpp
	../host_geometry.cpp
	../host_operations_complex.cpp
	../host_operations_gaugefield.cpp
	../host_operations_su3matrix.cpp
	../host_operations_matrix.cpp
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
	einhard)
add_test(M_gpu m_gpu)

add_executable(m_eo_gpu
	M_eo_GPU.cpp
	../host_geometry.cpp
	../host_operations_complex.cpp
	../host_operations_gaugefield.cpp
	../host_operations_su3matrix.cpp
	../host_operations_matrix.cpp
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
	einhard)
add_test(M_eo_gpu m_eo_gpu)

# ignore some tests
configure_file(CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake @ONLY)

# Test the inverter for correct inversion
add_test(Inverter1_CPU_TM_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1)
add_test(Inverter3_GPU_TM_EO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_3)
add_test(Inverter5_CPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_5)
add_test(Inverter7_GPU_TM ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_7)

add_test(Heatbath1_CPU_DP ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHeatbath.py ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_input_1)
add_test(Heatbath3_GPU_DP ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHeatbath.py ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_input_3)
add_test(Heatbath5_CPU_SP ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHeatbath.py ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_input_5)
add_test(Heatbath7_GPU_SP ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testHeatbath.py ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/heatbath_test_input_7)
