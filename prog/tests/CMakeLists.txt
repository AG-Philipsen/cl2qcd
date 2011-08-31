message(STATUS "Building Tests")

set(Boost_USE_STATIC_LIBS   OFF)
find_package(Boost 1.34.0 REQUIRED COMPONENTS unit_test_framework)

find_package(PythonInterp REQUIRED)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

add_executable(openclCompiler opencl_compiler.cpp ../logger.cpp ../exceptions.cpp)
target_link_libraries(openclCompiler ${Boost_LIBRARIES} ${OPENCL_LIBRARIES} einhard)
add_test(OpenclCompiler openclCompiler)

add_executable(su2su3extend
	su2su3extend.cpp
	../host_geometry.cpp
	../host_operations_complex.cpp
	../host_operations_gaugefield.cpp
	../host_operations_su3matrix.cpp
	../host_operations_matrix.cpp
	../gaugefield.cpp
	../inputparameters.cpp
	../host_readgauge.cpp
	../host_random.cpp
	../opencl.cpp
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

# Test the inverter for correct inversion
add_test(Inverter1 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_1)
add_test(Inverter2 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_2)
add_test(Inverter3 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_3)
add_test(Inverter4 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_4)
add_test(Inverter5 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_5)
add_test(Inverter6 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_6)
add_test(Inverter7 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_7)
add_test(Inverter8 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testInverter.py ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_ref ${CMAKE_CURRENT_SOURCE_DIR}/inverter_test_input_8)
