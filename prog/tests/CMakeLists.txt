message(STATUS "Building Tests")

set(Boost_USE_STATIC_LIBS   OFF)
find_package(Boost 1.34.0 REQUIRED COMPONENTS unit_test_framework)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

add_executable(openclCompiler opencl_compiler.cpp ../logger.cpp)
target_link_libraries(openclCompiler ${Boost_LIBRARIES} ${OPENCL_LIBRARIES} einhard)
add_test(OpenclCompiler openclCompiler)

