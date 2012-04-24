add_library(ranlux ranlxs.c ranlxd.c)

add_executable(testlx testlx.c)
target_link_libraries(testlx ranlux)

add_test(NAME Ranlux COMMAND testlx)
set_tests_properties(Ranlux PROPERTIES PASS_REGULAR_EXPRESSION "ranlxs and ranlxd work correctly on this machine")
