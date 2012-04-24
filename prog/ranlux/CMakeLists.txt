add_library(ranlux ranlxs.c ranlxd.c)

add_executable(testlx testlx.c)
target_link_libraries(testlx ranlux)

