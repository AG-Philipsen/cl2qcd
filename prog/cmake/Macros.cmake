# Utility functions for specifying the build of OpTiMaL
#
# Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

# Create a list of files just as you'd use set.
# Will automatically replace all filenames by absolute ones
macro(set_abs_paths DEST)
	set(FILES "${ARGV}")
	list(REMOVE_AT FILES 0) # get rid of dest value
	foreach(_FILE ${FILES})
		get_filename_component(_TMP ${_FILE} ABSOLUTE)
		list(APPEND ${DEST} ${_TMP})
	endforeach()
endmacro()

# Add adds one ore more modules to the library
# The module must declare a library called MODULE in the
# subdirectory called MODULE
macro(add_modules DEST MODULE)
	set(_MODULES "${ARGV}")
	list(REMOVE_AT _MODULES 0) # get rid of dest value
	foreach(_MODULE ${_MODULES})
		add_subdirectory(${_MODULE})
		set_property(GLOBAL APPEND PROPERTY DOC_SOURCE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/${_MODULE}")
	endforeach()
	target_link_libraries(${DEST} ${_MODULES})
endmacro()

# Add a unit test
# NAME name for the test to use in ctest
# EXE name for the executable to use
# The test is expectet to have only one source
# file called EXE.cpp and will be linked automatically
macro(add_unit_test NAME EXE)
	add_executable("${EXE}" "${EXE}.cpp")
	target_link_libraries("${EXE}" optimal)
	add_test("${NAME}" "${EXE}")
endmacro()
