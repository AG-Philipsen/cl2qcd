# Utility functions for specifying the build of OpTiMaL
#
# Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>,
# 	        2014 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>,
#		     Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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
# The test is expected to have only one source
# file called EXE.cpp and will be linked automatically
macro(add_unit_test NAME EXE)
	add_executable("${EXE}" "${EXE}.cpp")
	target_link_libraries("${EXE}" optimal)
	add_test("${NAME}" "${EXE}" "${ARGN}")
endmacro()

# Add a unit test
# NAME name for the test to use in ctest
# EXE name for the executable to use
# LIB name of the library to link to
# The test is expected to have only one source
# file called EXE.cpp and will be linked automatically
# ARGN includes all the arguments passed to the macro beyond EXE
#      and they should all be libraries needed by the executable
macro(add_unit_test_withLibraries NAME EXE)
	add_executable("${EXE}" "${EXE}.cpp")
	target_link_libraries("${EXE}" ${ARGN} ${Boost_LIBRARIES})
	add_test("${NAME}" "${EXE}")
endmacro()
