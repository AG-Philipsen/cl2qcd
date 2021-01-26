# Utility functions for specifying the build of CL2QCD
#
# Copyright (c) 2012,2013 Matthias Bach
# Copyright (c) 2013,2015,2018,2021 Alessandro Sciarra
# Copyright (c) 2014 Christopher Pinke
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.


# Look for OpenCL cpp bindings in OpenCL_INCLUDE_DIRS and relatively to
# libraries (in case the user uses AMD SDK which places its headers in
# a non-standard way).
function(find_OpenCL_cpp_bindings)
    if(OpenCL_HAS_CPP_BINDINGS)
        message(STATUS "Bindings already found, no further research")
        return()
    endif(OpenCL_HAS_CPP_BINDINGS)
    get_filename_component(OpenCL_LIB_DIRECTORY ${OpenCL_LIBRARIES} PATH)
    get_filename_component(OpenCL_CPP_BINDINGS_CANDIDATE ${OpenCL_LIB_DIRECTORY}/../../include ABSOLUTE)
    find_path(OpenCL_CPP_INCLUDE_DIRS CL/cl.hpp PATHS ${OpenCL_INCLUDE_DIRS} ${OpenCL_CPP_BINDINGS_CANDIDATE})
    if(OpenCL_CPP_INCLUDE_DIRS)
        set(OpenCL_HAS_CPP_BINDINGS TRUE)
        list(APPEND OpenCL_INCLUDE_DIRS ${OpenCL_CPP_INCLUDE_DIRS})
        # This is often the same, so clean up
        list(REMOVE_DUPLICATES OpenCL_INCLUDE_DIRS)
    endif(OpenCL_CPP_INCLUDE_DIRS)
    unset(OpenCL_CPP_INCLUDE_DIRS CACHE)
endfunction()

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

# Add adds one ore more modules to the library.
#
# This function has to be called with the keywords followed by their value(s)
#   TO followed by the name of the library
#   MODULES followed by the names of the modules to add to the library
#   INTERFACE (without value) if the library is an interface one
#
# The module must declare a library called MODULE in the subdirectory called MODULE
# (here MODULE is the name passed after the MODULES keyword)
function(add_modules)
    set(options INTERFACE)
    set(oneValueArgs TO)
    set(multiValueArgs MODULES)
    cmake_parse_arguments(add_modules "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    #Abort if this function was incorrectly called
    if(NOT "${add_modules_UNPARSED_ARGUMENTS}" STREQUAL "")
        message(FATAL_ERROR "CMake function \"add_modules\" called incorrectly (\"${add_modules_UNPARSED_ARGUMENTS}\" unrecognised), aborting!")
    endif()
    if("${add_modules_TO}" STREQUAL "" OR "${add_modules_MODULES}" STREQUAL "")
        message(FATAL_ERROR "CMake function \"add_modules\" called incorrectly (TO or MODULES missing), aborting!")
    endif()

    # Add the sub-directories and link as requested
    foreach(_MODULE ${add_modules_MODULES})
        add_subdirectory(${_MODULE})
        set_property(GLOBAL APPEND PROPERTY DOC_SOURCE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/${_MODULE}")
    endforeach()
    if(add_modules_INTERFACE)
        target_link_libraries(${add_modules_TO} INTERFACE ${add_modules_MODULES})
    else()
        target_link_libraries(${add_modules_TO} ${add_modules_MODULES})
    endif()
endfunction()

# According to
#   https://cmake.org/cmake/help/v3.5/command/cmake_parse_arguments.html#command:cmake_parse_arguments
# it is relatively simple to create a user-defined function which
# is able to deal with the task of creating test executable linked
# to one or more specified libraries and, optionally, running it with
# optionally specified command line options.
#
# This function has to be called with the keywords followed by their value(s)
#   EXECUTABLE followed by the name of the executable
#   SOURCE_FILES followed by the file(s) to be compiled
#   LIBRARIES followed by the library(ies) needed for the linking
#   NAME followed by the ctest label for the test
#   COMMAND_LINE_OPTIONS followed by command line options for the test

# Additional mutually-exclusive keywords to be specified without value are
#   CREATE_ONLY to decide not to add a test case to ctest
#   ADD_ONLY to decide not to create the executable
#
# The idea behind this function, is to use as ctest name the folder path to
# the _test.cpp file ending with the test file without the '_test.cpp' suffix.
# In this way, in most of the cases, the EXECUTABLE and the SOURCE_FILES can
# be automatically deduced (see NOTE below). Of course, the user has always the
# possibility to specify all to achieve some particular goal.
#
# NOTE: The SOURCE_FILES, only if omitted, is supposed to be a single source file which
#       is either the EXECUTABLE with the .cpp extension added if the EXECUTABLE was specifed
#       or which is obtained from the NAME by
#         1) removing any occurence of '_REC12', '_CPU', '_GPU' from it and
#         2) removing anything before the last '/', slash included.
#
#       If the EXECUTABLE key was not specified, it is deduced from NAME by
#         1) removing any occurence of '_REC12', '_CPU', '_GPU' from it and
#         2) replacing '/' with '_'.
#
function(add_unit_test)
    set(options CREATE_ONLY ADD_ONLY)
    set(oneValueArgs NAME EXECUTABLE)
    set(multiValueArgs COMMAND_LINE_OPTIONS LIBRARIES SOURCE_FILES)
    cmake_parse_arguments(add_unit_test "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    #Abort if this function was incorrectly called
    if(NOT "${add_unit_test_UNPARSED_ARGUMENTS}" STREQUAL "")
        message( FATAL_ERROR "CMake function \"add_unit_test\" called incorrectly (\"${add_unit_test_UNPARSED_ARGUMENTS}\" unrecognised), aborting!" )
    endif(NOT "${add_unit_test_UNPARSED_ARGUMENTS}" STREQUAL "")
    if(add_unit_test_ADD_ONLY AND add_unit_test_CREATE_ONLY)
        message( FATAL_ERROR "CMake function \"add_unit_test\" called incorrectly, ADD_ONLY and CREATE_ONLY specified together!" )
    endif(add_unit_test_ADD_ONLY AND add_unit_test_CREATE_ONLY)

    #If the executable has to be produced, the SOURCE_FILES must be deduced if not provided
    if(NOT add_unit_test_ADD_ONLY)
        if( "${add_unit_test_SOURCE_FILES}" STREQUAL "" )
            if( NOT "${add_unit_test_EXECUTABLE}" STREQUAL "" )
                string(CONCAT add_unit_test_SOURCE_FILES "${add_unit_test_EXECUTABLE}" ".cpp")
            else( NOT "${add_unit_test_EXECUTABLE}" STREQUAL "" )
                if("${add_unit_test_NAME}" STREQUAL "" )
                    message( FATAL_ERROR "SOURCE_FILES not deducible since NAME not givent to CMake function \"add_unit_test\", aborting!" )
                else("${add_unit_test_NAME}" STREQUAL "" )
                    string(REGEX REPLACE "_(REC12|CPU|GPU)" "" add_unit_test_SOURCE_FILES "${add_unit_test_NAME}")
                    string(REGEX REPLACE ".*/" "" add_unit_test_SOURCE_FILES "${add_unit_test_SOURCE_FILES}")
                    string(APPEND add_unit_test_SOURCE_FILES "_test.cpp")
                endif()
            endif()
        endif()
    endif()

    #The EXECUTABLE must be always available and it must be deduced if not provided
    if( "${add_unit_test_EXECUTABLE}" STREQUAL "" )
        if( "${add_unit_test_NAME}" STREQUAL "" )
            message( FATAL_ERROR "EXECUTABLE name not specified and not deducible in CMake function \"add_unit_test\", aborting!" )
        else( "${add_unit_test_NAME}" STREQUAL "" )
            string(REGEX REPLACE "_(REC12|CPU|GPU)" "" add_unit_test_EXECUTABLE "${add_unit_test_NAME}")
            string(REPLACE "/" "_" add_unit_test_EXECUTABLE "${add_unit_test_EXECUTABLE}")
            string(APPEND add_unit_test_EXECUTABLE "_test")
        endif()
    endif()

    #Final steps: create executable and/or ctest case
    if(NOT add_unit_test_ADD_ONLY)
        add_executable("${add_unit_test_EXECUTABLE}" ${add_unit_test_SOURCE_FILES})
        target_link_libraries("${add_unit_test_EXECUTABLE}" ${add_unit_test_LIBRARIES} ${Boost_LIBRARIES})
    endif()
    if(NOT add_unit_test_CREATE_ONLY)
        if( "${add_unit_test_NAME}" STREQUAL "" )
            message( FATAL_ERROR "NAME <name> not specified in CMake function \"add_unit_test\", but asked to add ctest test!" )
        endif()
        add_test("${add_unit_test_NAME}" "${add_unit_test_EXECUTABLE}" ${add_unit_test_COMMAND_LINE_OPTIONS})
    endif()
endfunction()
