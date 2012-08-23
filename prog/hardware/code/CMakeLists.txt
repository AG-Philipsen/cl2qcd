# Build specifications for the contents of the code package
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

#
# Collect all code in this part of the source tree
#

add_library(code
	opencl_module.cpp
)

#
# Definition of tets
#

add_unit_test(hardware::code::OpenclModule opencl_module_test)

