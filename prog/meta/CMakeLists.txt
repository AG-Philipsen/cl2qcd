# Build specifications for the contents of the Meta package
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

#
# Collect all code in this part of the source tree
#

add_library(meta
	inputparameters.cpp
	util.cpp
	util_io.cpp
	counter.cpp
	config_file_normalizer.cpp
	type_ops.cpp
	../logger.cpp
)
add_dependencies(meta git_commit_id_header)

#
# Definition of tets
#

configure_file(test_input_1 test_input_1 COPYONLY)
configure_file(test_input_2 test_input_2 COPYONLY)
configure_file(test_input_aliases test_input_aliases COPYONLY)

add_unit_test(meta::Inputparameters inputparameters_test)
add_unit_test(meta::ConfigFileNormalizer config_file_normalizer_test)
add_unit_test(meta::TypeOps type_ops_test)
