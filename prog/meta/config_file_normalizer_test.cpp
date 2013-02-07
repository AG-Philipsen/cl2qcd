/** @file
 * Testcases for the meta::ConfigFileNormalizer class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "config_file_normalizer.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE meta::ConfigFileNormalizer
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <stdexcept>

using meta::ConfigFileNormalizer;

static void check(const ConfigFileNormalizer& normalizer, const std::string& input_filename, const std::string reference_value);

BOOST_AUTO_TEST_CASE(invalid_file)
{
	const ConfigFileNormalizer cand;
	BOOST_CHECK_THROW(cand("foo.bla"), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_noop)
{
	const ConfigFileNormalizer cand;

	check(cand, "config_file_normalizer_test_empty", "");
	check(cand, "config_file_normalizer_test_input1", "foo\n"
                                                      "bla=blub\n"
                                                      "foo=bla\n"
                                                      "old=new\n"
                                                      "new=old\n"
                                                      "  old=new\n"
                                                      "  new=old\n"
                                                      "other = more\n");
}

BOOST_AUTO_TEST_CASE(simple)
{
	ConfigFileNormalizer cand;
	cand.add_alias("foo", "blubber");
	cand.add_alias("other", "another");

	check(cand, "config_file_normalizer_test_empty", "");
	check(cand, "config_file_normalizer_test_input1", "blubber\n"
                                                      "bla=blub\n"
                                                      "blubber=bla\n"
                                                      "old=new\n"
                                                      "new=old\n"
                                                      "  old=new\n"
                                                      "  new=old\n"
                                                      "another = more\n");
}

BOOST_AUTO_TEST_CASE(match_as_val)
{
	ConfigFileNormalizer cand;
	cand.add_alias("old", "legacy");

	check(cand, "config_file_normalizer_test_empty", "");
	check(cand, "config_file_normalizer_test_input1", "foo\n"
                                                      "bla=blub\n"
                                                      "foo=bla\n"
                                                      "legacy=new\n"
                                                      "new=old\n"
                                                      "  legacy=new\n"
                                                      "  new=old\n"
                                                      "other = more\n");
}

BOOST_AUTO_TEST_CASE(dont_match_partial)
{
	ConfigFileNormalizer cand;
	cand.add_alias("ol", "ups");

	check(cand, "config_file_normalizer_test_empty", "");
	check(cand, "config_file_normalizer_test_input1", "foo\n"
                                                      "bla=blub\n"
                                                      "foo=bla\n"
                                                      "old=new\n"
                                                      "new=old\n"
                                                      "  old=new\n"
                                                      "  new=old\n"
                                                      "other = more\n");
}

static void check(const ConfigFileNormalizer& normalizer, const std::string& input_filename, const std::string reference_value)
{
	const std::string full_input_path = std::string(SOURCEDIR) + "/meta/" + input_filename;
	BOOST_CHECK_EQUAL(normalizer(full_input_path), reference_value);

	std::ifstream input_file(full_input_path.c_str());
	BOOST_CHECK_EQUAL(normalizer(input_file), reference_value);
}
