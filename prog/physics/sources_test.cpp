/** @file
 * Tests for functions working with sources.
 */

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics
#include <boost/test/unit_test.hpp>

#include "sources.hpp"
#include <sstream>

void test_sources(std::string type, int num_sources)
{
	using namespace physics::lattices;

	std::stringstream tmp;
	tmp << "--num_sources=";
	tmp << num_sources;
	std::string n_sources_string = tmp.str();
	std::string sourcetype_string = std::string("--sourcetype=") + type;
	const char * _params[] = {"foo", n_sources_string.c_str(), sourcetype_string.c_str()};
	meta::Inputparameters params(3, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	auto sources = create_sources(system, prng);

	BOOST_REQUIRE_EQUAL(params.get_num_sources(), static_cast<const int>(sources.size()));

	release_spinorfields(sources);
}

BOOST_AUTO_TEST_CASE(sources)
{
	test_sources("point", 15);
	test_sources("volume", 2);
	test_sources("timeslice", 3);
	test_sources("zslice", 1);
}
