/** @file
 * Ranlux PRNG unit test
 *
 * Copyright 2012, 2013, 2015 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "prng.hpp"
#include "../hardware/code/testUtilities.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::PRNG
#include <boost/test/unit_test.hpp>

using namespace physics;

void verifyBothBuffersAreEquallyLarge( const hardware::buffers::PRNGBuffer* buf1, const hardware::buffers::PRNGBuffer* buf2 )
{
	size_t const buf1_bytes = buf1->get_bytes();
	size_t const buf2_bytes = buf2->get_bytes();
	BOOST_REQUIRE_EQUAL(buf1_bytes, buf2_bytes);
}

hardware::buffers::PRNGBuffer::prng_state_t * createPrngState( const size_t buffer_size)
{
	return new hardware::buffers::PRNGBuffer::prng_state_t[buffer_size];
}

void verifyBuffersAreDifferent( const hardware::buffers::PRNGBuffer* buf1, const hardware::buffers::PRNGBuffer* buf2 )
{
	size_t const buf_size = buf1->get_bytes() / sizeof(hardware::buffers::PRNGBuffer::prng_state_t);

	auto prng1_state = createPrngState( buf_size);
	auto prng2_state = createPrngState( buf_size);

	buf1->dump(prng1_state);
	buf2->dump(prng2_state);

	// check that sufficiently small blocks don't match
	for(size_t i = 0; i < buf_size; ++i) {
		BOOST_CHECK_NE(prng1_state[i], prng2_state[i]);
	}

	delete[] prng1_state;
	delete[] prng2_state;
}

BOOST_AUTO_TEST_CASE(initialization)
{
	std::vector<std::string> parameterStrings {"--host_seed=13"};
	auto parameters = createParameters(parameterStrings).release();
	auto system = new hardware::System(*parameters);
	PRNG prng(*system);

	std::vector<std::string> parameterStrings2 {"--host_seed=14"};
	auto parameters2 = createParameters(parameterStrings2).release();
	auto system2 = new hardware::System(*parameters2);
	PRNG prng2(*system2);

	BOOST_CHECK_NE(prng.get_double(), prng2.get_double());

	logger.info() << "Now checking buffers...";
	for(size_t i = 0; i < prng.get_buffers().size(); ++i) {
		verifyBothBuffersAreEquallyLarge( prng.get_buffers().at(i), prng2.get_buffers().at(i) );
		verifyBuffersAreDifferent( prng.get_buffers().at(i), prng2.get_buffers().at(i) );
		logger.info() << "Checked buffer " << i;
	}
	logger.info() << "...done";
}

BOOST_AUTO_TEST_CASE(store_and_resume)
{
	std::vector<std::string> parameterStrings {"--host_seed=5"};
	auto parameters = createParameters(parameterStrings).release();
	auto system = new hardware::System(*parameters);

	PRNG prng(*system);
	prng.store("tmp.prngstate");

	double prng1_res = prng.get_double();

	std::vector<std::string> parameterStrings2 {"--initial_prng_state=tmp.prngstate"};
	auto parameters2 = createParameters(parameterStrings2).release();
	auto system2 = new hardware::System(*parameters2);
	PRNG prng2(*system2);

	BOOST_CHECK_EQUAL(prng1_res, prng2.get_double());

	logger.info() << "Now checking buffers...";
	for(size_t i = 0; i < prng.get_buffers().size(); ++i) {
		auto buf1 = prng.get_buffers().at(i);
		auto buf2 = prng2.get_buffers().at(i);

		char* prng1_state = new char[buf1->get_bytes()];
		buf1->dump(reinterpret_cast<hardware::buffers::PRNGBuffer::prng_state_t *>(prng1_state));

		char* prng2_state = new char[buf2->get_bytes()];
		buf2->dump(reinterpret_cast<hardware::buffers::PRNGBuffer::prng_state_t *>(prng2_state));

		BOOST_CHECK_EQUAL_COLLECTIONS(prng1_state, prng1_state + buf1->get_bytes(), prng2_state, prng2_state + buf2->get_bytes());
		logger.info() << "Checked buffer " << i;

		delete[] prng1_state;
		delete[] prng2_state;
	}
	logger.info() << "...done";
}
