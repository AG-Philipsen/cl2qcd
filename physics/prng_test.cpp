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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::PRNG
#include <boost/test/unit_test.hpp>

#include "../interfaceImplementations/physicsParameters.hpp"
#include "../interfaceImplementations/hardwareParameters.hpp"
#include "../interfaceImplementations/openClKernelParameters.hpp"

using namespace physics;

BOOST_AUTO_TEST_SUITE(build)

	BOOST_AUTO_TEST_CASE(brokenInputFile_brokenTag)
	{
		const char * _params[] = {"foo", "--initial_prng_state=prngstate_brokenTag"};
		meta::Inputparameters parameters(2, _params);
		physics::PrngParametersImplementation prngParameters(parameters);
	    hardware::HardwareParametersImplementation hP(&parameters);
	    hardware::code::OpenClKernelParametersImplementation kP(parameters);
	    hardware::System system(hP, kP);
		BOOST_CHECK_THROW( PRNG prng(system, &prngParameters) , std::invalid_argument );
	}

BOOST_AUTO_TEST_SUITE_END()

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
	const char * _params[] = {"foo", "--host_seed=13"};
	meta::Inputparameters parameters(2, _params);
	physics::PrngParametersImplementation prngParameters(parameters);
    hardware::HardwareParametersImplementation hP(&parameters);
    hardware::code::OpenClKernelParametersImplementation kP(parameters);
    hardware::System system(hP, kP);
	PRNG prng(system, &prngParameters);

	const char * _params2[] = {"foo", "--host_seed=14"};
	meta::Inputparameters parameters2(2, _params2);
	physics::PrngParametersImplementation prngParameters2(parameters2);
    hardware::HardwareParametersImplementation hP2(&parameters2);
    hardware::code::OpenClKernelParametersImplementation kP2(parameters2);
    hardware::System system2(hP2, kP2);
	PRNG prng2(system2, &prngParameters2);

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
	const char * _params[] = {"foo", "--host_seed=46"};
	meta::Inputparameters parameters(2, _params);
	physics::PrngParametersImplementation prngParameters(parameters);
    hardware::HardwareParametersImplementation hP(&parameters);
    hardware::code::OpenClKernelParametersImplementation kP(parameters);
    hardware::System system(hP, kP);
	PRNG prng(system, &prngParameters);
	prng.store("tmp.prngstate");

	double tmp = prng.get_double();

	const char * _params2[] = {"foo", "--initial_prng_state=tmp.prngstate"};
	meta::Inputparameters parameters2(2, _params2);
	physics::PrngParametersImplementation prngParameters2(parameters2);
    hardware::HardwareParametersImplementation hP2(&parameters2);
    hardware::code::OpenClKernelParametersImplementation kP2(parameters2);
    hardware::System system2(hP2, kP2);
	PRNG prng2(system2, &prngParameters2);

	double tmp2 = prng2.get_double();

	BOOST_CHECK_EQUAL(tmp, tmp2);

	BOOST_REQUIRE_EQUAL(prng == prng2, true);
}
