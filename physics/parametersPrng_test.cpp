/**
 * Copyright 2015 Christopher Pinke
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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::ParametersPrng
#include <boost/test/unit_test.hpp>

#include "parametersPrng.hpp"

BOOST_AUTO_TEST_CASE( init )
{
	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);

	physics::PrngParametersImplementation test(&params);

	BOOST_CHECK_EQUAL(test.getHostSeed(), params.get_host_seed());
	BOOST_CHECK_EQUAL(test.getInitialPrngStateFilename(), params.get_initial_prng_state());
	BOOST_CHECK_EQUAL(test.useSameRandomNumbers(), params.get_use_same_rnd_numbers() );
	BOOST_CHECK_EQUAL(test.getNamePostfix(), params.get_prng_postfix() );
	BOOST_CHECK_EQUAL(test.getNamePrefix(), params.get_prng_prefix() );
	BOOST_CHECK_EQUAL(test.getNumberOfDigitsInName(), params.get_config_number_digits() );
}



