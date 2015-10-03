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

#include "openClCode.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::OpenClCode
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(implementByMeansOfMetaInputparameters)
{
	const char * argv []  = {"foo"};
	meta::Inputparameters fullParameters{1, argv};
	hardware::OpenClCode_fromMetaInputparameters codeBuilder( fullParameters );

//	BOOST_REQUIRE_EQUAL( hardwareParameters.useGpu(), fullParameters.get_use_gpu() );

}
