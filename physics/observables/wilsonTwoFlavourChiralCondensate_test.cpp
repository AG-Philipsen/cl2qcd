/** @file
 * Unit test for the physics::observables::wilson::TwoFlavourChiralCondensate class
 *
 * Copyright 2014,Christopher Pinke
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

#include "wilsonTwoFlavourChiralCondensate.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::observables::wilson::TwoFlavourChiralCondensate
#include <boost/test/unit_test.hpp>

#include "../../host_functionality/logger.hpp"
#include <stdexcept>

BOOST_AUTO_TEST_SUITE( BUILD )

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
  const char * _params[] = {"foo"};
  meta::Inputparameters params(1, _params);

  BOOST_REQUIRE_NO_THROW(physics::observables::wilson::TwoFlavourChiralCondensate tester(&params) );
}

BOOST_AUTO_TEST_SUITE_END()

