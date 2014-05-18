/** @file
 * Unit test for the physics::gaugeObservables class
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

#include "gaugeObservables.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::gaugeObservables
#include <boost/test/unit_test.hpp>

#include "../../host_functionality/logger.hpp"
#include "../../meta/type_ops.hpp"
#include <stdexcept>

#include "../lattices/gaugefield.hpp"

BOOST_AUTO_TEST_SUITE( BUILD )

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
  const char * _params[] = {"foo"};
  meta::Inputparameters params(1, _params);

  BOOST_REQUIRE_NO_THROW(physics::gaugeObservables tester(&params) );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( PLAQUETTE  )

BOOST_AUTO_TEST_CASE( PLAQUETTE_1 )
{
  const char * _params[] = {"foo", "--startcondition=cold"};
  meta::Inputparameters params(2, _params);

  physics::gaugeObservables tester(&params);
  hardware::System system(params);
  physics::PRNG prng(system);
  physics::lattices::Gaugefield gaugefield(system, prng);

  tester.measurePlaquette(&gaugefield);
  BOOST_CHECK_CLOSE(tester.getPlaquette(), 1., 1e-8);
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_2 )
{
  double referenceValue = 0.57107711169452713;
  double testPrecision = 1e-8;

  const char * _params[] = {"foo", "--startcondition=continue", "--sourcefile=conf.00200", "--nt=4"};
  meta::Inputparameters params(4, _params);

  physics::gaugeObservables tester(&params);
  hardware::System system(params);
  physics::PRNG prng(system);
  physics::lattices::Gaugefield gaugefield(system, prng);

  tester.measurePlaquette(&gaugefield);
  BOOST_CHECK_CLOSE(tester.getPlaquette(), referenceValue, testPrecision);
}

BOOST_AUTO_TEST_SUITE_END()

