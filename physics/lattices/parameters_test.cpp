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
#define BOOST_TEST_MODULE physics::lattice::Parameters
#include <boost/test/unit_test.hpp>

#include "parameters.hpp"

std::unique_ptr<const meta::Inputparameters> createDefaultMetaInputparameters()
{
	const char * _params[] = {"foo"};
	return std::unique_ptr<meta::Inputparameters>(new meta::Inputparameters(1, _params) );
}

BOOST_AUTO_TEST_CASE(testLatticeObjectParameters)
{
	auto params = createDefaultMetaInputparameters();
	LatticeObjectParametersImplementation test(&(*params));

	BOOST_CHECK_EQUAL(test.getNs(), 4);
	BOOST_CHECK_EQUAL(test.getNt(), 8);
	BOOST_CHECK_EQUAL(test.getPrecision(), 64);
	BOOST_CHECK_EQUAL(test.ignoreChecksumErrorsInIO(), false);
	BOOST_CHECK_EQUAL(test.getNumberOfElements(), 4*4*4*8*4);
	BOOST_CHECK_EQUAL(test.getKappa(), 0.125);
	BOOST_CHECK_EQUAL(test.getMu(), 0.006);
	BOOST_CHECK_EQUAL(test.getBeta(), 4.);
	BOOST_CHECK_EQUAL(test.getStartcondition(), common::startcondition::cold_start);
	BOOST_CHECK_EQUAL(test.getNamePostfix(), "");
	BOOST_CHECK_EQUAL(test.getNamePrefix(), "conf.");
	BOOST_CHECK_EQUAL(test.getNumberOfDigitsInName(), 5);
	BOOST_CHECK_EQUAL(test.getSmearingSteps(), 0);
	BOOST_CHECK_EQUAL(test.getSourcefileName(), "conf.00000");
}

BOOST_AUTO_TEST_CASE(testSpinorfieldParameters)
{
    auto params = createDefaultMetaInputparameters();
    SpinorfieldParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getNs(), 4);
    BOOST_CHECK_EQUAL(test.getNt(), 8);
}

BOOST_AUTO_TEST_CASE(testGaugemomentaParameters)
{
    auto params = createDefaultMetaInputparameters();
    GaugemomentaParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getNs(), 4);
    BOOST_CHECK_EQUAL(test.getNt(), 8);
}

BOOST_AUTO_TEST_CASE(testStaggeredfieldEoParameters)
{
    auto params = createDefaultMetaInputparameters();
    StaggaredfieldEoParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getNs(), 4);
    BOOST_CHECK_EQUAL(test.getNt(), 8);
}

BOOST_AUTO_TEST_CASE(testRootedStaggeredfieldEoParameters)
{
    auto params = createDefaultMetaInputparameters();
    RootedStaggaredfieldEoParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getNs(), 4);
    BOOST_CHECK_EQUAL(test.getNt(), 8);
    BOOST_CHECK_EQUAL(test.getMetropolisRationalApproximationOrder(), 15);
    BOOST_CHECK_EQUAL(test.getMolecularDynamicsRationalApproximationOrder(), 8);
}
