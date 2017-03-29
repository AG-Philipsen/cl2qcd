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
#define BOOST_TEST_MODULE physics::lattice::parametersInterface
#include <boost/test/unit_test.hpp>

#include "latticesParameters.hpp"

static std::unique_ptr<const meta::Inputparameters> createDefaultMetaInputparameters()
{
	const char * _params[] = {"foo"};
	return std::unique_ptr<meta::Inputparameters>(new meta::Inputparameters(1, _params) );
}

BOOST_AUTO_TEST_CASE(testLatticeObjectParameters)
{
	auto params = createDefaultMetaInputparameters();
	physics::lattices::GaugefieldParametersImplementation test(&(*params));

	BOOST_CHECK_EQUAL(test.getNs(), params->get_nspace());
	BOOST_CHECK_EQUAL(test.getNt(), params->get_ntime());
	BOOST_CHECK_EQUAL(test.getPrecision(), params->get_precision());
	BOOST_CHECK_EQUAL(test.ignoreChecksumErrorsInIO(), params->get_ignore_checksum_errors());
	BOOST_CHECK_EQUAL(test.getNumberOfElements(), std::pow(params->get_nspace(), 3.)*params->get_ntime()*NDIM);
	BOOST_CHECK_EQUAL(test.getKappa(), params->get_kappa());
	BOOST_CHECK_EQUAL(test.getMu(), params->get_mu());
	BOOST_CHECK_EQUAL(test.getBeta(), params->get_beta());
	BOOST_CHECK_EQUAL(test.getStartcondition(), params->get_startcondition());
	BOOST_CHECK_EQUAL(test.getNamePostfix(), params->get_config_postfix());
	BOOST_CHECK_EQUAL(test.getNamePrefix(), params->get_config_prefix());
	BOOST_CHECK_EQUAL(test.getNumberOfDigitsInName(), params->get_config_number_digits());
	BOOST_CHECK_EQUAL(test.getSmearingSteps(), params->get_rho_iter());
	BOOST_CHECK_EQUAL(test.getSourcefileName(), params->get_sourcefile());
}

BOOST_AUTO_TEST_CASE(testGaugemomentaParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::lattices::GaugemomentaParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getNs(), params->get_nspace());
    BOOST_CHECK_EQUAL(test.getNt(), params->get_ntime());
    BOOST_CHECK_EQUAL(test.getNumberOfElements(), std::pow(params->get_nspace(), 3.)*params->get_ntime()*NDIM);
}

BOOST_AUTO_TEST_CASE(testSpinorfieldParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::lattices::SpinorfieldParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getNs(), params->get_nspace());
    BOOST_CHECK_EQUAL(test.getNt(), params->get_ntime());
    BOOST_CHECK_EQUAL(test.getNumberOfElements(), std::pow(params->get_nspace(), 3.)*params->get_ntime());
}

BOOST_AUTO_TEST_CASE(testStaggeredfieldEoParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::lattices::StaggeredfieldEoParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getNumberOfElements(), std::pow(params->get_nspace(), 3.)*params->get_ntime());
}

BOOST_AUTO_TEST_CASE(testRootedStaggeredfieldEoParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::lattices::RootedStaggeredfieldEoParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getNumberOfElements(), std::pow(params->get_nspace(), 3.)*params->get_ntime());
    BOOST_CHECK_EQUAL(test.getMetropolisRationalApproximationOrder(), params->get_metro_approx_ord());
    BOOST_CHECK_EQUAL(test.getMolecularDynamicsRationalApproximationOrder(), params->get_md_approx_ord());
    BOOST_CHECK_EQUAL(test.getNumberOfPseudofermions(), params->get_num_pseudofermions());
}
