/** @file
 * additionalParameters_test.cpp
 *
 * Copyright 2016 Alessandro Sciarra
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
#define BOOST_TEST_MODULE physics::additionalParameters
#include <boost/test/unit_test.hpp>

#include "physicsParameters.hpp"

static std::unique_ptr<const meta::Inputparameters> createDefaultMetaInputparameters()
{
    const char * _params[] = {"foo"};
    return std::unique_ptr<meta::Inputparameters>(new meta::Inputparameters(1, _params) );
}

BOOST_AUTO_TEST_SUITE(testWilsonAdditionalParameters)

    BOOST_AUTO_TEST_CASE(withoutMassPreconditioning)
    {
        auto params = createDefaultMetaInputparameters();
        physics::WilsonAdditionalParameters test(*params, false);

        BOOST_CHECK_EQUAL(test.getKappa(), params->get_kappa());
        BOOST_CHECK_EQUAL(test.getMubar(), meta::get_mubar(*params));
        BOOST_REQUIRE_THROW(test.getMass(), Print_Error_Message);
        BOOST_REQUIRE_THROW(test.getConservative(), Print_Error_Message);
    }

    BOOST_AUTO_TEST_CASE(withMassPreconditioning)
    {
        auto params = createDefaultMetaInputparameters();
        physics::WilsonAdditionalParameters test(*params, true);

        BOOST_CHECK_EQUAL(test.getKappa(), params->get_kappa_mp());
        BOOST_CHECK_EQUAL(test.getMubar(), meta::get_mubar_mp(*params));
        BOOST_REQUIRE_THROW(test.getMass(), Print_Error_Message);
        BOOST_REQUIRE_THROW(test.getConservative(), Print_Error_Message);
    }

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(testStaggeredAdditionalParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::StaggeredAdditionalParameters test(*params);

    BOOST_REQUIRE_THROW(test.getKappa(), Print_Error_Message);
    BOOST_REQUIRE_THROW(test.getMubar(), Print_Error_Message);
    BOOST_CHECK_EQUAL(test.getMass(), params->get_mass());
    BOOST_CHECK_EQUAL(test.getConservative(), params->get_conservative());
}

BOOST_AUTO_TEST_CASE(testFermionParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::FermionParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getNs(), params->get_nspace());
    BOOST_CHECK_EQUAL(test.getNt(), params->get_ntime());
    BOOST_CHECK_EQUAL(test.getNumberOfElements(), std::pow(params->get_nspace(), 3.)*params->get_ntime());
    BOOST_CHECK_EQUAL(test.getFermionicActionType(), params->get_fermact());
    BOOST_CHECK_EQUAL(test.useMergedFermionicKernels(), params->get_use_merge_kernels_fermion());
}

BOOST_AUTO_TEST_CASE(testFermionEoParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::FermionEoParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getFermionicActionType(), params->get_fermact());
    BOOST_CHECK_EQUAL(test.useMergedFermionicKernels(), params->get_use_merge_kernels_fermion());
}

BOOST_AUTO_TEST_CASE(testSourcesParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::SourcesParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.placeSourcesOnHost(), params->get_place_sources_on_host());
    BOOST_CHECK_EQUAL(test.getSourceType(), params->get_sourcetype());
    BOOST_CHECK_EQUAL(test.getSourceT(), params->get_source_t());
    BOOST_CHECK_EQUAL(test.getSourceZ(), params->get_source_z());

}

BOOST_AUTO_TEST_CASE(testPrngParameters)
{
    const char * _params[] = {"foo"};
    meta::Inputparameters params(1, _params);

    physics::PrngParametersImplementation test(params);

    BOOST_CHECK_EQUAL(test.getHostSeed(), params.get_host_seed());
    BOOST_CHECK_EQUAL(test.getInitialPrngStateFilename(), params.get_initial_prng_state());
    BOOST_CHECK_EQUAL(test.useSameRandomNumbers(), params.get_use_same_rnd_numbers() );
    BOOST_CHECK_EQUAL(test.getNamePostfix(), params.get_prng_postfix() );
    BOOST_CHECK_EQUAL(test.getNamePrefix(), params.get_prng_prefix() );
    BOOST_CHECK_EQUAL(test.getNumberOfDigitsInName(), params.get_config_number_digits() );
}
