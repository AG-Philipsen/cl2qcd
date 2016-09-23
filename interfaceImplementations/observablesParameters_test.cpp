/**
 * Copyright 2015 Alessandro Sciarra, Christopher Czaban
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
#define BOOST_TEST_MODULE physics::observables::parametersInterface
#include <boost/test/unit_test.hpp>

#include "observablesParameters.hpp"

static std::unique_ptr<const meta::Inputparameters> createDefaultMetaInputparameters()
{
    const char * _params[] = {"foo"};
    return std::unique_ptr<meta::Inputparameters>(new meta::Inputparameters(1, _params) );
}

BOOST_AUTO_TEST_CASE(testGaugeObservablesParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::observables::GaugeObservablesParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.measureRectangles(), params->get_measure_rectangles());
    BOOST_CHECK_EQUAL(test.measureTransportCoefficientKappa(), params->get_measure_transportcoefficient_kappa());
    BOOST_CHECK_EQUAL(test.printToScreen(), params->get_print_to_screen());
    BOOST_CHECK_EQUAL(boost::lexical_cast<std::string>(test.getBeta()), boost::lexical_cast<std::string>(params->get_beta()));
    BOOST_CHECK_EQUAL(test.getTransportCoefficientKappaFilename(), params->get_transportcoefficientKappaFilename());
    BOOST_CHECK_EQUAL(test.getRectanglesFilename(), params->get_rectanglesFilename());
    BOOST_CHECK_EQUAL(test.getGaugeObservablesFilename("conf.00000"), meta::get_gauge_obs_file_name(*params, "conf.00000"));
    BOOST_CHECK_EQUAL(test.getTemporalPlaquetteNormalization(), meta::get_tplaq_norm(*params));
    BOOST_CHECK_EQUAL(test.getSpatialPlaquetteNormalization(), meta::get_splaq_norm(*params));
    BOOST_CHECK_EQUAL(test.getPlaquetteNormalization(), meta::get_plaq_norm(*params));
    BOOST_CHECK_EQUAL(test.getSpatialVolume(), meta::get_volspace(*params));
    BOOST_CHECK_EQUAL(test.getPolyakovLoopNormalization(), meta::get_poly_norm(*params));
}

BOOST_AUTO_TEST_CASE(testStaggeredChiralCondensateObservablesParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::observables::StaggeredChiralCondensateParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getNumberOfSources(), params->get_num_sources());
    BOOST_CHECK_EQUAL(boost::lexical_cast<std::string>(test.getMass()), boost::lexical_cast<std::string>(params->get_mass()));
    BOOST_CHECK_EQUAL(test.measurePbp(), params->get_measure_pbp());
    BOOST_CHECK_EQUAL(test.getSolverPrecision(), params->get_solver_prec());
    BOOST_CHECK_EQUAL(test.getNumberOfTastes(), params->get_num_tastes());
    BOOST_CHECK_EQUAL(test.getPbpNumberOfMeasurements(), params->get_pbp_measurements());
    BOOST_CHECK_EQUAL(test.get4dVolume(), meta::get_vol4d(*params));
    BOOST_CHECK_EQUAL(test.getPbpFilename("conf.00000"), meta::get_ferm_obs_pbp_file_name(*params, "conf.00000"));
}

BOOST_AUTO_TEST_CASE(testWilsonChiralCondensateObservablesParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::observables::WilsonTwoFlavourChiralCondensateParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getFermionicActionType(), params->get_fermact());
    BOOST_CHECK_EQUAL(test.getPbpVersion(), params->get_pbp_version());
    BOOST_CHECK_EQUAL(test.measurePbp(), params->get_measure_pbp());
    BOOST_CHECK_EQUAL(test.getNumberOfSources(), params->get_num_sources());
    BOOST_CHECK_EQUAL(boost::lexical_cast<std::string>(test.getKappa()), boost::lexical_cast<std::string>(params->get_kappa()));
    BOOST_CHECK_EQUAL(test.get4dVolume(), meta::get_vol4d(*params));
    BOOST_CHECK_EQUAL(test.useEvenOdd(), params->get_use_eo());
    BOOST_CHECK_EQUAL(test.getPbpFilename("conf.00000"), meta::get_ferm_obs_pbp_file_name(*params, "conf.00000"));
    BOOST_CHECK_EQUAL(boost::lexical_cast<std::string>(test.getMubar()), boost::lexical_cast<std::string>(meta::get_mubar(*params)));
}

BOOST_AUTO_TEST_CASE(testWilsonCorrelatorsObservablesParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::observables::WilsonTwoFlavourCorrelatorsParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.printToScreen(), params->get_print_to_screen());
    BOOST_CHECK_EQUAL(test.getCorrelatorDirection(), params->get_corr_dir());
    BOOST_CHECK_EQUAL(test.getSourceType(), params->get_sourcetype());
    BOOST_CHECK_EQUAL(test.getNs(), params->get_nspace());
    BOOST_CHECK_EQUAL(test.getNt(), params->get_ntime());
    BOOST_CHECK_EQUAL(test.getCorrelatorFilename("conf.00000"), meta::get_ferm_obs_corr_file_name(*params, "conf.00000"));
    BOOST_CHECK_EQUAL(test.placeSourcesOnHost(), params->get_place_sources_on_host());
    BOOST_CHECK_EQUAL(test.getNumberOfSources(), params->get_num_sources());
}

BOOST_AUTO_TEST_CASE(testStaggeredCorrelatorsObservablesParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::observables::StaggeredTwoFlavourCorrelatorsParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.printToScreen(), params->get_print_to_screen());
//    BOOST_CHECK_EQUAL(test.getCorrelatorDirection(), params->get_corr_dir());
    BOOST_CHECK_EQUAL(test.getSourceType(), params->get_sourcetype());
//    BOOST_CHECK_EQUAL(test.getNs(), params->get_nspace());
    BOOST_CHECK_EQUAL(test.getNt(), params->get_ntime());
    BOOST_CHECK_EQUAL(test.getCorrelatorFilename("conf.00000"), meta::get_ferm_obs_corr_file_name(*params, "conf.00000"));
//    BOOST_CHECK_EQUAL(test.placeSourcesOnHost(), params->get_place_sources_on_host());
    BOOST_CHECK_EQUAL(test.getNumberOfSources(), params->get_num_sources());
    BOOST_CHECK_EQUAL(test.getSolverPrecision(), params->get_solver_prec());
}




