/** @file
 * algorithms interfaces tests
 *
 * Copyright 2016 Alessandro Sciarra, Christopher Czaban
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
#define BOOST_TEST_MODULE physics::fermionmatrix::parametersInterface
#include <boost/test/unit_test.hpp>

#include "algorithmsParameters.hpp"

static std::unique_ptr<const meta::Inputparameters> createDefaultMetaInputparameters()
{
    const char * _params[] = {"foo"};
    return std::unique_ptr<meta::Inputparameters>(new meta::Inputparameters(1, _params) );
}

BOOST_AUTO_TEST_CASE(testSolversParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::algorithms::SolversParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getCgIterationBlockSize(), params->get_cg_iteration_block_size());
    BOOST_CHECK_EQUAL(test.getCgMax(), params->get_cgmax());
    BOOST_CHECK_EQUAL(test.getCgMinimumIterationCount(), params->get_cg_minimum_iteration_count());
    BOOST_CHECK_EQUAL(test.getCgUseAsyncCopy(), params->get_cg_use_async_copy());
    BOOST_CHECK_EQUAL(test.getIterRefresh(), params->get_iter_refresh());
    BOOST_CHECK_EQUAL(test.getSolver(), params->get_solver());
    BOOST_CHECK_EQUAL(test.getUseMergeKernelsSpinor(), params->get_use_merge_kernels_spinor());
}

BOOST_AUTO_TEST_CASE(testForcesParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::algorithms::ForcesParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getFermact(), params->get_fermact());
    BOOST_CHECK_EQUAL(test.getForcePreconditioning(), params->get_force_prec());
    BOOST_CHECK_EQUAL(test.getRhoIterations(), params->get_rho_iter());
    BOOST_CHECK_EQUAL(test.getSolver(), params->get_solver());
    BOOST_CHECK_EQUAL(test.getUseGaugeOnly(), params->get_use_gauge_only());
    BOOST_CHECK_EQUAL(test.getUseRectangles(), meta::get_use_rectangles(*params));
    BOOST_CHECK_EQUAL(test.getUseSmearing(), params->get_use_smearing());
}

BOOST_AUTO_TEST_CASE(testMinMaxEigenvalueParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::algorithms::MinMaxEigenvalueParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getFindMinMaxIterationBlockSize(), params->get_findminmax_iteration_block_size());
    BOOST_CHECK_EQUAL(test.getFindMinMaxMaxValue(), params->get_findminmax_max());
}

BOOST_AUTO_TEST_CASE(testInversionParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::algorithms::InversionParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getFermact(), params->get_fermact());
    BOOST_CHECK_EQUAL(test.getSolver(), params->get_solver());
    BOOST_CHECK_EQUAL(test.getSolverPrec(), params->get_solver_prec());
    BOOST_CHECK_EQUAL(test.getUseEo(), params->get_use_eo());
    BOOST_CHECK_EQUAL(test.getUseSmearing(), params->get_use_smearing());
}

BOOST_AUTO_TEST_CASE(testIntegratorParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::algorithms::IntegratorParametersImplementation test(*params);

    for(int i=0; i<3; i++){
        BOOST_CHECK_EQUAL(test.getIntegrationSteps(i), params->get_integrationsteps(i));
        BOOST_CHECK_EQUAL(test.getIntegrator(i), params->get_integrator(i));
        BOOST_CHECK_EQUAL(test.getLambda(i), params->get_lambda(i));
    }
    BOOST_CHECK_EQUAL(test.getNumTimescales(), params->get_num_timescales());
    BOOST_CHECK_EQUAL(test.getTau(), params->get_tau());
    BOOST_CHECK_EQUAL(test.getUseMp(), params->get_use_mp());
}

BOOST_AUTO_TEST_CASE(testMolecularDynamicsParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::algorithms::MolecularDynamicsImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getSolverPrec(), params->get_solver_prec());
}

BOOST_AUTO_TEST_CASE(testMetropolisParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::algorithms::MetropolisParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getC0(), meta::get_c0(*params));
    BOOST_CHECK_EQUAL(test.getC1(), meta::get_c1(*params));
    BOOST_CHECK_EQUAL(test.getFermact(), params->get_fermact());
    BOOST_CHECK_EQUAL(test.getRectanglesNormalization(), meta::get_rect_norm(*params));
    BOOST_CHECK_EQUAL(test.getSolver(), params->get_solver());
    BOOST_CHECK_EQUAL(test.getSolverPrec(), params->get_solver_prec());
    BOOST_CHECK_EQUAL(test.getUseGaugeOnly(), params->get_use_gauge_only());
    BOOST_CHECK_EQUAL(test.getUseMp(), params->get_use_mp());
    BOOST_CHECK_EQUAL(test.getUseRectangles(), meta::get_use_rectangles(*params));
}

BOOST_AUTO_TEST_CASE(testHmcParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::algorithms::HmcParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getBeta(), params->get_beta());
    BOOST_CHECK_EQUAL(test.getUseEo(), params->get_use_eo());
    BOOST_CHECK_EQUAL(test.getUseGaugeOnly(), params->get_use_gauge_only());
    BOOST_CHECK_EQUAL(test.getUseMp(), params->get_use_mp());
}

BOOST_AUTO_TEST_CASE(testRhmcParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::algorithms::RhmcParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getBeta(), params->get_beta());
    BOOST_CHECK_EQUAL(test.getConservative(), params->get_conservative());
    BOOST_CHECK_EQUAL(test.getFindMinMaxPrec(), params->get_findminmax_prec());
    BOOST_CHECK_EQUAL(test.getMass(), params->get_mass());
    BOOST_CHECK_EQUAL(test.getUseGaugeOnly(), params->get_use_gauge_only());
    BOOST_CHECK_EQUAL(test.getUseMp(), params->get_use_mp());
    BOOST_CHECK_EQUAL(test.getUseEo(), params->get_use_eo());
}









