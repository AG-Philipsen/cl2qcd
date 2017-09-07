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
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/.
 */

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::interfaceHandler
#include <boost/test/unit_test.hpp>

#include "interfacesHandler.hpp"
//#include "../physics/interfacesHandler.hpp"
#include "../physics/lattices/gaugefield.hpp"
#include "../physics/lattices/gaugemomenta.hpp"
#include "../physics/lattices/spinorfield.hpp"
#include "../physics/lattices/staggeredfield_eo.hpp"
#include "../physics/lattices/rooted_staggeredfield_eo.hpp"
#include "../physics/fermionmatrix/fermionmatrix.hpp"

static std::unique_ptr<const meta::Inputparameters> createDefaultMetaInputparameters()
{
    const char * _params[] = {"foo"};
    return std::unique_ptr<meta::Inputparameters>(new meta::Inputparameters(1, _params) );
}

BOOST_AUTO_TEST_CASE(testInterfaceHandler)
{
    auto params = createDefaultMetaInputparameters();
    physics::InterfacesHandlerImplementation test{*params};
    const physics::lattices::GaugefieldParametersImplementation gaugefieldParametersImplementation{params.get()};
    const physics::lattices::GaugemomentaParametersImplementation gaugemomentaParametersImplementation{*params};
    const physics::lattices::SpinorfieldParametersImplementation spinorfieldParametersImplementation{*params};
    const physics::lattices::SpinorfieldEoParametersImplementation spinorfieldEoParametersImplementation{*params};
    const physics::lattices::StaggeredfieldEoParametersImplementation staggaredfieldEoParametersImplementation{*params};
    const physics::lattices::RootedStaggeredfieldEoParametersImplementation rootedStaggaredfieldEoParametersImplementation{*params};
    const physics::fermionmatrix::FermionmatrixParametersImplementation fermionmatrixParametersImplementation{*params};
    const physics::FermionParametersImplementation fermionParametersImplementation{*params};
    const physics::FermionEoParametersImplementation fermionEoParametersImplementation{*params};
    const physics::observables::GaugeObservablesParametersImplementation gaugeObservablesParametersImplementation{*params};
    const physics::observables::WilsonTwoFlavourChiralCondensateParametersImplementation wilsonTwoFlavourChiralCondensateParametersImplementation{*params};
    const physics::observables::StaggeredChiralCondensateParametersImplementation staggeredChiralCondensateParametersImplementation{*params};
    const physics::observables::WilsonTwoFlavourCorrelatorsParametersImplementation wilsonTwoFlavourCorrelatorsParametersImplementation{*params};
    const physics::observables::StaggeredTwoFlavourCorrelatorsParametersImplementation staggeredTwoFlavourCorrelatorsParametersImplementation{*params};
    const physics::algorithms::SolversParametersImplementation solversParametersImplementation{*params};
    const physics::algorithms::MinMaxEigenvalueParametersImplementation minMaxEigenvalueParametersImplementation{*params};
    const physics::algorithms::ForcesParametersImplementation forcesParametersImplementation{*params};
    const physics::algorithms::InversionParametersImplementation inversionParametersImplementation{*params};
    const physics::algorithms::IntegratorParametersImplementation integratorParametersImplementation{*params};
    const physics::algorithms::MolecularDynamicsImplementation molecularDynamicsImplementation{*params};
    const physics::algorithms::MetropolisParametersImplementation metropolisParametersImplementation{*params};
    const physics::algorithms::HmcParametersImplementation hmcParametersImplementation{*params};
    const physics::algorithms::RhmcParametersImplementation rhmcParametersImplementation{*params};
    const physics::SourcesParametersImplementation sourcesParametersImplementation{*params};
    const physics::WilsonAdditionalParameters wilsonAdditionalParameters{*params, true};
    const physics::StaggeredAdditionalParameters staggeredAdditionalParameters{*params};

    BOOST_CHECK( typeid(gaugefieldParametersImplementation) == typeid(test.getInterface<physics::lattices::Gaugefield>()) );
    BOOST_CHECK( typeid(gaugemomentaParametersImplementation) == typeid(test.getInterface<physics::lattices::Gaugemomenta>()) );
    BOOST_CHECK( typeid(spinorfieldParametersImplementation) == typeid(test.getInterface<physics::lattices::Spinorfield>()) );
    BOOST_CHECK( typeid(spinorfieldEoParametersImplementation) == typeid(test.getInterface<physics::lattices::Spinorfield_eo>()) );
    BOOST_CHECK( typeid(staggaredfieldEoParametersImplementation) == typeid(test.getInterface<physics::lattices::Staggeredfield_eo>()) );
    BOOST_CHECK( typeid(rootedStaggaredfieldEoParametersImplementation) == typeid(test.getInterface<physics::lattices::Rooted_Staggeredfield_eo>()) );
    BOOST_CHECK( typeid(fermionmatrixParametersImplementation) == typeid(test.getInterface<physics::fermionmatrix::M>()) );
    BOOST_CHECK( typeid(fermionmatrixParametersImplementation) == typeid(test.getInterface<physics::fermionmatrix::Qplus>()) );
    BOOST_CHECK( typeid(fermionmatrixParametersImplementation) == typeid(test.getInterface<physics::fermionmatrix::Qminus>()) );
    BOOST_CHECK( typeid(fermionParametersImplementation) == typeid(test.getInterface<physics::fermionmatrix::QplusQminus>()) );
    BOOST_CHECK( typeid(fermionEoParametersImplementation) == typeid(test.getInterface<physics::fermionmatrix::Aee>()) );
    BOOST_CHECK( typeid(fermionEoParametersImplementation) == typeid(test.getInterface<physics::fermionmatrix::Aee_minus>()) );
    BOOST_CHECK( typeid(fermionEoParametersImplementation) == typeid(test.getInterface<physics::fermionmatrix::Qplus_eo>()) );
    BOOST_CHECK( typeid(fermionEoParametersImplementation) == typeid(test.getInterface<physics::fermionmatrix::Qminus_eo>()) );
    BOOST_CHECK( typeid(fermionEoParametersImplementation) == typeid(test.getInterface<physics::fermionmatrix::QplusQminus_eo>()) );
    BOOST_CHECK( typeid(gaugeObservablesParametersImplementation) == typeid(test.getGaugeObservablesParametersInterface()) );
    BOOST_CHECK( typeid(wilsonTwoFlavourChiralCondensateParametersImplementation) == typeid(test.getWilsonTwoFlavourChiralCondensateParametersInterface()) );
    BOOST_CHECK( typeid(staggeredChiralCondensateParametersImplementation) == typeid(test.getStaggeredChiralCondensateParametersInterface()) );
    BOOST_CHECK( typeid(wilsonTwoFlavourCorrelatorsParametersImplementation) == typeid(test.getWilsonTwoFlavourCorrelatorsParametersInterface()) );
    BOOST_CHECK( typeid(staggeredTwoFlavourCorrelatorsParametersImplementation) == typeid(test.getStaggeredTwoFlavourCorrelatorsParametersInterface()) );
    BOOST_CHECK( typeid(solversParametersImplementation) == typeid(test.getSolversParametersInterface()) );
    BOOST_CHECK( typeid(minMaxEigenvalueParametersImplementation) == typeid(test.getMinMaxEigenvalueParametersInterface()) );
    BOOST_CHECK( typeid(forcesParametersImplementation) == typeid(test.getForcesParametersInterface()) );
    BOOST_CHECK( typeid(inversionParametersImplementation) == typeid(test.getInversionParemetersInterface()) );
    BOOST_CHECK( typeid(integratorParametersImplementation) == typeid(test.getIntegratorParametersInterface()) );
    BOOST_CHECK( typeid(molecularDynamicsImplementation) == typeid(test.getMolecularDynamicsInterface()) );
    BOOST_CHECK( typeid(metropolisParametersImplementation) == typeid(test.getMetropolisParametersInterface()) );
    BOOST_CHECK( typeid(hmcParametersImplementation) == typeid(test.getHmcParametersInterface()) );
    BOOST_CHECK( typeid(rhmcParametersImplementation) == typeid(test.getRhmcParametersInterface()) );
    BOOST_CHECK( typeid(sourcesParametersImplementation) == typeid(test.getSourcesParametersInterface()) );
    BOOST_CHECK( typeid(wilsonAdditionalParameters) == typeid(test.getAdditionalParameters<physics::lattices::Spinorfield>(true)) );
    BOOST_CHECK( typeid(wilsonAdditionalParameters) == typeid(test.getAdditionalParameters<physics::lattices::Spinorfield>(false)) );
    BOOST_CHECK( typeid(wilsonAdditionalParameters) == typeid(test.getAdditionalParameters<physics::lattices::Spinorfield_eo>(true)) );
    BOOST_CHECK( typeid(wilsonAdditionalParameters) == typeid(test.getAdditionalParameters<physics::lattices::Spinorfield_eo>(false)) );
    BOOST_CHECK( typeid(staggeredAdditionalParameters) == typeid(test.getAdditionalParameters<physics::lattices::Staggeredfield_eo>()) );
    BOOST_CHECK( typeid(staggeredAdditionalParameters) == typeid(test.getAdditionalParameters<physics::lattices::Rooted_Staggeredfield_eo>()) );
}
