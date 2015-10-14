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
    const physics::lattices::StaggaredfieldEoParametersImplementation staggaredfieldEoParametersImplementation{*params};
    const physics::lattices::RootedStaggaredfieldEoParametersImplementation rootedStaggaredfieldEoParametersImplementation{*params};
    const physics::fermionmatrix::FermionmatrixParametersImplementation fermionmatrixParametersImplementation{*params};
    const physics::observables::GaugeObservablesParametersImplementation gaugeObservablesParametersImplementation{*params};
    const physics::observables::WilsonTwoFlavourChiralCondensateParametersImplementation wilsonTwoFlavourChiralCondensateParametersImplementation{*params};
    const physics::observables::StaggeredChiralCondensateParametersImplementation staggeredChiralCondensateParametersImplementation{*params};
    const physics::observables::WilsonTwoFlavourCorrelatorsParametersImplementation wilsonTwoFlavourCorrelatorsParametersImplementation{*params};

    BOOST_CHECK( typeid(gaugefieldParametersImplementation) == typeid(test.getGaugefieldParametersInterface()) );
    BOOST_CHECK( typeid(gaugemomentaParametersImplementation) == typeid(test.getGaugemomentaParametersInterface()) );
    BOOST_CHECK( typeid(spinorfieldParametersImplementation) == typeid(test.getSpinorfieldParametersInterface()) );
    BOOST_CHECK( typeid(staggaredfieldEoParametersImplementation) == typeid(test.getStaggeredfieldEoParametersInterface()) );
    BOOST_CHECK( typeid(rootedStaggaredfieldEoParametersImplementation) == typeid(test.getRootedStaggeredfieldEoParametersInterface()) );
    BOOST_CHECK( typeid(fermionmatrixParametersImplementation) == typeid(test.getFermionmatrixParametersInterface()) );
    BOOST_CHECK( typeid(gaugeObservablesParametersImplementation) == typeid(test.getGaugeObservablesParametersInterface()) );
    BOOST_CHECK( typeid(wilsonTwoFlavourChiralCondensateParametersImplementation) == typeid(test.getWilsonTwoFlavourChiralCondensateParametersInterface()) );
    BOOST_CHECK( typeid(staggeredChiralCondensateParametersImplementation) == typeid(test.getStaggeredChiralCondensateParametersInterface()) );
    BOOST_CHECK( typeid(wilsonTwoFlavourCorrelatorsParametersImplementation) == typeid(test.getWilsonTwoFlavourCorrelatorsCondensateParametersInterface()) );
}
