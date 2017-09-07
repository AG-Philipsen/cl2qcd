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
#define BOOST_TEST_MODULE physics::fermionmatrix::parametersInterface
#include <boost/test/unit_test.hpp>

#include "fermionmatrixParameters.hpp"

static std::unique_ptr<const meta::Inputparameters> createDefaultMetaInputparameters()
{
    const char * _params[] = {"foo"};
    return std::unique_ptr<meta::Inputparameters>(new meta::Inputparameters(1, _params) );
}

BOOST_AUTO_TEST_CASE(testFermionmatrixParameters)
{
    auto params = createDefaultMetaInputparameters();
    physics::fermionmatrix::FermionmatrixParametersImplementation test(*params);

    BOOST_CHECK_EQUAL(test.getFermionicActionType(), params->get_fermact());
    BOOST_CHECK_EQUAL(test.useMergedFermionicKernels(), params->get_use_merge_kernels_fermion());
}
