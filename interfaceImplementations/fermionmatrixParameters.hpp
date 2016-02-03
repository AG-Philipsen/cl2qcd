/** @file
 * fermionmatrixParameters.hpp
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

#pragma once

#include "../meta/inputparameters.hpp"
#include "../physics/fermionmatrix/fermionmatrixInterfaces.hpp"

namespace physics {
    namespace fermionmatrix {

        class FermionmatrixParametersImplementation : public FermionmatrixParametersInterface {
            public:
                FermionmatrixParametersImplementation() = delete;
                FermionmatrixParametersImplementation(const meta::Inputparameters& paramsIn)
                    : parameters(paramsIn)
                {
                }
                virtual ~FermionmatrixParametersImplementation()
                {
                }
                common::action getFermionicActionType() const override
                {
                    return parameters.get_fermact();
                }
                bool useMergedFermionicKernels() const override
                {
                    return parameters.get_use_merge_kernels_fermion();
                }
            private:
                const meta::Inputparameters& parameters;
        };

        class FermionmatrixStaggeredParametersImplementation : public FermionmatrixStaggeredParametersInterface {
            public:
                FermionmatrixStaggeredParametersImplementation() = delete;
                FermionmatrixStaggeredParametersImplementation(const meta::Inputparameters& parametersIn)
                : parameters(parametersIn)
                {
                }
                virtual ~FermionmatrixStaggeredParametersImplementation()
                {
                }
            private:
                const meta::Inputparameters& parameters;
        };

    }
}

