/** @file
 * Definition of fermionmatrix operations.
 *
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

#pragma once

#include "../lattices/latticesInterfaces.hpp"

namespace physics {

    namespace fermionmatrix {

        class FermionmatrixParametersInterface {
            public:
                virtual ~FermionmatrixParametersInterface(){}
                virtual common::action getFermionicActionType() const = 0;
                virtual bool useMergedFermionicKernels() const = 0;
        };

        class FermionmatrixStaggeredParametersInterface {
            public:
                virtual ~FermionmatrixStaggeredParametersInterface() = 0;
        };
        //Pure virtual destructors must be implemented outside the class! (inline for multiple inclusion of header)
        inline FermionmatrixStaggeredParametersInterface::~FermionmatrixStaggeredParametersInterface(){}

    }

    //TODO: Think to a better place where to put these interfaces
    class FermionParametersInterface : public fermionmatrix::FermionmatrixParametersInterface, public lattices::SpinorfieldParametersInterface {
        public:
            virtual ~FermionParametersInterface(){}
    };

    class FermionEoParametersInterface : public fermionmatrix::FermionmatrixParametersInterface, public lattices::SpinorfieldEoParametersInterface {
        public:
            virtual ~FermionEoParametersInterface(){}
    };

    class FermionStaggeredEoParametersInterface : public fermionmatrix::FermionmatrixStaggeredParametersInterface, public lattices::StaggeredfieldEoParametersInterface {
        public:
            virtual ~FermionStaggeredEoParametersInterface(){}
    };

}


