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

#include "../../meta/inputparameters.hpp"

namespace physics {
    namespace fermionmatrix {

        class FermionmatrixParametersInterface {
            public:
                virtual ~FermionmatrixParametersInterface()
                {
                }
                virtual common::action getFermionicActionType() const = 0;
                virtual bool useMergedFermionicKernels() const = 0;
        };

        class FermionmatrixStaggeredParametersInterface {
            public:
                virtual ~FermionmatrixStaggeredParametersInterface() = 0;
        };
        //Pure virtual destructors must be implemented outside the class! (inline for multiple inclusion of header)
        inline FermionmatrixStaggeredParametersInterface::~FermionmatrixStaggeredParametersInterface()
        {
        }

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

//TODO: Think to a better place where to implement these interfaces
#include "../lattices/latticesInterfaces.hpp"

namespace physics {

    class FermionParametersInterface : public fermionmatrix::FermionmatrixParametersInterface, public lattices::SpinorfieldParametersInterface {
        public:
            virtual ~FermionParametersInterface()
            {
            }
    };

    class FermionEoParametersInterface : public fermionmatrix::FermionmatrixParametersInterface, public lattices::SpinorfieldEoParametersInterface {
        public:
            virtual ~FermionEoParametersInterface()
            {
            }
    };

    class FermionStaggeredEoParametersInterface : public fermionmatrix::FermionmatrixStaggeredParametersInterface, public lattices::StaggeredfieldEoParametersInterface {
        public:
            virtual ~FermionStaggeredEoParametersInterface()
            {
            }
    };


    class FermionParametersImplementation final : public FermionParametersInterface,
                                                  private lattices::SpinorfieldParametersImplementation,
                                                  private fermionmatrix::FermionmatrixParametersImplementation {
        public:
            FermionParametersImplementation() = delete;
            FermionParametersImplementation(const meta::Inputparameters& parametersIn)
                    : lattices::SpinorfieldParametersImplementation(parametersIn), fermionmatrix::FermionmatrixParametersImplementation(parametersIn)
            {
            }
            unsigned getNt() const override
            {
                return lattices::SpinorfieldParametersImplementation::getNt();
            }
            unsigned getNs() const override
            {
                return lattices::SpinorfieldParametersImplementation::getNs();
            }
            unsigned getNumberOfElements() const override
            {
                return lattices::SpinorfieldParametersImplementation::getNumberOfElements();
            }
            common::action getFermionicActionType() const override
            {
                return fermionmatrix::FermionmatrixParametersImplementation::getFermionicActionType();
            }
            bool useMergedFermionicKernels() const override
            {
                return fermionmatrix::FermionmatrixParametersImplementation::useMergedFermionicKernels();
            }
    };

    class FermionEoParametersImplementation final : public FermionEoParametersInterface,
                                                    private lattices::SpinorfieldEoParametersImplementation,
                                                    private fermionmatrix::FermionmatrixParametersImplementation {
        public:
            FermionEoParametersImplementation() = delete;
            FermionEoParametersImplementation(const meta::Inputparameters& parametersIn)
                    : lattices::SpinorfieldEoParametersImplementation(parametersIn), fermionmatrix::FermionmatrixParametersImplementation(parametersIn)
            {
            }
            common::action getFermionicActionType() const override
            {
                return fermionmatrix::FermionmatrixParametersImplementation::getFermionicActionType();
            }
            bool useMergedFermionicKernels() const override
            {
                return fermionmatrix::FermionmatrixParametersImplementation::useMergedFermionicKernels();
            }
    };

    class FermionStaggeredEoParametersImplementation: public FermionStaggeredEoParametersInterface,
                                                      private lattices::StaggeredfieldEoParametersImplementation,
                                                      private fermionmatrix::FermionmatrixStaggeredParametersImplementation {
        public:
            FermionStaggeredEoParametersImplementation() = delete;
            FermionStaggeredEoParametersImplementation(const meta::Inputparameters& parametersIn)
                    : lattices::StaggeredfieldEoParametersImplementation(parametersIn), fermionmatrix::FermionmatrixStaggeredParametersImplementation(parametersIn)
            {
            }
            unsigned getNumberOfElements() const override
            {
                return lattices::StaggeredfieldEoParametersImplementation::getNumberOfElements();
            }
    };

}

