/** @file
 * Declaration of the physics::lattices::Gaugefield class
 *
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

#pragma once

#include "../../common_header_files/types.h"
#include <iostream>

namespace physics {
    namespace lattices {

        class LatticeObjectParametersInterface {
            public:
                virtual ~LatticeObjectParametersInterface(){}
                virtual unsigned getNs() const = 0;
                virtual unsigned getNt() const = 0;
                virtual unsigned getPrecision() const = 0;
                virtual bool ignoreChecksumErrorsInIO() const = 0;
                virtual unsigned getNumberOfElements() const = 0;
                virtual double getKappa() const = 0;
                virtual double getMu() const = 0;
                virtual double getBeta() const = 0;
                virtual common::startcondition getStartcondition() const = 0;
                virtual std::string getNamePrefix() const = 0;
                virtual std::string getNamePostfix() const = 0;
                virtual unsigned getNumberOfDigitsInName() const = 0;
                virtual unsigned getSmearingSteps() const = 0;
                virtual std::string getSourcefileName() const = 0;
        };

        class GaugemomentaParametersInterface {
            public:
                virtual ~GaugemomentaParametersInterface(){}
                virtual unsigned getNs() const = 0;
                virtual unsigned getNt() const = 0;
                virtual unsigned getNumberOfElements() const = 0;
        };

        class SpinorfieldParametersInterface {
            public:
                virtual ~SpinorfieldParametersInterface(){}
                virtual unsigned getNs() const = 0;
                virtual unsigned getNt() const = 0;
                virtual unsigned getNumberOfElements() const = 0;
        };

        class StaggaredfieldEoParametersInterface {
            public:
                virtual ~StaggaredfieldEoParametersInterface(){}
                virtual unsigned getNs() const = 0;
                virtual unsigned getNt() const = 0;
                virtual unsigned getNumberOfElements() const = 0;
        };

        class RootedStaggaredfieldEoParametersInterface {
            public:
                virtual ~RootedStaggaredfieldEoParametersInterface(){}
                virtual unsigned getMetropolisRationalApproximationOrder() const = 0;
                virtual unsigned getMolecularDynamicsRationalApproximationOrder() const = 0;
        };

    }
}

#include "../../meta/inputparameters.hpp"
#include "../../meta/util.hpp"

#include <iostream>

namespace physics {
    namespace lattices {

        class LatticeObjectParametersImplementation : public LatticeObjectParametersInterface {
            public:
                LatticeObjectParametersImplementation() = delete;
                LatticeObjectParametersImplementation(const meta::Inputparameters * paramsIn)
                        : parameters(paramsIn)
                {
                }
                virtual ~LatticeObjectParametersImplementation()
                {
                }
                virtual unsigned getNs() const override
                {
                    return parameters->get_nspace();
                }
                virtual unsigned getNt() const override
                {
                    return parameters->get_ntime();
                }
                virtual unsigned getPrecision() const override
                {
                    return parameters->get_precision();
                }
                virtual bool ignoreChecksumErrorsInIO() const override
                {
                    return parameters->get_ignore_checksum_errors();
                }
                virtual unsigned getNumberOfElements() const override
                {
                    return meta::get_vol4d(*parameters) * NDIM;
                }
                virtual double getKappa() const override
                {
                    return parameters->get_kappa();
                }
                virtual double getMu() const override
                {
                    return parameters->get_mu();
                }
                virtual double getBeta() const override
                {
                    return parameters->get_beta();
                }
                virtual common::startcondition getStartcondition() const override
                {
                    return parameters->get_startcondition();
                }
                virtual std::string getNamePrefix() const override
                {
                    return parameters->get_config_prefix();
                }
                virtual std::string getNamePostfix() const override
                {
                    return parameters->get_config_postfix();
                }
                virtual unsigned getNumberOfDigitsInName() const override
                {
                    return parameters->get_config_number_digits();
                }
                virtual unsigned getSmearingSteps() const override
                {
                    return parameters->get_rho_iter();
                }
                virtual std::string getSourcefileName() const override
                {
                    return parameters->get_sourcefile();
                }

            private:
                const meta::Inputparameters * parameters;
        };

        class GaugemomentaParametersImplementation final : public GaugemomentaParametersInterface {
            public:
                GaugemomentaParametersImplementation() = delete;
                GaugemomentaParametersImplementation(const meta::Inputparameters& paramsIn)
                        : parameters(paramsIn)
                {
                }
                ~GaugemomentaParametersImplementation()
                {
                }
                unsigned getNt() const override
                {
                    return parameters.get_ntime();
                }
                unsigned getNs() const override
                {
                    return parameters.get_nspace();
                }
                unsigned getNumberOfElements() const override
                {
                    return meta::get_vol4d(parameters.get_ntime(), parameters.get_nspace()) * NDIM;
                }
            private:
                const meta::Inputparameters& parameters;
        };

        class SpinorfieldParametersImplementation final : public SpinorfieldParametersInterface {
            public:
                SpinorfieldParametersImplementation() = delete;
                SpinorfieldParametersImplementation(const meta::Inputparameters& paramsIn)
                        : parameters(paramsIn)
                {
                }
                ~SpinorfieldParametersImplementation()
                {
                }
                unsigned getNt() const override
                {
                    return parameters.get_ntime();
                }
                unsigned getNs() const override
                {
                    return parameters.get_nspace();
                }
                unsigned getNumberOfElements() const override
                {
                    return meta::get_vol4d(parameters.get_ntime(), parameters.get_nspace());
                }
            private:
                const meta::Inputparameters& parameters;
        };

        class StaggaredfieldEoParametersImplementation final : public StaggaredfieldEoParametersInterface {
            public:
                StaggaredfieldEoParametersImplementation() = delete;
                StaggaredfieldEoParametersImplementation(const meta::Inputparameters& paramsIn)
                        : parameters(paramsIn)
                {
                }
                ~StaggaredfieldEoParametersImplementation()
                {
                }
                unsigned getNt() const override
                {
                    return parameters.get_ntime();
                }
                unsigned getNs() const override
                {
                    return parameters.get_nspace();
                }
                unsigned getNumberOfElements() const override
                {
                    return meta::get_vol4d(parameters.get_ntime(), parameters.get_nspace());
                }
            private:
                const meta::Inputparameters& parameters;
        };

        class RootedStaggaredfieldEoParametersImplementation final : public StaggaredfieldEoParametersInterface,
                public RootedStaggaredfieldEoParametersInterface {
            public:
                RootedStaggaredfieldEoParametersImplementation() = delete;
                RootedStaggaredfieldEoParametersImplementation(const meta::Inputparameters& paramsIn)
                        : parameters(paramsIn)
                {
                }
                ~RootedStaggaredfieldEoParametersImplementation()
                {
                }
                unsigned getNt() const override
                {
                    return parameters.get_ntime();
                }
                unsigned getNs() const override
                {
                    return parameters.get_nspace();
                }
                unsigned getMetropolisRationalApproximationOrder() const override
                {
                    return parameters.get_metro_approx_ord();
                }
                unsigned getMolecularDynamicsRationalApproximationOrder() const override
                {
                    return parameters.get_md_approx_ord();
                }
                unsigned getNumberOfElements() const override
                {
                    return meta::get_vol4d(parameters.get_ntime(), parameters.get_nspace());
                }
            private:
                const meta::Inputparameters& parameters;
        };

    }
}
