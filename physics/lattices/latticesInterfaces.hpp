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
#include <string>

namespace physics {
    namespace lattices {

        class GaugefieldParametersInterface {
            public:
                virtual ~GaugefieldParametersInterface(){}
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

        class SpinorfieldEoParametersInterface {
            public:
                virtual ~SpinorfieldEoParametersInterface() = 0;
        };
        //Pure virtual destructors must be implemented outside the class! (inline for multiple inclusion of header)
        inline SpinorfieldEoParametersInterface::~SpinorfieldEoParametersInterface(){}

        class StaggeredfieldEoParametersInterface {
            public:
                virtual ~StaggeredfieldEoParametersInterface(){}
                virtual unsigned getNumberOfElements() const = 0;
        };

        class RootedStaggeredfieldEoParametersInterface : public StaggeredfieldEoParametersInterface {
            public:
                virtual ~RootedStaggeredfieldEoParametersInterface(){}
                virtual unsigned getMetropolisRationalApproximationOrder() const = 0;
                virtual unsigned getMolecularDynamicsRationalApproximationOrder() const = 0;
                virtual unsigned getNumberOfPseudofermions() const = 0;
        };

    }
}

