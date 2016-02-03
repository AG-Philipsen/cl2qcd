/** @file
 * algorithms interfaces declaration
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

#pragma once

#include "../../common_header_files/types.h"

namespace physics{
	namespace algorithms{

	    class SolversParametersInterface {
	        public:
	            virtual ~SolversParametersInterface() {}
	            virtual unsigned getCgMax() const = 0;
	            virtual unsigned getIterRefresh() const = 0;
	            virtual common::solver getSolver() const = 0;
	            virtual unsigned getCgIterationBlockSize() const = 0;
	            virtual unsigned getCgMinimumIterationCount() const = 0;
	            virtual bool getCgUseAsyncCopy() const = 0;
	            virtual bool getUseMergeKernelsSpinor() const = 0;
	    };

        class ForcesParametersInterface{
            public:
                virtual ~ForcesParametersInterface(){}
                virtual common::action getFermact() const = 0;
                virtual double getForcePreconditioning() const = 0;
                virtual unsigned getRhoIterations() const = 0;
                virtual common::solver getSolver() const = 0;
                virtual bool getUseSmearing() const = 0;
                virtual bool getUseGaugeOnly() const = 0;
                virtual bool getUseRectangles() const = 0;
        };

        class MinMaxEigenvalueParametersInterface {
            public:
                virtual ~MinMaxEigenvalueParametersInterface(){}
                virtual unsigned getFindMinMaxIterationBlockSize() const = 0;
                virtual unsigned getFindMinMaxMaxValue() const = 0;
        };

		class InversionParemetersInterface{
            public:
                virtual ~InversionParemetersInterface(){}
                virtual common::action getFermact() const = 0;
                virtual common::solver getSolver() const = 0;
                virtual double getSolverPrec() const = 0;
                virtual bool getUseEo() const = 0;
                virtual bool getUseSmearing() const = 0;
		};

		class IntegratorParametersInterface{
            public:
                virtual ~IntegratorParametersInterface(){}
                virtual unsigned getIntegrationSteps(const unsigned) const = 0;
                virtual common::integrator getIntegrator(const unsigned) const = 0;
                virtual double getLambda(const unsigned) const = 0;
                virtual unsigned getNumTimescales() const = 0;
                virtual double getTau() const = 0;
                virtual bool getUseMp() const = 0;
		};

        class MolecularDynamicsInterface{
            public:
                virtual ~MolecularDynamicsInterface(){}
                virtual double getSolverPrec() const = 0;
        };

		class MetropolisParametersInterface{
            public:
                virtual ~MetropolisParametersInterface(){}
                virtual double getC0() const = 0;
                virtual double getC1() const = 0;
                virtual common::action getFermact() const = 0;
                virtual size_t getRectanglesNormalization() const = 0;
                virtual common::solver getSolver() const = 0;
                virtual double getSolverPrec() const = 0;
                virtual bool getUseGaugeOnly() const = 0;
                virtual bool getUseMp() const = 0;
                virtual bool getUseRectangles() const = 0;
		};

		class HmcParametersInterface{
            public:
                virtual ~HmcParametersInterface(){}
                virtual double getBeta() const = 0;
                virtual bool getUseEo() const = 0;
                virtual bool getUseGaugeOnly() const = 0;
                virtual bool getUseMp() const = 0;
		};

		class RhmcParametersInterface{
            public:
                virtual ~RhmcParametersInterface(){}
                virtual double getBeta() const = 0;
                virtual bool getConservative() const = 0;
                virtual double getFindMinMaxPrec() const = 0;
                virtual double getMass() const = 0;
                virtual bool getUseGaugeOnly() const = 0;
                virtual bool getUseMp() const = 0;
                virtual bool getUseEo() const = 0;
		};

	}
}

