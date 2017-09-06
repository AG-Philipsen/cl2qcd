/** @file
 * Interfaces for physics::observables
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

namespace physics {
    namespace observables {

        class GaugeObservablesParametersInterface {
            public:
                virtual ~GaugeObservablesParametersInterface(){}
                virtual bool measureRectangles() const = 0;
                virtual bool measureTransportCoefficientKappa() const = 0;
                virtual bool printToScreen() const = 0;
                virtual hmc_float getBeta() const = 0;
                virtual std::string getTransportCoefficientKappaFilename() const = 0;
                virtual std::string getRectanglesFilename() const = 0;
                virtual std::string getGaugeObservablesFilename(std::string) const = 0;
                virtual unsigned getTemporalPlaquetteNormalization() const = 0;
                virtual unsigned getSpatialPlaquetteNormalization() const = 0;
                virtual unsigned getPlaquetteNormalization() const = 0;
                virtual unsigned getSpatialVolume() const = 0;
                virtual unsigned getPolyakovLoopNormalization() const = 0;
        };

        class WilsonTwoFlavourChiralCondensateParametersInterface {
            public:
                virtual ~WilsonTwoFlavourChiralCondensateParametersInterface(){}
                virtual common::action getFermionicActionType() const = 0;
                virtual common::pbp_version getPbpVersion() const = 0;
                virtual bool measurePbp() const = 0;
                virtual int getNumberOfSources() const = 0;
                virtual hmc_float getKappa() const = 0;
                virtual hmc_float getMubar() const = 0;
                virtual unsigned get4dVolume() const = 0;
                virtual bool useEvenOdd() const = 0;
                virtual std::string getPbpFilename(std::string configurationName) const = 0;
        };

        class StaggeredChiralCondensateParametersInterface {
            public:
                virtual ~StaggeredChiralCondensateParametersInterface(){}
                virtual bool measurePbp() const = 0;
                virtual int getNumberOfSources() const = 0;
                virtual hmc_float getMass() const = 0;
                virtual hmc_float getSolverPrecision() const = 0;
                virtual hmc_float getNumberOfTastes() const = 0;
                virtual unsigned get4dVolume() const = 0;
                virtual std::string getPbpFilename(std::string configurationName) const = 0;
                virtual unsigned getPbpNumberOfMeasurements() const = 0;
        };

        class WilsonTwoFlavourCorrelatorsParametersInterface {
            public:
                virtual ~WilsonTwoFlavourCorrelatorsParametersInterface(){}
                virtual bool printToScreen() const = 0;
                virtual void printInformationOfFlavourDoubletCorrelator(std::ostream* of = nullptr) const = 0;
                virtual unsigned getCorrelatorDirection() const = 0;
                virtual common::sourcetypes getSourceType() const = 0;
                virtual unsigned getNs() const = 0;
                virtual unsigned getNt() const = 0;
                virtual std::string getCorrelatorFilename(std::string currentConfigurationName) const = 0;
                virtual bool placeSourcesOnHost() const = 0;
                virtual int getNumberOfSources() const = 0;
        };

        class StaggeredTwoFlavourCorrelatorsParametersInterface {
			public:
				virtual ~StaggeredTwoFlavourCorrelatorsParametersInterface(){}
				virtual bool printToScreen() const = 0;
				virtual void printInformationOfFlavourDoubletCorrelator(std::ostream* of = nullptr) const = 0;
//                virtual unsigned getCorrelatorDirection() const = 0;
                virtual common::sourcetypes getSourceType() const = 0;
//                virtual unsigned getNs() const = 0;
                virtual unsigned getNt() const = 0;
				virtual std::string getCorrelatorFilename(std::string currentConfigurationName) const = 0;
				virtual hmc_float getSolverPrecision() const = 0;
//                virtual bool placeSourcesOnHost() const = 0;
                virtual int getNumberOfSources() const = 0;
        };


    }
}
