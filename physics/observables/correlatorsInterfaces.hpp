/** @file
 * Interfaces for physics::observables
 *
 * Copyright 2015 Christopher Czaban
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

#include <cstring>
#include "../../meta/inputparameters.hpp"
#include "../../meta/util.hpp"
#include "../../meta/parametersSources.hpp"


namespace physics {
	namespace observables {

		class WilsonTwoFlavourCorrelatorsParametersInterface {
			public:
				virtual ~WilsonTwoFlavourCorrelatorsParametersInterface() {};
				virtual bool printToScreen() const = 0;
				virtual unsigned getCorrelatorDirection() const = 0;
				virtual meta::ParametersSources::sourcetypes getSourceType() const = 0;
				virtual unsigned getNs() const = 0;
				virtual unsigned getNt() const = 0;
				virtual std::string getFermionObservablesCorrleatorFileName(std::string currentConfigurationName) const = 0;
				virtual bool getPlaceOfSourcesOnHost() const = 0;
				virtual	int getNumberOfSources() const = 0;
		};

		class WilsonTwoFlavorCorrelatorsParametersImplementation final: public WilsonTwoFlavourCorrelatorsParametersInterface{
			public:
				WilsonTwoFlavorCorrelatorsParametersImplementation() = delete;
				WilsonTwoFlavorCorrelatorsParametersImplementation(meta::Inputparameters& parametersIn) : parameters(parametersIn) {};
				~WilsonTwoFlavorCorrelatorsParametersImplementation();
				bool printToScreen () const override
				{
					return parameters.get_print_to_screen();
				}
				unsigned getCorrelatorDirection() const override
				{
					return parameters.get_corr_dir();
				}
				meta::ParametersSources::sourcetypes getSourceType() const override
				{
					return parameters.get_sourcetype();
				}
				unsigned getNs() const override
				{
					return parameters.get_nspace();
				}
				unsigned getNt() const override
				{
					return parameters.get_ntime();
				}
				std::string getFermionObservablesCorrelatorFileName(std::string currentConfigurationName) const override
				{
					return meta::get_ferm_obs_corr_file_name(parameters, currentConfigurationName);
				}
				bool getPlaceOfSourcesOnHost() const override
				{
					return parameters.get_place_sources_on_host();
				}
				int getNumberOfSources() const override
				{
					return parameters.get_num_sources();
				}

			private:
				const meta::Inputparameters& parameters;
		};
	}
}
