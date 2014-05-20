/** @file
 * physics::observables::wilson::TwoFlavourChiralCondensate class
 *
 * Copyright 2014 Christopher Pinke
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

#ifndef WILSONTWOFLAVOURCHIRALCONDENSATE_HPP_
#define WILSONTWOFLAVOURCHIRALCONDENSATE_HPP_

#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"
#include <fstream>
#include <cmath>
#include "../../meta/inputparameters.hpp"

#include "../algorithms/inversion.hpp"
#include "../algorithms/flavour_doublet.hpp"
#include "../sources.hpp"

namespace physics{

	namespace observables{

		namespace wilson{

			class TwoFlavourChiralCondensate
			{
			public:
				TwoFlavourChiralCondensate(const meta::Inputparameters * parametersIn)
				{
	  			parameters = parametersIn;
					if(! parameters->get_measure_pbp() ) 
					{
						throw std::logic_error("Chiral condensate calculation disabled in parameter setting. Aborting...");
					}
				}
				TwoFlavourChiralCondensate() = delete;

				//perhaps move this out of class?
				double measureChiralCondensate(const physics::lattices::Gaugefield * gaugefield, int iteration)
				{
					auto system = gaugefield->getSystem();
					auto prng = gaugefield->getPrng();

					std::string currentConfigurationName = "replace";
					filenameForChiralCondensateData = meta::get_ferm_obs_pbp_file_name(*parameters, currentConfigurationName);
					int sourceNumber = 0;

					for (; sourceNumber < parameters->get_num_sources(); sourceNumber++) {
						auto sources = physics::create_sources(*system, *prng, 1);
						auto result = physics::lattices::create_spinorfields(*system, sources.size());
						physics::algorithms::perform_inversion(&result, gaugefield, sources, *system);
						physics::algorithms::flavour_doublet_chiral_condensate(result, sources, filenameForChiralCondensateData, gaugefield->get_parameters_source().trajectorynr_source, *system);
						physics::lattices::release_spinorfields(result);
						physics::lattices::release_spinorfields(sources);
					}

					return 0.;
				}
	
			private:
				const meta::Inputparameters * parameters;
				double chiralCondensate;
				std::ofstream outputToFile;
				std::string filenameForChiralCondensateData;
				
				void writeChiralCondensateToFile(int iter,  const std::string& filename);
      };
    }
  }
}
#endif
