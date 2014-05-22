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

namespace physics{

	namespace observables{

		namespace wilson{

			class TwoFlavourChiralCondensate
			{
			public:
			        TwoFlavourChiralCondensate(const physics::lattices::Gaugefield * gaugefield, std::string configurationName = "what");
				TwoFlavourChiralCondensate() = delete;
				~TwoFlavourChiralCondensate();
				
				std::vector<double> getChiralCondensate();
				void measureChiralCondensate(const physics::lattices::Gaugefield * gaugefield);
				void writeChiralCondensateToFile();
				
			private:
				const physics::lattices::Gaugefield * gaugefield;
				const meta::Inputparameters * parameters;
				const hardware::System * system;
				const physics::PRNG * prng;
				int trajectoryNumber;
				std::vector<double> chiralCondensate;
				std::ofstream outputToFile;
				std::string filenameForChiralCondensateData;
			        std::string configurationName;
				
				void checkInputparameters();
				double norm_std() const ;
			        double norm_tm() const;
				double flavourChiralCondensate_std(const physics::lattices::Spinorfield* phi, const physics::lattices::Spinorfield* xi);
			        double flavour_doublet_chiral_condensate_tm(const physics::lattices::Spinorfield* phi);
				void openFileForWriting();
				void flavour_doublet_chiral_condensate(const physics::lattices::Spinorfield* inverted, const physics::lattices::Spinorfield* sources);
			};

		  std::vector<double> measureChiralCondensateAndWriteToFile(const physics::lattices::Gaugefield * gaugefield, std::string configurationName);
    }
  }
}
#endif
