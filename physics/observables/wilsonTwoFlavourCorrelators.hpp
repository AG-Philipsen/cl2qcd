/** @file
 * physics::observables::wilson::TwoFlavourCorrelators class
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

#ifndef WILSONTWOFLAVOURCORRELATORS_HPP_
#define WILSONTWOFLAVOURCORRELATORS_HPP_

#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"
#include <fstream>
#include <cmath>
#include "../../meta/inputparameters.hpp"

namespace physics{

	namespace observables{

		namespace wilson{

			class TwoFlavourCorrelators
			{
			public:
				
			private:
				const physics::lattices::Gaugefield * gaugefield;
				const meta::Inputparameters * parameters;
				const hardware::System * system;
				const physics::PRNG * prng;
				int trajectoryNumber;
				std::vector<double> correlator;
				std::ofstream outputToFile;
				std::string filenameForCorrelatorData;
				std::string configurationName;
			};
    }
  }
}
#endif
