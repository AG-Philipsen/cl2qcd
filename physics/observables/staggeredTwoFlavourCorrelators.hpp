/* staggeredTwoFlavourCorrelators.hpp
 *
 *@file
 * Implementation of the staggered two flavour correlator calculation.
 *
 * Copyright 2016 Alessandro Sciarra, Tim Breitenfelder
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


/*Right now this file was just created for reasons of general clarity, the full content will follow soon*/

#ifndef STAGGEREDTWOFLAVOURCORRELATORS_HPP_
#define STAGGEREDTWOFLAVOURCORRELATORS_HPP_


#include "../lattices/gaugefield.hpp"
#include "../lattices/staggeredfield_eo.hpp"
#include "observablesInterfaces.hpp"
#include "../interfacesHandler.hpp"

namespace physics{
	namespace observables{
		namespace staggered{

		    std::vector<hmc_float> calculatePseudoscalarCorrelator(const std::pair<std::vector<physics::lattices::Staggeredfield_eo*>, std::vector<physics::lattices::Staggeredfield_eo*> >&,
		                                                           physics::InterfacesHandler&);

		    void measurePseudoscalarCorrelatorOnGaugefieldAndWriteToFile(const physics::lattices::Gaugefield&, std::string, physics::InterfacesHandler&);

		}
	}
}

#endif /* STAGGEREDTWOFLAVOURCORRELATORS_HPP_ */
