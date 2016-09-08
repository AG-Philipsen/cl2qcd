/* staggeredTwoFlavourCorrelators.cpp
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


 /*Right now this file was just created for reasons of clarity, the full content will follow soon*/

#include "staggeredTwoFlavourCorrelators.hpp"

#include <fstream>
#include <cmath>
#include "../lattices/staggeredfield_eo.hpp"

#include "../sources.hpp"
#include "../algorithms/inversion.hpp"

#include <cassert>
#include "../../meta/util.hpp"
#include "../lattices/util.hpp"
#include "../lattices/swappable.hpp"
#include "../../hardware/device.hpp"
#include "../../hardware/code/correlator.hpp"
#include "../interfacesHandler.hpp"

static void calculate_correlator(const std::string& type, const std::vector<const hardware::buffers::Plain<hmc_float>*>& results, physics::lattices::Staggeredfield_eo* corr,
        physics::lattices::Staggerdfield_eo* source, const hardware::System& system,
        const physics::observables::StaggeredTwoFlavourCorrelatorsParametersInterface& parametersInterface, physics::InterfacesHandler& interfacesHandler)
{

	throw Print_Error_Message("method implementation is in process but not yet finished!");

}

static void calculate_correlator(const std::string& type, const std::vector<const hardware::buffers::Plain<hmc_float>*>& results,
                                 physics::lattices::Staggeredfield_eo* corr1, physics::lattices::Staggeredfield_eo* source1,
                                 physics::lattices::Staggeredfield_eo* corr2, physics::lattices::Staggeredfield_eo* source2,
                                 physics::lattices::Staggeredfield_eo* corr3, physics::lattices::Staggeredfield_eo* source3,
                                 const physics::observables::StaggeredTwoFlavourCorrelatorsParametersInterface& params)
{

	throw Print_Error_Message("method implementation is in process but not yet finished!");

}






void physics::observables::staggered::measureTwoFlavourCorrelatorsOnGaugefield(const physics::lattices::Gaugefield * gaugefield, std::string currentConfigurationName, physics::InterfacesHandler & interfacesHandler)
{

	throw Print_Error_Message("method implementation is in process but not yet finished!");

}

