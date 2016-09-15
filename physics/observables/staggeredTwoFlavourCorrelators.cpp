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


 /*Right now this file has not the full functionality, it is still in progress*/

#include "staggeredTwoFlavourCorrelators.hpp"

#include <fstream>
#include <cmath>
#include "../lattices/staggeredfield_eo.hpp"

#include "../sources.hpp"
#include "../algorithms/inversion.hpp"

#include <cassert>
#include "../../meta/util.hpp"
#include "../lattices/util.hpp"
#include "../../hardware/device.hpp"
#include "../../hardware/code/correlator.hpp"
#include "../interfacesHandler.hpp"

static void writeCorrelatorToFile(const std::string filename, std::vector<hmc_float> correlator)
{
    //Open, write (like in Wilson) and CLOSE!!
    throw Print_Error_Message("Function writeCorrelatorToFile not implemented yet!");
}

static std::vector<physics::lattices::Staggeredfield_eo*> createAndInvertSources(const hardware::System& system, const physics::PRNG& prng, const size_t n_sources)
{
    /*
     * 1) Create the sources using the function
     *         physics::create_staggered_sources
     *    which I added for you and tested. The test fails since the point source kernel is missing!
     *
     * 2) For loop on the number of sources and for each source
     *     - Create what is needed to call CG-M used as a standard CG
     *     - Call to CG-M to invert the source. The calculation that has been done here has to reflect the
     *       calculation we did together with the even odd decomposed field. I think here you will have to
     *       to treat the even and/or the odd part only, depending on which lattice site the source was. If
     *       you do not know, come by and we see.
     *
     * 3) Return the proper object
     */
    throw Print_Error_Message("Function createAndInvertSources not implemented yet!");
}

std::vector<hmc_float> physics::observables::staggered::calculatePseudoscalarCorrelator(const std::vector<physics::lattices::Staggeredfield_eo*>& invertedSources,
                                                                                        const hardware::System& system, physics::InterfacesHandler& interfacesHandler)
{
    //Like calculate_correlator_componentwise in Wilson!
    throw Print_Error_Message("Function calculatePseudoscalarCorrelator not implemented yet!");
}

void physics::observables::staggered::measurePseudoscalarCorrelatorOnGaugefieldAndWriteToFile(const physics::lattices::Gaugefield * gaugefield,
                                                                                                    std::string currentConfigurationName,
                                                                                                    physics::InterfacesHandler & interfacesHandler)
{

    auto system = gaugefield->getSystem();
    const physics::observables::StaggeredTwoFlavourCorrelatorsParametersInterface& parametersInterface = interfacesHandler.getStaggeredTwoFlavourCorrelatorsParametersInterface();
    auto prng = gaugefield->getPrng();

    std::string filenameForCorrelatorData = parametersInterface.getCorrelatorFilename(currentConfigurationName);

    /*
     * Here we will call some functions to get the job done:
     *  - Create and invert sources -> static
     *  - Calculate the pseudoscalar correlator -> this will be in the hpp as well to be tested
     *  - Write result to file -> static
     */
}

