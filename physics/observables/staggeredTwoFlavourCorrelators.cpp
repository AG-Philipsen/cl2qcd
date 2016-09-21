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
#include <cassert>

#include "../lattices/staggeredfield_eo.hpp"
#include "../sources.hpp"
#include "../algorithms/inversion.hpp"
#include "../../meta/util.hpp"
#include "../lattices/util.hpp"
#include "../../hardware/device.hpp"
#include "../../hardware/code/correlator.hpp"
#include "../interfacesHandler.hpp"
#include "../algorithms/solver_shifted.hpp"

static void writeCorrelatorToFile(const std::string filename, std::vector<hmc_float> correlator)
{
    //Open, write (like in Wilson) and CLOSE!!
    throw Print_Error_Message("Function writeCorrelatorToFile not implemented yet!");
}

static std::vector<physics::lattices::Staggeredfield_eo*> createAndInvertSources(const hardware::System& system, const physics::PRNG& prng, const size_t numberOfSources,
																				 physics::InterfacesHandler & interfacesHandler)
{
	// step 1a)
	// creating the sources

	const std::vector<physics::lattices::Staggeredfield_eo *>  sources = physics::create_staggered_sources(system, prng, numberOfSources, interfacesHandler);

    // step 1b)
    // Ask the InterfaceHandler for the sourceInterface
	// to decide whether the point source is on a even or odd site

    const physics::SourcesParametersInterface & sourcesParameters = interfacesHandler.getSourcesParametersInterface();

    const int tpos = SourcesParametersInterface.getSourceT(sources);
    const int xpos = SourcesParametersInterface.getSourceX(sources);
    const int ypos = SourcesParametersInterface.getSourceY(sources);
    const int zpos = SourcesParametersInterface.getSourceZ(sources);

    bool sourceOnEvenSite;
    if((tpos+xpos+ypos+zpos)%2 == 0)
    {
    	sourceOnEvenSite = TRUE;
    }
    else
    {
    	sourceOnEvenSite = FALSE;
    }

    // step 2):
    // invert every single source

    const physics::AdditionalParameters& additionalParameters = interfacesHandler.getAdditionalParameters<physics::lattices::Staggeredfield_eo>();
    const hmc_float mass = additionalParameters.getMass();
    std::vector<hmc_float> sigma(1, 0.0);

    for (unsigned int interationsNumb = 0; iterationsNumb < numberOfSources; ++iterationsNumb)
    {

    	const std::vector<physics::lattices::Staggeredfield_eo *>  source = sources[iterationsNumb];


    	if (sourceOnEvenSite == TRUE)
    	{
			// Case Even
			// upper_left =  EVEN
    		// creator of MdagM even or odd unclear

			// calculate phi_e = ((MdagM)_ee)^-1 * S_e
			physics::fermionmatrix::MdagM_eo MdagM(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>(), upper_left);
			std::vector<std::shared_ptr<Staggeredfield_eo> > phi_e; //This is the type to be used in the inverter
			cg_m(phi_e, MdagM, gf, sigma, source, system, interfacesHandler, parametersInterface.getSolverPrecision(), additionalParameters);

			// calculate chi_e = m * phi_e
			Staggeredfield_eo chi_e(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
			physics::lattices::sax(&chi_e, mass, phi_e);

			// calculate chi_o = - D_oe * phi_e
			Staggeredfield_eo chi_o(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
			physics::fermionmatrix::D_KS_eo Doe(system, interfacesHandler.getInterface<physics::fermionmatrix::D_KS_eo>(), ODD);
			Doe(&chi_o, gf, phi_e);
			physics::lattices::sax(&chi_o, -1.0 , chi_o);

			// overwrite inverted source for output
			sources[iterationsNumb]=chi_o;
    	}
    	else
    	{
			// Case Odd:
			// upper_left = ODD
    		// creator of MdagM even or odd unclear

			// calculate phi_o = ((MdagM)_oo)^-1 * S_o
			physics::fermionmatrix::MdagM_eo MdagM(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>(), upper_left);
			std::vector<std::shared_ptr<Staggeredfield_eo> > phi_o; //This is the type to be used in the inverter
			cg_m(phi_o, MdagM, gf, sigma, source, system, interfacesHandler, parametersInterface.getSolverPrecision(), additionalParameters);

			// calculate chi_e = -D_eo * phi_o
			Staggeredfield_eo chi_o(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
			physics::lattices::sax(&chi_o, mass, phi_o);

			// calculate chi_o = m * phi_o
			Staggeredfield_eo chi_e(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
			physics::fermionmatrix::D_KS_eo Deo(system, interfacesHandler.getInterface<physics::fermionmatrix::D_KS_eo>(), EVEN);
			Deo(&chi_e, gf, phi_o);
			physics::lattices::sax(&chi_e, -1.0 , chi_e);

			// overwrite inverted source for output
			sources[iterationsNumb]=chi_e;
    	}

    }


    // throw Print_Error_Message("Function createAndInvertSources not implemented yet!");

}

//-----------------------------------------------------------------------------------------------------------------------

// List of Tasks:

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

//------------------------------------------------------------------------------------------------------------------------


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

    /* static std::vector<physics::lattices::Staggeredfield_eo*> invertedSources;
     * invertedSources = createAndInvertSources(system, prng, numberOfSources, interfacesHandler);
     *
     * std::vector<hmc_float> correlator;
     * correlator = physics::observables::staggered::calculatePseudoscalarCorrelator(invertedSources, system, interfacesHandler);
     *
     * const std::string filename = {'PS-staggered_Pion Correlator Measurement Data'}
     * writeCorrelatorToFile(filename, correlator);
     */



    /*
     * Here we will call some functions to get the job done:
     *  - Create and invert sources -> static
     *  - Calculate the pseudoscalar correlator -> this will be in the hpp as well to be tested
     *  - Write result to file -> static
     */
}

