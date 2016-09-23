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
#include "../../hardware/code/correlator_staggered.hpp"

static void writeCorrelatorToFile(const std::string filename, std::vector<hmc_float> correlator,
                                  const physics::observables::StaggeredTwoFlavourCorrelatorsParametersInterface& parametersInterface)
{
    //TODO: Remove this use namespace!!
    using namespace std;

	ofstream of;
	of.open (filename.c_str(), ios_base::app);
	if(!of.is_open())
	{
		throw File_Exception(filename);
	}

	if(parametersInterface.printToScreen())
		parametersInterface.printInformationOfFlavourDoubletCorrelator();

	parametersInterface.printInformationOfFlavourDoubletCorrelator(&of);

	logger.info() << "staggered pseudoscalar pion correlator:" ;
	for(size_t j = 0; j < correlator.size(); j++)
	{
		logger.info() << j << "\t" << scientific << setprecision(14) << correlator[j];
		of << scientific << setprecision(14) << "0 1\t" << j << "\t" << correlator[j] << endl;
	}

	of << endl;
	of.close();

}

static std::pair<std::vector<physics::lattices::Staggeredfield_eo*>, std::vector<physics::lattices::Staggeredfield_eo*> >
createAndInvertSources(const hardware::System& system, const physics::PRNG& prng, const physics::lattices::Gaugefield& gaugefield,
                       const int numberOfSources, physics::InterfacesHandler& interfacesHandler)
{
	// step 1a)
	// creating the sources

	const std::vector<physics::lattices::Staggeredfield_eo *>  sources = physics::create_staggered_sources(system, prng, numberOfSources, interfacesHandler);

    // step 1b)
    // Ask the InterfaceHandler for the sourceInterface
	// to decide whether the point source is on a even or odd site

    const physics::SourcesParametersInterface & sourcesParameters = interfacesHandler.getSourcesParametersInterface();

    const int tpos = sourcesParameters.getSourceT();
    const int xpos = sourcesParameters.getSourceX();
    const int ypos = sourcesParameters.getSourceY();
    const int zpos = sourcesParameters.getSourceZ();

    const bool sourceOnEvenSite = ((tpos+xpos+ypos+zpos)%2 == 0) ? true : false;

    // step 2):
    // invert every single source

    std::vector<physics::lattices::Staggeredfield_eo*> invertedSourcesEvenParts;
    std::vector<physics::lattices::Staggeredfield_eo*> invertedSourcesOddParts;
    const physics::AdditionalParameters& additionalParameters = interfacesHandler.getAdditionalParameters<physics::lattices::Staggeredfield_eo>();
    const physics::observables::StaggeredTwoFlavourCorrelatorsParametersInterface& staggeredTwoFlavourCorrelatorsParametersInterface = interfacesHandler.getStaggeredTwoFlavourCorrelatorsParametersInterface();
    const hmc_float mass = additionalParameters.getMass();
    std::vector<hmc_float> sigma(1, 0.0);

    for (int iterationNumber = 0; iterationNumber < numberOfSources; ++iterationNumber)
    {
    	if(sourceOnEvenSite)
    	{
    	    bool upper_left =  EVEN;

			// calculate phi_e = ((MdagM)_ee)^-1 * S_e
			physics::fermionmatrix::MdagM_eo MdagM(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>(), upper_left);
			std::vector<std::shared_ptr<physics::lattices::Staggeredfield_eo> > phi_e; //This is the type to be used in the CG-M
			phi_e.emplace_back(std::make_shared<physics::lattices::Staggeredfield_eo>(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>()));
			physics::algorithms::solvers::cg_m(phi_e, MdagM, gaugefield, sigma, *(sources[iterationNumber]), system, interfacesHandler, staggeredTwoFlavourCorrelatorsParametersInterface.getSolverPrecision(), additionalParameters);

			// calculate chi_e = m * phi_e
			physics::lattices::Staggeredfield_eo* chi_e = new physics::lattices::Staggeredfield_eo(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
			physics::lattices::sax(chi_e, mass, *(phi_e[0]));

			// calculate chi_o = - D_oe * phi_e
			physics::lattices::Staggeredfield_eo* chi_o = new physics::lattices::Staggeredfield_eo(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
			physics::fermionmatrix::D_KS_eo Doe(system, interfacesHandler.getInterface<physics::fermionmatrix::D_KS_eo>(), ODD);
			Doe(chi_o, gaugefield, *(phi_e[0]));
			physics::lattices::sax(chi_o, -1.0 , *chi_o);

            // set inverted sources for output
            invertedSourcesEvenParts.push_back(chi_e);
            invertedSourcesOddParts.push_back(chi_o);
    	}
    	else
    	{
    	    bool upper_left = ODD;

			// calculate phi_o = ((MdagM)_oo)^-1 * S_o
    	    physics::fermionmatrix::MdagM_eo MdagM(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>(), upper_left);
    	    std::vector<std::shared_ptr<physics::lattices::Staggeredfield_eo> > phi_o; //This is the type to be used in the CG-M
    	    phi_o.emplace_back(std::make_shared<physics::lattices::Staggeredfield_eo>(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>()));
    	    physics::algorithms::solvers::cg_m(phi_o, MdagM, gaugefield, sigma, *(sources[iterationNumber]), system, interfacesHandler, staggeredTwoFlavourCorrelatorsParametersInterface.getSolverPrecision(), additionalParameters);

    	    // calculate chi_e = -D_eo * phi_o
    	    physics::lattices::Staggeredfield_eo* chi_e = new physics::lattices::Staggeredfield_eo(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
            physics::fermionmatrix::D_KS_eo Deo(system, interfacesHandler.getInterface<physics::fermionmatrix::D_KS_eo>(), EVEN);
            Deo(chi_e, gaugefield, *(phi_o[0]));
            physics::lattices::sax(chi_e, -1.0 , *chi_e);

			// calculate chi_o = m * phi_o
            physics::lattices::Staggeredfield_eo* chi_o = new physics::lattices::Staggeredfield_eo(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
            physics::lattices::sax(chi_o, mass, *(phi_o[0]));

            // set inverted sources for output
            invertedSourcesEvenParts.push_back(chi_e);
            invertedSourcesOddParts.push_back(chi_o);
    	}

    }

    if((int)invertedSourcesEvenParts.size() != numberOfSources || (int)invertedSourcesOddParts.size() != numberOfSources)
        throw Print_Error_Message("Something went really wrong in \"createAndInvertSources\" since the number of inverted sources does not match the number of sources!");


    physics::lattices::release_staggeredfields_eo(sources);
    return std::make_pair(invertedSourcesEvenParts, invertedSourcesOddParts);

}



std::vector<hmc_float> physics::observables::staggered::calculatePseudoscalarCorrelator(
        const std::pair<std::vector<physics::lattices::Staggeredfield_eo*>, std::vector<physics::lattices::Staggeredfield_eo*> >& invertedSources,
        const hardware::System& system, physics::InterfacesHandler& interfacesHandler)
{
	//Like calculate_correlator_componentwise in Wilson!
	//Until now, system not used in this function!


	const physics::observables::StaggeredTwoFlavourCorrelatorsParametersInterface& staggeredTwoFlavourCorrelatorsParametersInterface = interfacesHandler.getStaggeredTwoFlavourCorrelatorsParametersInterface();
    const std::vector<physics::lattices::Staggeredfield_eo*> invertedSourcesOnEvenSites = invertedSources.first;
    const std::vector<physics::lattices::Staggeredfield_eo*> invertedSourcesOnOddSites = invertedSources.second;

    if(invertedSourcesOnEvenSites.size()!=invertedSourcesOnOddSites.size())
    {
    	throw Print_Error_Message("# of invertedSourcesOnEvenSites is unequal to # of invertedSourcesOnOddSites, big mistake in invert sources");
    }

    auto firstInvertedevenSource = invertedSourcesOnEvenSites.at(0);
    auto evenSiteBuffers = firstInvertedevenSource->get_buffers();

    auto firstInvertedOddSource = invertedSourcesOnOddSites.at(0);
    auto oddSiteBuffers = firstInvertedOddSource->get_buffers();

    //Test on single device
	if(evenSiteBuffers.size() != 1)
	{
		throw Print_Error_Message("Multiple GPU calculation not available in staggered case");
	}

	if(oddSiteBuffers.size() != 1)
	{
			throw Print_Error_Message("Multiple GPU calculation not available in staggered case");
	}

	static size_t numberOfCorrelatorEntries = staggeredTwoFlavourCorrelatorsParametersInterface.getNt();

	//Allocate result (and not results) as const hardware::buffers::Plain<hmc_float>*
	const hardware::buffers::Plain<hmc_float>* correlatorResult;

	const hardware::Device *device = evenSiteBuffers->get_device();
	correlatorResult = new hardware::buffers::Plain<hmc_float>(numberOfCorrelatorEntries, device);
	correlatorResult->clear();

//	if(correlatorResult.getSourceType() != common::sourcetypes::point)
//	{
//		throw Print_Error_Message("staggered correlator calculation only for pointsources implemented");
//	}

	//for loop with calculate correlator inside

	//const std::string& type = "ps";  //other pions not implemented yet.


	for(size_t i = 0; i < invertedSourcesOnEvenSites.size(); i++)
	{

		auto evenSiteinvertedSource = invertedSourcesOnEvenSites.at(i);
		auto evenSiteinvertedSourceBuffer = evenSiteinvertedSource->get_buffers();
		auto oddSiteinvertedSource = invertedSourcesOnOddSites.at(i);
		auto oddSiteinvertedSourceBuffer = oddSiteinvertedSource->get_buffers();

		auto code = correlatorResult->get_device()->getCorrelatorStaggeredCode();
		//hardware::code::Correlator_staggered::pseudoScalarCorrelator
		code->pseudoScalarCorrelator(correlatorResult, evenSiteinvertedSourceBuffer[i]);
		code->pseudoScalarCorrelator(correlatorResult, oddSiteinvertedSourceBuffer[i]);

	}

	//Create host_result, dump the data from device directly into it without loop and sum (+=)
	// Tim : I think the loop over the number of correlator entries is needed, but not the loop over the number of buffers...

	std::vector<hmc_float> host_correlatorResult(numberOfCorrelatorEntries);
	for(size_t i = 0; i < numberOfCorrelatorEntries; ++i)
	{
		host_correlatorResult[i] = 0.;
	}

	auto result = correlatorResult;
	std::vector<hmc_float> out(numberOfCorrelatorEntries);
	result->dump(out.data());
	for(size_t i = 0; i < numberOfCorrelatorEntries; ++i)
	{
		logger.trace() << out[i];
		host_correlatorResult[i] += out[i];
	}
	delete result;
	return host_correlatorResult;

    //throw Print_Error_Message("Function calculatePseudoscalarCorrelator not implemented yet!");
}

void physics::observables::staggered::measurePseudoscalarCorrelatorOnGaugefieldAndWriteToFile(
		const physics::lattices::Gaugefield& gaugefield,std::string currentConfigurationName,
		physics::InterfacesHandler & interfacesHandler)
{

    auto system = gaugefield.getSystem();
    const physics::observables::StaggeredTwoFlavourCorrelatorsParametersInterface& parametersInterface = interfacesHandler.getStaggeredTwoFlavourCorrelatorsParametersInterface();
    auto prng = gaugefield.getPrng();

    /*
     * TODO: Make test on correlator direction and throw if not t-dir.
     */

    std::string filenameForCorrelatorData = parametersInterface.getCorrelatorFilename(currentConfigurationName);

    std::pair<std::vector<physics::lattices::Staggeredfield_eo*>, std::vector<physics::lattices::Staggeredfield_eo*> > invertedSources;
    invertedSources = createAndInvertSources(*system, *prng, gaugefield, parametersInterface.getNumberOfSources(), interfacesHandler);

    std::vector<hmc_float> correlator = physics::observables::staggered::calculatePseudoscalarCorrelator(invertedSources, *system, interfacesHandler);

    writeCorrelatorToFile(filenameForCorrelatorData, correlator, parametersInterface);

    physics::lattices::release_staggeredfields_eo(invertedSources.first);
    physics::lattices::release_staggeredfields_eo(invertedSources.second);

    /*
     * Here we will call some functions to get the job done:
     *  - Create and invert sources -> static
     *  - Calculate the pseudoscalar correlator -> this will be in the hpp as well to be tested
     *  - Write result to file -> static
     */
}

