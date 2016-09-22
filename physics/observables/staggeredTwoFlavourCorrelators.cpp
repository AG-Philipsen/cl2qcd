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

	//How does ParameterInterface know that we are in StaggeredTwoFlavourCorrelatorsParametersInterface ??
	//parametersInterface.printInformationOfFlavourDoubletCorrelator(&of);



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


	// Wilson code:
	/*
	 * static void flavour_doublet_correlators(const std::vector<physics::lattices::Spinorfield*>& result, const std::vector<physics::lattices::Spinorfield*>& sources,
                                        std::string corr_fn, const hardware::System& system,
                                        const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& parametersInterface, physics::InterfacesHandler & interfacesHandler)
		{
			using namespace std;

			ofstream of(corr_fn.c_str(), ios_base::app);
			if(!of.is_open()) {
			  throw File_Exception(corr_fn);
			}

			auto result_ps = physics::observables::wilson::calculate_correlator("ps", result, sources, system, interfacesHandler);
			auto result_sc = physics::observables::wilson::calculate_correlator("sc", result, sources, system, interfacesHandler);
			auto result_vx = physics::observables::wilson::calculate_correlator("vx", result, sources, system, interfacesHandler);
			auto result_vy = physics::observables::wilson::calculate_correlator("vy", result, sources, system, interfacesHandler);
			auto result_vz = physics::observables::wilson::calculate_correlator("vz", result, sources, system, interfacesHandler);
			auto result_ax = physics::observables::wilson::calculate_correlator("ax", result, sources, system, interfacesHandler);
			auto result_ay = physics::observables::wilson::calculate_correlator("ay", result, sources, system, interfacesHandler);
			auto result_az = physics::observables::wilson::calculate_correlator("az", result, sources, system, interfacesHandler);

			if(parametersInterface.printToScreen())
				parametersInterface.printInformationOfFlavourDoubletCorrelator();

			parametersInterface.printInformationOfFlavourDoubletCorrelator(&of);

			// @todo One could also implement to write all results on screen if wanted
			//the pseudo-scalar (J=0, P=1)
			logger.info() << "pseudo scalar correlator:" ;
			for(size_t j = 0; j < result_ps.size(); j++) {
				logger.info() << j << "\t" << scientific << setprecision(14) << result_ps[j];
				of << scientific << setprecision(14) << "0 1\t" << j << "\t" << result_ps[j] << endl;
			}

			//the scalar (J=0, P=0)
			for(size_t j = 0; j < result_sc.size(); j++) {
				of << scientific << setprecision(14) << "0 0\t" << j << "\t" << result_sc[j] << endl;
			}

			//the vector (J=1, P=1)
			if(result_vx.size() != result_vy.size() || result_vx.size() != result_vz.size()) {
				throw Print_Error_Message("Internal error: Vector correlators are not of equal length");
			}
			for(size_t j = 0; j < result_vx.size(); j++) {
				of << scientific << setprecision(14) << "1 1\t" << j << "\t" << (result_vx[j] + result_vy[j] + result_vz[j]) / 3. << "\t" << result_vx[j] << "\t" << result_vy[j] << "\t" << result_vz[j] << endl;
			}

			//the axial vector (J=1, P=0)
			if(result_ax.size() != result_ay.size() || result_ax.size() != result_az.size()) {
				throw Print_Error_Message("Internal error: Vector correlators are not of equal length");
			}
			for(size_t j = 0; j < result_ax.size(); j++) {
				of << scientific << setprecision(14) << "1 0\t" << j << "\t" << (result_ax[j] + result_ay[j] + result_az[j]) / 3. << "\t" << result_ax[j] << "\t" << result_ay[j] << "\t" << result_az[j] << endl;
			}

			//the avps correlator
			if (parametersInterface.getCorrelatorDirection() == 0)
			  {
				auto result_avps = physics::observables::wilson::calculate_correlator("avps", result, sources, system, interfacesHandler);
				for(size_t j = 0; j < result_avps.size(); j++) {
				  of << scientific << setprecision(14) << "1 0 0 1\t" << j << "\t" << result_avps[j] << endl;
				}
			  }
			of << endl;
		}
	 *
	 *
	 *
	 */





    //throw Print_Error_Message("Function writeCorrelatorToFile not implemented yet!");
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

	/*
	 * auto code = results[i]->get_device()->getCorrelatorCode();
	 *  code->correlator(code->get_correlator_kernel(type), results[i], corr_bufs[i], window_bufs[i]);
	 *
	 *
	 *
	 * Unclear where function correlator is defined since eclipse gives an syntax error msg for the auto code ... line
	 *
	 */








	/* Wilson code: (just for overview)
	 *
	 * static std::vector<hmc_float> calculate_correlator_componentwise(const std::string& type, const std::vector<physics::lattices::Spinorfield*>& corr,
                                                                 const std::vector<physics::lattices::Spinorfield*>& sources, const hardware::System& system,
                                                                 const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& parametersInterface,
                                                                 physics::InterfacesHandler & interfacesHandler)
		{
			// assert single device
			auto first_corr = corr.at(0);
			try_swap_in(first_corr);
			auto first_field_buffers = first_corr->get_buffers();
			const size_t num_buffers = first_field_buffers.size();
			const size_t num_corr_entries = get_num_corr_entries(parametersInterface);

			// for each source
			if(corr.size() != sources.size()) {
				throw std::invalid_argument("Correlated and source fields need to be of the same size.");
			}

			std::vector<const hardware::buffers::Plain<hmc_float>*> results(num_buffers);
			for(size_t i = 0; i < num_buffers; ++i) {
				auto device = first_field_buffers[i]->get_device();
				results[i] = new hardware::buffers::Plain<hmc_float>(num_corr_entries, device);
				results[i]->clear();
			}

			for(size_t i = 0; i < corr.size(); i++) {
				calculate_correlator(type, results, corr.at(i), sources.at(i), system, parametersInterface, interfacesHandler);
			}

			std::vector<hmc_float> host_result(num_corr_entries);
			for(size_t i = 0; i < num_corr_entries; ++i) {
				host_result[i] = 0.;
			}
			for(auto result: results) {
				std::vector<hmc_float> out(num_corr_entries);
				result->dump(out.data());
				for(size_t i = 0; i < num_corr_entries; ++i) {
					logger.trace() << out[i];
					host_result[i] += out[i];
				}
				delete result;
			}
			return host_result;
		}
	 *
	 *
	 *
	 *
	 *
	 */




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

