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

#include "wilsonTwoFlavourCorrelators.hpp"

#include <fstream>
#include <cmath>
#include "../lattices/spinorfield.hpp"

//TODO: refactor. Consider also how to handle interfaces (at the moment both the handler and the correlator interfaces are given around,
//      but from the first one can get the second). One should do something like in the Wilson pbp class in the .cpp file (namely have as
//      private member the interface of the correlator and as additional argument of some methods the interfacesHandler).
//class TwoFlavourCorrelators
//{
//public:
//
//private:
//  const physics::lattices::Gaugefield& gaugefield;
//  const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& parametersInterface;
//  const hardware::System& system;
//  const physics::PRNG& prng;
//  int trajectoryNumber;
//  std::vector<double> correlator;
//  std::ofstream outputToFile;
//  std::string filenameForCorrelatorData;
//  std::string configurationName;
//};

#include "../sources.hpp"
#include "../algorithms/inversion.hpp"

#include <cassert>
#include <fstream>
#include "../../meta/util.hpp"
#include "../lattices/util.hpp"
#include "../lattices/swappable.hpp"
#include "../../hardware/device.hpp"
#include "../../hardware/code/correlator.hpp"
#include "../interfacesHandler.hpp"

static size_t get_num_corr_entries(const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& parametersInterface);
static void calculate_correlator(const std::string& type, const std::vector<const hardware::buffers::Plain<hmc_float>*>& results, physics::lattices::Spinorfield* corr,
                                 physics::lattices::Spinorfield* source, const hardware::System& system,
                                 const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& parametersInterface, physics::InterfacesHandler& interfacesHandler);
static void calculate_correlator(const std::string& type, const std::vector<const hardware::buffers::Plain<hmc_float>*>& results,
                                 physics::lattices::Spinorfield* corr1, physics::lattices::Spinorfield* source1,
                                 physics::lattices::Spinorfield* corr2, physics::lattices::Spinorfield* source2,
                                 physics::lattices::Spinorfield* corr3, physics::lattices::Spinorfield* source3,
                                 physics::lattices::Spinorfield* corr4, physics::lattices::Spinorfield* source4,
                                 const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& params);

static void flavour_doublet_correlators(const std::vector<physics::lattices::Spinorfield*>& result, const std::vector<physics::lattices::Spinorfield*>& sources,
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

static void calculate_correlator(const std::string& type, const std::vector<const hardware::buffers::Plain<hmc_float>*>& results,
                                 physics::lattices::Spinorfield* corr, physics::lattices::Spinorfield* source, const hardware::System& system,
                                 const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& parametersInterface,
                                 physics::InterfacesHandler& interfacesHandler)
{
	try_swap_in(corr);
	try_swap_in(source);

	// assert single device
	auto corr_bufs = corr->get_buffers();
	auto source_bufs = source->get_buffers();

	size_t num_bufs = results.size();
	if(num_bufs != source_bufs.size() || num_bufs != corr_bufs.size()) {
		throw std::invalid_argument("The arguments are using different devices.");
	}

	// the ps_z kernel needs to have the source windowed...
	if(num_bufs > 1 && parametersInterface.getSourceType() != common::sourcetypes::point && type == "ps" && parametersInterface.getCorrelatorDirection() == 3) {
		physics::lattices::Spinorfield window(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		auto window_bufs = window.get_buffers();
		for(size_t i_window = 0; i_window < num_bufs; ++i_window) {
			fill_window(&window, *source, i_window);
			for(size_t i = 0; i < num_bufs; ++i) {
				auto code = results[i]->get_device()->getCorrelatorCode();
				code->correlator(code->get_correlator_kernel(type), results[i], corr_bufs[i], window_bufs[i]);
			}
		}
	} else {
		for(size_t i = 0; i < num_bufs; ++i) {
			auto code = results[i]->get_device()->getCorrelatorCode();
			if(parametersInterface.getSourceType() == common::sourcetypes::point) {
				code->correlator(code->get_correlator_kernel(type), results[i], corr_bufs[i]);
			} else {
				code->correlator(code->get_correlator_kernel(type), results[i], corr_bufs[i], source_bufs[i]);
			}
		}
	}

	try_swap_out(corr);
	try_swap_out(source);
}

static void calculate_correlator(const std::string& type, const std::vector<const hardware::buffers::Plain<hmc_float>*>& results,
                                 physics::lattices::Spinorfield* corr1, physics::lattices::Spinorfield* source1,
                                 physics::lattices::Spinorfield* corr2, physics::lattices::Spinorfield* source2,
                                 physics::lattices::Spinorfield* corr3, physics::lattices::Spinorfield* source3,
                                 physics::lattices::Spinorfield* corr4, physics::lattices::Spinorfield* source4,
                                 const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& parametersInterface)
{
	// assert single device
	try_swap_in(corr1);
	try_swap_in(source1);
	auto corr1_bufs = corr1->get_buffers();
	auto source1_bufs = source1->get_buffers();
	try_swap_in(corr2);
	try_swap_in(source2);
	auto corr2_bufs = corr2->get_buffers();
	auto source2_bufs = source2->get_buffers();
	try_swap_in(corr3);
	try_swap_in(source3);
	auto corr3_bufs = corr3->get_buffers();
	auto source3_bufs = source3->get_buffers();
	try_swap_in(corr4);
	try_swap_in(source4);
	auto corr4_bufs = corr4->get_buffers();
	auto source4_bufs = source4->get_buffers();

	size_t num_bufs = results.size();
	if(num_bufs != source1_bufs.size() || num_bufs != corr1_bufs.size()
	   || num_bufs != source2_bufs.size() || num_bufs != corr2_bufs.size()
	   || num_bufs != source3_bufs.size() || num_bufs != corr3_bufs.size()
	   || num_bufs != source4_bufs.size() || num_bufs != corr4_bufs.size() ) {
		throw std::invalid_argument("The arguments are using different devices.");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
        auto code = results[i]->get_device()->getCorrelatorCode();
		if(parametersInterface.getSourceType() == common::sourcetypes::point) {
			code->correlator(code->get_correlator_kernel(type), results[i], corr1_bufs[i], corr2_bufs[i], corr3_bufs[i], corr4_bufs[i]);
		} else {
			code->correlator(code->get_correlator_kernel(type), results[i], corr1_bufs[i], source1_bufs[i], corr2_bufs[i], source2_bufs[i], corr3_bufs[i], source3_bufs[i], corr4_bufs[i], source4_bufs[i]);
		}
	}

	try_swap_out(corr1);
	try_swap_out(source1);
	try_swap_out(corr2);
	try_swap_out(source2);
	try_swap_out(corr3);
	try_swap_out(source3);
	try_swap_out(corr4);
	try_swap_out(source4);
}

static std::vector<hmc_float> calculate_correlator_componentwise(const std::string& type, const std::vector<physics::lattices::Spinorfield*>& corr,
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

static std::vector<hmc_float> calculate_correlator_colorwise(const std::string& type, const std::vector<physics::lattices::Spinorfield*>& corr,
                                                             const std::vector<physics::lattices::Spinorfield*>& sources,
                                                             const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& parametersInterface)
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

	// TODO adjust correlator kernels!
	for(size_t i = 0; i < corr.size(); i += 4) {
		calculate_correlator(type, results, corr.at(i), sources.at(i), corr.at(i + 1), sources.at(i + 1), corr.at(i + 2), sources.at(i + 2),
		                     corr.at(i + 3), sources.at(i + 3), parametersInterface);
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

static size_t get_num_corr_entries(const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& parametersInterface)
{
    switch (parametersInterface.getCorrelatorDirection()) {
        case 0 :
            return parametersInterface.getNt();
        case 3 :
            return parametersInterface.getNs();
        default :
            std::stringstream errmsg;
            errmsg << "Correlator direction " << parametersInterface.getCorrelatorDirection() << " has not been implemented.";
            throw Print_Error_Message(errmsg.str());
    }
}


std::vector<hmc_float> physics::observables::wilson::calculate_correlator(const std::string& type, const std::vector<physics::lattices::Spinorfield*>& corr,
                                                                          const std::vector<physics::lattices::Spinorfield*>& sources,
                                                                          const hardware::System& system, physics::InterfacesHandler& interfacesHandler)
{
    const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& parametersInterface = interfacesHandler.getWilsonTwoFlavourCorrelatorsParametersInterface();
	if(type == "ps"  || type == "avps" ) {
		return calculate_correlator_componentwise(type, corr, sources, system, parametersInterface, interfacesHandler);
	} else if (type == "sc" || type == "vx" || type == "vy" || type == "vz" || type == "ax" || type == "ay" || type == "az") {
		return calculate_correlator_colorwise(type, corr, sources, parametersInterface);
	} {
		throw Print_Error_Message("Correlator calculation has not been implemented for " + type, __FILE__, __LINE__);
	}
}

void physics::observables::wilson::measureTwoFlavourDoubletCorrelatorsOnGaugefieldAndWriteToFile(const physics::lattices::Gaugefield * gaugefield,
                                                                                                   std::string currentConfigurationName,
                                                                                                   physics::InterfacesHandler & interfacesHandler)
{
    auto system = gaugefield->getSystem();
    const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& parametersInterface = interfacesHandler.getWilsonTwoFlavourCorrelatorsParametersInterface();
	auto prng = gaugefield->getPrng();

	std::string filenameForCorrelatorData = parametersInterface.getCorrelatorFilename(currentConfigurationName);
	// for the correlator calculation, all sources are needed on the device
	const std::vector<physics::lattices::Spinorfield*> sources = physics::create_swappable_sources(*system, *prng, parametersInterface.getNumberOfSources(), interfacesHandler);
	const std::vector<physics::lattices::Spinorfield*> result = physics::lattices::create_swappable_spinorfields(*system, sources.size(), interfacesHandler, parametersInterface.placeSourcesOnHost());
	swap_out(sources);
	swap_out(result);
	physics::algorithms::perform_inversion(&result, gaugefield, sources, *system, interfacesHandler);
	logger.info() << "Finished inversion. Starting measurements.";
	swap_in(sources);
	swap_in(result);
	flavour_doublet_correlators(result, sources, filenameForCorrelatorData, *system, parametersInterface, interfacesHandler);
	release_spinorfields(result);
	release_spinorfields(sources);
}







