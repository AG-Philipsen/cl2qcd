/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#include "inverterExecutable.h"

inline inverterExecutable::inverterExecutable(int argc, const char* argv[]) : measurementExecutable(argc, argv)
{
	initializationTimer.reset();
	printParametersToScreenAndFile();
	initializationTimer.add();
}

inline void inverterExecutable::writeInverterLogfile()
{
	outputToFile.open(filenameForInverterLogfile);
	if (outputToFile.is_open()) {
		meta::print_info_inverter(ownName, &outputToFile, parameters);
		outputToFile.close();
	} else {
		throw File_Exception(filenameForInverterLogfile);
	}
}

void inverterExecutable::printParametersToScreenAndFile()
{
	meta::print_info_inverter(ownName, parameters);
	writeInverterLogfile();
}

inline void inverterExecutable::measureTwoFlavourDoubletCorrelatorsOnGaugefield() {
	filenameForTwoFlavourDoubletCorrelatorData = meta::get_ferm_obs_corr_file_name(parameters, currentConfigurationName);
	// for the correlator calculation, all sources are needed on the device
	const std::vector<physics::lattices::Spinorfield*> sources = physics::create_swappable_sources(*system, *prng, parameters.get_num_sources());
	const std::vector<physics::lattices::Spinorfield*> result = physics::lattices::create_swappable_spinorfields(*system, sources.size(), parameters.get_place_sources_on_host());
	swap_out(sources);
	swap_out(result);
	physics::algorithms::perform_inversion(&result, gaugefield, sources, *system);
	logger.info() << "Finished inversion. Starting measurements.";
	swap_in(sources);
	swap_in(result);
	physics::algorithms::flavour_doublet_correlators(result, sources, filenameForTwoFlavourDoubletCorrelatorData, *system);
	release_spinorfields(result);
	release_spinorfields(sources);
}

inline void inverterExecutable::measureTwoFlavourDoubletChiralCondensateOnGaugefield() {
	filenameForTwoFlavourDoubletChiralCondensateData = meta::get_ferm_obs_pbp_file_name(parameters, currentConfigurationName);
	int sourceNumber = 0;
	for (; sourceNumber < parameters.get_num_sources(); sourceNumber++) {
		auto sources = physics::create_sources(*system, *prng, 1);
		auto result = physics::lattices::create_spinorfields(*system, sources.size());
		physics::algorithms::perform_inversion(&result, gaugefield, sources, *system);
		physics::algorithms::flavour_doublet_chiral_condensate(result, sources, filenameForTwoFlavourDoubletChiralCondensateData, gaugefield->get_parameters_source().trajectorynr_source, *system);
		release_spinorfields(result);
		release_spinorfields(sources);
	}
}

inline void inverterExecutable::performApplicationSpecificMeasurements() {
	logger.info() << "Measure fermionic observables on configuration: " << currentConfigurationName;
	if (parameters.get_print_to_screen()) {
		print_gaugeobservables(*gaugefield, gaugefield->get_parameters_source().trajectorynr_source);
	}
	if (parameters.get_measure_correlators()) {
		measureTwoFlavourDoubletCorrelatorsOnGaugefield();
	}
	if (parameters.get_measure_pbp()) {
		measureTwoFlavourDoubletChiralCondensateOnGaugefield();
	}
}

int main(int argc, const char* argv[])
{
	try{
		inverterExecutable inverterInstance(argc, argv);
		inverterInstance.performMeasurements();
	} //try
	//exceptions from Opencl classes
	catch (Opencl_Error& e) {
		logger.fatal() << e.what();
		exit(1);
	} catch (File_Exception& fe) {
		logger.fatal() << "Could not open file: " << fe.get_filename();
		logger.fatal() << "Aborting.";
		exit(1);
	} catch (Print_Error_Message& em) {
		logger.fatal() << em.what();
		exit(1);
	} catch (Invalid_Parameters& es) {
		logger.fatal() << es.what();
		exit(1);
	}

	return 0;
}

