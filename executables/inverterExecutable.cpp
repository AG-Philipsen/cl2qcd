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

#include "../physics/observables/wilsonTwoFlavourChiralCondensate.hpp"
inverterExecutable::inverterExecutable(int argc, const char* argv[]) : measurementExecutable(argc, argv)
{
	initializationTimer.reset();
	printParametersToScreenAndFile();
	initializationTimer.add();
}

void inverterExecutable::writeInverterLogfile()
{
  outputToFile.open(filenameForLogfile,
		    std::ios::out | std::ios::app);
	if (outputToFile.is_open()) {
		meta::print_info_inverter(&outputToFile, parameters);
		outputToFile.close();
	} else {
		throw File_Exception(filenameForLogfile);
	}
}

void inverterExecutable::printParametersToScreenAndFile()
{
	meta::print_info_inverter(parameters);
	writeInverterLogfile();
}

void inverterExecutable::measureTwoFlavourDoubletCorrelatorsOnGaugefield() {
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

void inverterExecutable::performApplicationSpecificMeasurements() {
	logger.info() << "Measure fermionic observables on configuration: " << currentConfigurationName;
	gaugeObservablesInstance.measureGaugeObservables(gaugefield, gaugefield->get_parameters_source().trajectorynr_source);
	if (parameters.get_measure_correlators()) {
		measureTwoFlavourDoubletCorrelatorsOnGaugefield();
	}
	if (parameters.get_measure_pbp()) {
	  physics::observables::wilson:: measureChiralCondensateAndWriteToFile(gaugefield, currentConfigurationName);
	}
}

