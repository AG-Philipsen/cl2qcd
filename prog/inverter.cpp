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

#include "inverter.h"

#include "physics/lattices/gaugefield.hpp"
#include "physics/lattices/spinorfield.hpp"
#include "physics/sources.hpp"
#include "physics/algorithms/flavour_doublet.hpp"
#include "physics/algorithms/inversion.hpp"

#include "meta/util.hpp"

class generalExecutable
{

public:
	generalExecutable(int argc, const char* argv[]) : parameters(argc, argv)
	{
		ownName = argv[0];
		totalRuntimeOfExecutable.reset();
		initializationTimer.reset();
		switchLogLevel(parameters.get_log_level());
		system = new hardware::System(parameters);
		initializationTimer.add();
	}
	~generalExecutable()
	{
		totalRuntimeOfExecutable.add();
		printRuntimeInformationToScreenAndFile();
	}

protected:
	const char* ownName;
	usetimer totalRuntimeOfExecutable;
	usetimer initializationTimer;
	usetimer performanceTimer;
	meta::Inputparameters parameters;
	hardware::System * system;
	ofstream outputToFile;
	const char* generalTimeOutputFilename = "general_time_output";

	void printRuntimeInformationToScreenAndFile()
	{
		printGeneralTimesToScreen();
		printGeneralTimesToFile();
		return;
	}

	void printGeneralTimesToScreen()
	{
		logger.info() << "## *******************************************************************";
		logger.info() << "## General Times [mus]:";
		logger.info() << "## *******************************************************************";
		logger.info() << "## Program Parts:\t" << setfill(' ') << setw(5) << "total" << '\t' << setw(5) << "perc";
		logger.info() << "## Total:\t" << setfill(' ') << setw(12) << totalRuntimeOfExecutable.getTime();
		logger.info() << "## Init.:\t" << setfill(' ') << setw(12) << initializationTimer.getTime() << '\t' << fixed << setw(5) << setprecision(1) << percent(initializationTimer.getTime(), totalRuntimeOfExecutable.getTime()) ;
		logger.info() << "## Perf.:\t" << setfill(' ') << setw(12) << performanceTimer.getTime() << '\t' << fixed << setw(5) << setprecision(1) << percent(performanceTimer.getTime(), totalRuntimeOfExecutable.getTime()) ;
		logger.info() << "## *******************************************************************";
		return;
	}
	void printGeneralTimesToFile()
	{
		logger.info() << "## writing general times to file: \"" << generalTimeOutputFilename << "\"";
		outputToFile.open(generalTimeOutputFilename);
		if(outputToFile.is_open()) {
			outputToFile  << "## *******************************************************************" << endl;
			outputToFile  << "## General Times [mus]:" << endl;
			outputToFile << "## Total\tInit\tPerformance" << endl;
			outputToFile  << totalRuntimeOfExecutable.getTime() << "\t" << initializationTimer.getTime() << '\t' << performanceTimer.getTime() << endl;
			outputToFile.close();
		} else {
			logger.warn() << "Could not open output file for general time output.";
		}
		return;
	}
};

class inverterExecutable : public generalExecutable
{
public:
	inverterExecutable(int argc, const char* argv[]) : generalExecutable(argc, argv)
	{
		initializationTimer.reset();
		prng = new physics::PRNG(*system);
		meta::print_info_inverter(ownName, parameters);
		writeInverterLogfile();
		setIterationVariables();
		initializationTimer.add();
	}

	~inverterExecutable()
	{
		if (parameters.get_profile_solver()) {
			writeProfilingDataToFile();
		}
	}

	void performMeasurements()
	{
		performanceTimer.reset();
		logger.trace() << "Perform inversion(s) on device..";
		for (iteration = iterationStart; iteration < iterationEnd; iteration += iterationIncrement)
		{
			performMeasurementsForSpecificIteration();
		}
		logger.trace() << "Inversion(s) done";
		performanceTimer.add();
	}

protected:
	physics::PRNG * prng;
	physics::lattices::Gaugefield * gaugefield;
	const std::string filenameForCurrentPrngState = "prng.inverter.save";
	const std::string filenameForInverterLogfile = "inverter.log";
	const std::string filenameForProfilingData = string(ownName) + string("_profiling_data");
	std::string filenameForTwoFlavourDoubletChiralCondensateData;
	std::string filenameForTwoFlavourDoubletCorrelatorData;
	std::string currentConfigurationName;
	std::fstream outputStreamForProfilingData;
	int iterationStart;
	int iterationEnd;
	int iterationIncrement;
	int iteration;

	void setIterationVariables()
	{
		iterationStart =
				(parameters.get_read_multiple_configs()) ?
						parameters.get_config_read_start() : 0;
		iterationEnd =
				(parameters.get_read_multiple_configs()) ?
						parameters.get_config_read_end() + 1 : 1;
		iterationIncrement =
				(parameters.get_read_multiple_configs()) ?
						parameters.get_config_read_incr() : 1;
	}

	void saveCurrentPrngStateToFile()
	{
		logger.info() << "saving current prng state to \"" << filenameForCurrentPrngState << "\"";
		prng->store(filenameForCurrentPrngState);
	}

	void initializeGaugefieldAccordingToIterationVariable()
	{
		currentConfigurationName = meta::create_configuration_name(parameters,iteration);
		gaugefield = new physics::lattices::Gaugefield(*system, *prng, currentConfigurationName);
	}

	void initializeGaugefieldAccordingToConfigurationGivenInSourcefileParameter()
	{
		currentConfigurationName = parameters.get_sourcefile();
		gaugefield = new physics::lattices::Gaugefield(*system, *prng);
	}

	void initializeGaugefield()
	{
		if (parameters.get_read_multiple_configs()) {
			initializeGaugefieldAccordingToIterationVariable();
		} else {
			initializeGaugefieldAccordingToConfigurationGivenInSourcefileParameter();
		}
	}

	void performMeasurementsForSpecificIteration()
	{
		initializeGaugefield();
		measureFermionicObservablesOnGaugefield();
		saveCurrentPrngStateToFile();
	}

	void writeInverterLogfile()
	{
		outputToFile.open(filenameForInverterLogfile);
		if (outputToFile.is_open()) {
			meta::print_info_inverter(ownName, &outputToFile, parameters);
			outputToFile.close();
		} else {
			logger.warn() << "Could not open log file for inverter.";
		}
	}

	void writeProfilingDataToFile()
	{
		outputStreamForProfilingData.open(filenameForProfilingData.c_str(),	std::ios::out | std::ios::app);
		if (outputStreamForProfilingData.is_open()) {
			meta::print_info_inverter(ownName, &outputStreamForProfilingData, parameters);
			outputStreamForProfilingData.close();
		} else {
			logger.warn() << "Could not open " << filenameForProfilingData;
		}
		print_solver_profiling(filenameForProfilingData);
	}

	void measureTwoFlavourDoubletCorrelatorsOnGaugefield()
	{
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

	void measureTwoFlavourDoubletChiralCondensateOnGaugefield()
	{
		filenameForTwoFlavourDoubletChiralCondensateData = meta::get_ferm_obs_pbp_file_name(parameters,	currentConfigurationName);
		int sourceNumber = 0;

		for (; sourceNumber < parameters.get_num_sources(); sourceNumber++)
		{
			auto sources = physics::create_sources(*system, *prng, 1);
			auto result = physics::lattices::create_spinorfields(*system, sources.size());
			physics::algorithms::perform_inversion(&result, gaugefield, sources, *system);
			physics::algorithms::flavour_doublet_chiral_condensate(result, sources, filenameForTwoFlavourDoubletChiralCondensateData, gaugefield->get_parameters_source().trajectorynr_source, *system);
			release_spinorfields(result);
			release_spinorfields(sources);
		}
	}

	void measureFermionicObservablesOnGaugefield() {
		logger.info() << "Measure fermionic observables on configuration: "	<< currentConfigurationName;
		if (parameters.get_print_to_screen()) {
			print_gaugeobservables(*gaugefield,	gaugefield->get_parameters_source().trajectorynr_source);
		}
		if (parameters.get_measure_correlators()) {
			measureTwoFlavourDoubletCorrelatorsOnGaugefield();
		}
		if (parameters.get_measure_pbp()) {
			measureTwoFlavourDoubletChiralCondensateOnGaugefield();
		}
	}
};

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

