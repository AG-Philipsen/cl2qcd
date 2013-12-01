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

class generalExecutable{
public:
	generalExecutable(int argc, const char* argv[]) : parameters(argc, argv){
		ownName = argv[0];
		totalRuntimeOfExecutable.reset();
		initializationTimer.reset();
		switchLogLevel(parameters.get_log_level());
		system = new hardware::System(parameters);
		initializationTimer.add();
	}
	~generalExecutable(){
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

	void printRuntimeInformationToScreenAndFile()
	{
		printGeneralTimesToScreen();
		printGeneralTimesToFile();
		return;
	}

	void printGeneralTimesToScreen(){
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
	void printGeneralTimesToFile(){
		logger.info() << "## writing general times to file: \"general_time_output\"";
		ofstream ofile;
		ofile.open("general_time_output");
		if(ofile.is_open()) {
			ofile  << "## *******************************************************************" << endl;
			ofile  << "## General Times [mus]:" << endl;
			ofile << "## Total\tInit\tPerformance" << endl;
			ofile  << totalRuntimeOfExecutable.getTime() << "\t" << initializationTimer.getTime() << '\t' << performanceTimer.getTime() << endl;
			ofile.close();
		} else {
			logger.warn() << "Could not open output file for general time output.";
		}
		return;
	}
};

class inverterExecutable : public generalExecutable {
public:
	inverterExecutable(int argc, const char* argv[]) : generalExecutable(argc, argv){
		initializationTimer.reset();
		prng = new physics::PRNG(*system);
		meta::print_info_inverter(ownName, parameters);
		writeInverterLogfile();
		setIterationVariables();
		initializationTimer.add();
	}

	~inverterExecutable(){
		if(parameters.get_profile_solver() ) {
			string profiling_out;
			profiling_out = string(ownName) + string("_profiling_data");
			fstream prof_file;
			prof_file.open(profiling_out.c_str(), std::ios::out | std::ios::app);
			if(prof_file.is_open()) {
				meta::print_info_inverter(ownName, &prof_file, parameters);
				prof_file.close();
			} else {
				logger.warn() << "Could not open " << profiling_out;
			}
			print_solver_profiling(profiling_out);
		}
	}

	void setIterationVariables() {
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

	void saveCurrentPrngStateToFile() {
		logger.info() << "saving current prng state to \"" << filenameForCurrentPrngState << "\"";
		prng->store(filenameForCurrentPrngState);
	}

	void performMeasurementsForSpecificConfiguration() {
		currentConfigurationName = meta::create_configuration_name(parameters,iteration);
		logger.info() << "Measure fermionic observables on configuration: " << currentConfigurationName;
		gaugefield = new physics::lattices::Gaugefield(*system, *prng, currentConfigurationName);
		measureFermionicObservablesOnGaugefield();
	}

	void performMeasurementsForConfigurationGivenInSourcefileParameter() {
		currentConfigurationName = parameters.get_sourcefile();
		logger.info() << "Measure fermionic observables on configuration: " << currentConfigurationName;
		gaugefield = new physics::lattices::Gaugefield(*system, *prng);
		measureFermionicObservablesOnGaugefield();
	}

	void performMeasurementsForSpecificIteration() {
		if (parameters.get_read_multiple_configs()) {
			performMeasurementsForSpecificConfiguration();
		} else {
			performMeasurementsForConfigurationGivenInSourcefileParameter();
		}
		saveCurrentPrngStateToFile();
	}

	void performMeasurements() {
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
	std::string currentConfigurationName;
	int iterationStart;
	int iterationEnd;
	int iterationIncrement;
	int iteration;

	void writeInverterLogfile() {
		ofstream ofile;
		ofile.open("inverter.log");
		if (ofile.is_open()) {
			meta::print_info_inverter(ownName, &ofile, parameters);
			ofile.close();
		} else {
			logger.warn() << "Could not open log file for inverter.";
		}
	}

	void measureFermionicObservablesOnGaugefield()
	{
		using namespace physics;
		using namespace physics::lattices;
		using namespace physics::algorithms;

		auto parameters = system->get_inputparameters();
		if(parameters.get_print_to_screen() ) {
			print_gaugeobservables(*gaugefield, gaugefield->get_parameters_source().trajectorynr_source);
		}
		if(parameters.get_measure_correlators() ) {
			// for the correlator calculation, all sources are needed on the device
			const std::vector<Spinorfield*> sources = create_swappable_sources(*system, *prng, parameters.get_num_sources());
			const std::vector<Spinorfield*> result = create_swappable_spinorfields(*system, sources.size(), parameters.get_place_sources_on_host());

			swap_out(sources);
			swap_out(result);

			perform_inversion(&result, gaugefield, sources, *system);

			logger.info() << "Finished inversion. Starting measurements.";

			swap_in(sources);
			swap_in(result);

			//get name for file to which correlators are to be stored
			std::string corr_fn = meta::get_ferm_obs_corr_file_name(parameters, currentConfigurationName);
			flavour_doublet_correlators(result, sources, corr_fn, *system);
			release_spinorfields(result);
			release_spinorfields(sources);
		}
		if(parameters.get_measure_pbp() ) {
			//get name for file to which pbp is to be stored
			std::string pbp_fn = meta::get_ferm_obs_pbp_file_name(parameters, currentConfigurationName);
			// the chiral condensate needs only one source at a time
			for(int i_sources = 0; i_sources < parameters.get_num_sources(); i_sources ++) {
				auto sources = create_sources(*system, *prng, 1);
				auto result = create_spinorfields(*system, sources.size());

				perform_inversion(&result, gaugefield, sources, *system);

				flavour_doublet_chiral_condensate(result, sources, pbp_fn, gaugefield->get_parameters_source().trajectorynr_source, *system);
				release_spinorfields(result);
				release_spinorfields(sources);
			}
		}
	}
};

int main(int argc, const char* argv[])
{
	try {
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

