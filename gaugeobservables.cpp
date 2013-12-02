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

#include "gaugeobservables.h"

#include "meta/util.hpp"
#include "physics/lattices/gaugefield.hpp"

int main(int argc, const char* argv[])
{
	using physics::lattices::Gaugefield;

	try {
		logger.info() << "This executable requires the following parameter value(s) to work properly:";
		logger.info() << "startcondition:\tcontinue";
		meta::Inputparameters parameters(argc, argv);

		//check settings
		if(parameters.get_startcondition() != meta::Inputparameters::start_from_source ) {
			logger.fatal() << "Found wrong startcondition! Aborting..";
			throw Invalid_Parameters("Found wrong startcondition!", "continue", parameters.get_startcondition());
		}

		switchLogLevel(parameters.get_log_level());

		meta::print_info_gaugeobservables(argv[0], parameters);

		ofstream ofile("gaugeobservables.log");
		if(ofile.is_open()) {
			meta::print_info_gaugeobservables(argv[0], &ofile, parameters);
			ofile.close();
		} else {
			logger.warn() << "Could not log file for gaugeobservables.";
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		init_timer.reset();

		const hardware::System system(parameters);
		physics::PRNG prng(system);

		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// gaugeobservables
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();

		const int iter_end = parameters.get_config_read_end();
		const int iter_start = parameters.get_config_read_start();
		const int iter_incr = parameters.get_config_read_incr();
		int iter = 0;

		logger.info() << "Measure gaugeobservables on device(s)... ";

		if(parameters.get_read_multiple_configs()) {
			//main loop
			for(iter = iter_start; iter < iter_end; iter += iter_incr) {
				const std::string config_name = meta::create_configuration_name(parameters, iter);
				logger.info() << "Measure gaugeobservables of configuration: " << config_name;
				const Gaugefield gaugefield(system, prng, config_name);
				if(parameters.get_print_to_screen() ) {
					print_gaugeobservables(gaugefield, iter);
				}
				const std::string gaugeout_name = get_gauge_obs_file_name(parameters, config_name);
				print_gaugeobservables(gaugefield, iter, gaugeout_name);
			}
		} else {
			//in this case only the config from the initialization is taken into account
			const std::string config_name = parameters.get_sourcefile();
			logger.info() << "Measure gaugeobservables of configuration: " << config_name;
			const Gaugefield gaugefield(system, prng, config_name);
			//@todo: adjust the "iter" here to be the number from the sourcefile!!
			if(parameters.get_print_to_screen() ) {
				print_gaugeobservables(gaugefield, iter);
			}
			const std::string gaugeout_name = get_gauge_obs_file_name(parameters, "");
			print_gaugeobservables(gaugefield, iter, gaugeout_name);
		}
		logger.info() << "... done";
		perform_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		total_timer.add();
		general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);

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
