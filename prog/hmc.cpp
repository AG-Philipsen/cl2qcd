#include "hmc.h"

#include "meta/util.hpp"

int main(int argc, const char* argv[])
{
	try {
		meta::Inputparameters parameters(argc, argv);
		switchLogLevel(parameters.get_log_level());

		meta::print_info_hmc(argv[0], parameters);

		ofstream ofile;
		ofile.open("hmc.log");
		if(ofile.is_open()) {
			meta::print_info_hmc(argv[0], &ofile, parameters);
			ofile.close();
		} else {
			logger.warn() << "Could not log file for hmc.";
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		init_timer.reset();
		hmc_observables obs;

		hardware::System system(parameters);
		physics::PRNG prng(system);
		Gaugefield_hmc gaugefield(&system);

		//use 1 task: the hmc-algorithm
		int numtasks = 1;
		if(parameters.get_device_count() == 2 )
			logger.warn() << "Only 1 device demanded by input file. All calculations performed on primary device.";

		cl_device_type primary_device = parameters.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU;

		logger.trace() << "Init gaugefield" ;
		gaugefield.init(numtasks, primary_device, prng);


		logger.info() << "Gaugeobservables:";
		gaugefield.print_gaugeobservables(0);
		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// hmc
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();
		/** @todo usage of solver_timer has to be checked. No output yet */
		usetimer solver_timer;

		//start from the iterationnumber from sourcefile
		//NOTE: this is 0 in case of cold or hot start
		int iter = gaugefield.get_parameters_source().trajectorynr_source;
		int hmc_iter = iter + parameters.get_hmcsteps();
		hmc_float acc_rate = 0.;
		int writefreq = parameters.get_writefrequency();
		int savefreq = parameters.get_savefrequency();

		logger.info() << "perform HMC on device(s)... ";

		//main hmc-loop
		for(; iter < hmc_iter; iter ++) {
			//generate new random-number for Metropolis step
			hmc_float rnd_number = prng.get_double();
			gaugefield.perform_hmc_step(&obs, iter, rnd_number, &solver_timer, prng);
			acc_rate += obs.accept;
			if( ( (iter + 1) % writefreq ) == 0 ) {
				std::string gaugeout_name = meta::get_hmc_obs_file_name(parameters, "");
				gaugefield.print_hmcobservables(obs, iter, gaugeout_name);
			} else if(parameters.get_print_to_screen() )
				gaugefield.print_hmcobservables(obs, iter);

			if( parameters.get_saveconfigs() == true && ( (iter + 1) % savefreq ) == 0 ) {
				// save gaugefield
				gaugefield.synchronize(0);
				gaugefield.save(iter + 1);

				// save prng
				std::ostringstream prng_file_name;
				prng_file_name << "prng.";
				prng_file_name.fill('0');
				prng_file_name.width(parameters.get_config_number_digits());
				prng_file_name << iter + 1;
				prng.store(prng_file_name.str());
			}
		}

		//always save the config on the last iteration
		gaugefield.synchronize(0);
		std::string outputfile = "conf.save";
		logger.info() << "saving current gaugefield to file \"" << outputfile << "\"";
		gaugefield.save(outputfile, iter + 1);
		outputfile = "prng.save";
		logger.info() << "saving current prng state to \"" << outputfile << "\"";
		prng.store(outputfile);

		logger.info() << "HMC done";
		logger.info() << "Acceptance rate: " << fixed <<  setprecision(1) << percent(acc_rate, parameters.get_hmcsteps()) << "%";
		perform_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		total_timer.add();
		general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// free variables
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		gaugefield.finalize();

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
