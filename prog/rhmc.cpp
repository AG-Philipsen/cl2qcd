/** @file
 * rhmc algorithm main()
 * 
 * (c) 2013 Alessandro Sciarra <sciarra@compeng.uni-frankfurt.de>
 */

#include "rhmc.h"

#include "meta/util.hpp"
#include <cmath>

static void print_hmcobservables(const hmc_observables& obs, int iter, const std::string& filename, const meta::Inputparameters& params);
static void print_hmcobservables(const hmc_observables& obs, int iter);

int main(int argc, const char* argv[])
{
	using physics::algorithms::perform_rhmc_step;

	try {
		meta::Inputparameters parameters(argc, argv);
		switchLogLevel(parameters.get_log_level());
		
		meta::print_info_rhmc(argv[0], parameters);

		ofstream ofile;
		ofile.open("rhmc.log");
		if(ofile.is_open()) {
			meta::print_info_hmc(argv[0], &ofile, parameters);
			ofile.close();
		} else {
			logger.warn() << "Could not log file for rhmc.";
		}
#if 0
		/////////////////////////////////////////////////////////////////////////////////
		// Initialization
		/////////////////////////////////////////////////////////////////////////////////
		init_timer.reset();
		hmc_observables obs;

		hardware::System system(parameters);
		physics::PRNG prng(system);

		logger.trace() << "Init gaugefield" ;
		physics::lattices::Gaugefield gaugefield(system, prng);

		logger.info() << "Gaugeobservables:";
		print_gaugeobservables(gaugefield, 0);
		init_timer.add();

		/////////////////////////////////////////////////////////////////////////////////
		// RHMC
		/////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();

		//start from the iterationnumber from sourcefile
		//NOTE: this is 0 in case of cold or hot start
		int iter = gaugefield.get_parameters_source().trajectorynr_source;
		const int hmc_iter = iter + parameters.get_hmcsteps();
		hmc_float acc_rate = 0.;
		const int writefreq = parameters.get_writefrequency();
		const int savefreq = parameters.get_savefrequency();

		logger.info() << "perform HMC on device(s)... ";

		//main hmc-loop
		for(; iter < hmc_iter; iter ++) {
			//generate new random-number for Metropolis step
			const hmc_float rnd_number = prng.get_double();

			obs = perform_hmc_step(&gaugefield, iter, rnd_number, prng, system);

			acc_rate += obs.accept;
			if( ( (iter + 1) % writefreq ) == 0 ) {
				std::string gaugeout_name = meta::get_hmc_obs_file_name(parameters, "");
				print_hmcobservables(obs, iter, gaugeout_name, parameters);
			} else if(parameters.get_print_to_screen() ) {
				print_hmcobservables(obs, iter);
			}

			if( savefreq != 0 && ( (iter + 1) % savefreq ) == 0 ) {
				// save gaugefield
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
		std::string outputfile = "conf.save";
		logger.info() << "saving current gaugefield to file \"" << outputfile << "\"";
		gaugefield.save(outputfile, iter);
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
#endif
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





static void print_hmcobservables(const hmc_observables& obs, const int iter, const std::string& filename, const meta::Inputparameters& params)
{
	const hmc_float exp_deltaH = std::exp(obs.deltaH);
	std::fstream hmcout(filename.c_str(), std::ios::out | std::ios::app);
	if(!hmcout.is_open()) throw File_Exception(filename);
	hmcout << iter << "\t";
	hmcout.width(8);
	hmcout.precision(15);
	//print plaquette (plaq, tplaq, splaq)
	hmcout << obs.plaq << "\t" << obs.tplaq << "\t" << obs.splaq;
	//print polyakov loop (re, im, abs)
	hmcout << "\t" << obs.poly.re << "\t" << obs.poly.im << "\t" << sqrt(obs.poly.re * obs.poly.re + obs.poly.im * obs.poly.im);
	//print deltaH, exp(deltaH), acceptance-propability, accept (yes or no)
	hmcout <<  "\t" << obs.deltaH << "\t" << exp_deltaH << "\t" << obs.prob << "\t" << obs.accept;
	//print number of iterations used in inversions with full and force precision
	/**
	 * @todo: The counters should be implemented once the solver class is used!"
	 * until then, only write "0"!
	 */
	int iter0 = 0;
	int iter1 = 0;
	hmcout << "\t" << iter0 << "\t" << iter1;
	if(params.get_use_mp() ) {
		hmcout << "\t" << iter0 << "\t" << iter1;
	}
	if(meta::get_use_rectangles(params) ) {
		//print rectangle value
		hmcout << "\t" << obs.rectangles;
	}
	hmcout << std::endl;
	hmcout.close();

	//print to screen
	print_hmcobservables(obs, iter);
}

static void print_hmcobservables(const hmc_observables& obs, const int iter)
{
	using namespace std;
	//short version of output, all obs are collected in the output file anyways...
	logger.info() << "\tHMC [OBS]:\t" << iter << setw(8) << setfill(' ') << "\t" << setprecision(15) << obs.plaq << "\t" << obs.poly.re << "\t" << obs.poly.im;
}