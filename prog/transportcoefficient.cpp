#include "general_header.h"

#include "meta/util.hpp"
#include "physics/prng.hpp"
#include "physics/lattices/gaugefield.hpp"
#include "physics/algorithms/heatbath.hpp"
#include "physics/algorithms/kappa_clover.hpp"
#include <fstream>

static void print_kappa(hmc_float kappa, int iter, std::string filename);

int main(int argc, const char* argv[])
{
	using physics::PRNG;
	using physics::lattices::Gaugefield;
	using physics::algorithms::heatbath;
	using physics::algorithms::kappa_clover;

	try {
		meta::Inputparameters parameters(argc, argv);
		switchLogLevel(parameters.get_log_level());

		meta::print_info_tkkappa(argv[0], parameters);

		fstream logfile;
		logfile.open("tk_kappa_hybrid.log", std::ios::out | std::ios::app);
		if(logfile.is_open()) {
			meta::print_info_tkkappa(argv[0], &logfile, parameters);
			logfile.close();
		} else {
			throw File_Exception("tk_kappa_hybrid.log");
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		init_timer.reset();

		hardware::System system(parameters);
		PRNG prng(system);
		Gaugefield gaugefield(system, prng);

		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Do the iterations
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();

		logger.trace() << "Start thermalization" ;
		int ntherm = parameters.get_thermalizationsteps();
		for(int i = 0; i < ntherm; i++) {
			heatbath(gaugefield, prng);
		}

		logger.info() << "Start hybrid heatbath and tk_kappa";
		//first output is considered to be zeroth iteration
		int iter = 0;
		std::string gaugeout_name = get_gauge_obs_file_name(parameters, "");
		print_gaugeobservables(gaugefield, iter);
		print_gaugeobservables(gaugefield, iter, gaugeout_name);
		iter++;

		//first iteration: whether we want to do auto-timing
		int nheat_frequency = parameters.get_writefrequency();
//		if(parameters.get_use_autotuning() == true) {
//			gaugefield.perform_tasks(parameters.get_writefrequency(), parameters.get_overrelaxsteps(), &nheat_frequency);
//		} else {
		for(int iter = 0; iter < nheat_frequency; iter++) {
			heatbath(gaugefield, prng, parameters.get_overrelaxsteps());
		}
		hmc_float kappa = kappa_clover(gaugefield, parameters.get_beta());
//		}
		print_gaugeobservables(gaugefield, iter);
		print_gaugeobservables(gaugefield, iter, gaugeout_name);
		print_kappa(kappa, iter, "kappa_clover.dat");

		for(iter = 2; iter < parameters.get_heatbathsteps() / nheat_frequency; iter++) {
			for(int iter = 0; iter < nheat_frequency; iter++) {
				heatbath(gaugefield, prng, parameters.get_overrelaxsteps());
			}
			kappa = kappa_clover(gaugefield, parameters.get_beta());
			print_gaugeobservables(gaugefield, iter);
			print_gaugeobservables(gaugefield, iter, gaugeout_name);
			print_kappa(kappa, iter, "kappa_clover.dat");
		}

		gaugefield.save("conf.save");
		logger.trace() << "... done";
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

static void print_kappa(hmc_float kappa, int iter, std::string filename)
{
	std::ofstream outfile(filename.c_str(), std::ios::app);
	if(!outfile.is_open()) throw File_Exception(filename);
	outfile.width(8);
	outfile << iter;
	outfile << "\t";
	outfile.precision(15);
	outfile << kappa << std::endl;
	outfile.close();
	return;
}
