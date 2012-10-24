#include "transportcoefficient.h"

#include "meta/util.hpp"

int main(int argc, const char* argv[])
{
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
		Gaugefield_heatbath_kappa gaugefield(&system);
		int numtasks = 2;

		// this is the device type for the heatbath
		cl_device_type primary_device = parameters.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU;

		gaugefield.init(numtasks, primary_device);

		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Do the iterations
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();

		logger.trace() << "Start thermalization" ;
		int ntherm = parameters.get_thermalizationsteps();
		if(ntherm > 0) gaugefield.perform_heatbath(ntherm, 0);

		logger.info() << "Start hybrid heatbath and tk_kappa";
		//first output is considered to be zeroth iteration
		int iter = 0;
		std::string gaugeout_name = get_gauge_obs_file_name(parameters, "");
		gaugefield.print_gaugeobservables(iter);
		gaugefield.print_gaugeobservables(iter, gaugeout_name);
		iter++;

		//first iteration: whether we want to do auto-timing
		int nheat_frequency = parameters.get_writefrequency();
//		if(parameters.get_use_autotuning() == true) {
//			gaugefield.perform_tasks(parameters.get_writefrequency(), parameters.get_overrelaxsteps(), &nheat_frequency);
//		} else {
		gaugefield.perform_tasks(nheat_frequency, parameters.get_overrelaxsteps());
//		}
		gaugefield.synchronize(0);
		gaugefield.print_gaugeobservables(iter);
		//    gaugefield.print_gaugeobservables_from_task(iter,0);
		//    gaugefield.print_gaugeobservables_from_task(iter,1);
		gaugefield.print_gaugeobservables(iter, gaugeout_name);
		gaugefield.print_kappa(iter, "kappa_clover.dat");

		for(iter = 2; iter < parameters.get_heatbathsteps() / nheat_frequency; iter++) {
			gaugefield.perform_tasks(parameters.get_writefrequency(), parameters.get_overrelaxsteps());
			gaugefield.synchronize(0);
			gaugefield.print_gaugeobservables(iter);
			//    gaugefield.print_gaugeobservables_from_task(iter,0);
			//    gaugefield.print_gaugeobservables_from_task(iter,1);
			gaugefield.print_gaugeobservables(iter, gaugeout_name);
			gaugefield.print_kappa(iter, "kappa_clover.dat");
		}

		gaugefield.save("conf.save");
		logger.trace() << "... done";
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
