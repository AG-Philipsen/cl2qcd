#include "transportcoefficient.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char* argv[])
{
	try {
		po::options_description desc("Allowed options");
		desc.add_options()
		("help,h", "Produce this help message")
		("input-file", po::value<std::string>(), "File containing the input parameters")
		("log-level", po::value<std::string>(), "Minimum output log level: ALL TRACE DEBUG INFO WARN ERROR FATAL OFF");
		po::positional_options_description pos_opts;
		pos_opts.add("input-file", 1);
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_opts).run(), vm);
		if( vm.count( "help" ) ) { // see http://stackoverflow.com/questions/5395503/required-and-optional-arguments-using-boost-library-program-options as to why this is done before po::notifiy(vm)
			std::cout << desc << '\n';
			return 0;
		}
		po::notify(vm); // checks whether all required arguments are set

		if(vm.count("log-level")) {
			switchLogLevel(vm["log-level"].as<std::string>());
		}

		if(!vm.count("input-file")) {
			logger.fatal() << "No input file specified. Please specify a file containing the input parameters.";
		}

		const char* inputfile = vm["input-file"].as<std::string>().c_str();
		inputparameters parameters;
		parameters.readfile(inputfile);
		parameters.print_info_tkkappa(argv[0]);

		//name of file to store gauge observables
		stringstream gaugeout_name;
		gaugeout_name << "gaugeobservables_beta" << parameters.get_beta();

		fstream logfile;
		logfile.open("tk_kappa_hybrid.log", std::ios::out | std::ios::app);
		if(logfile.is_open()) {
			parameters.print_info_tkkappa(argv[0], &logfile);
			logfile.close();
		} else {
			throw File_Exception("tk_kappa_hybrid.log");
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		init_timer.reset();

		Gaugefield_heatbath_kappa gaugefield;
		int numtasks = 2;

		// this is the device type for the heatbath
		cl_device_type primary_device_type = CL_DEVICE_TYPE_GPU;

		gaugefield.init(numtasks, primary_device_type, &parameters);

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
		gaugefield.print_gaugeobservables(iter);
		gaugefield.print_gaugeobservables(iter, gaugeout_name.str());
		iter++;

		//first iteration: whether we want to do auto-timing
		int nheat_frequency = parameters.get_writefrequency();
		if(parameters.get_use_autotuning() == true) {
			gaugefield.perform_tasks(parameters.get_writefrequency(), parameters.get_overrelaxsteps(), &nheat_frequency);
		} else {
			gaugefield.perform_tasks(nheat_frequency, parameters.get_overrelaxsteps());
		}
		gaugefield.synchronize(0);
		gaugefield.print_gaugeobservables(iter);
		//    gaugefield.print_gaugeobservables_from_task(iter,0);
		//    gaugefield.print_gaugeobservables_from_task(iter,1);
		gaugefield.print_gaugeobservables(iter, gaugeout_name.str());
		gaugefield.print_kappa(iter, "kappa_clover.dat");

		for(iter = 2; iter < parameters.get_heatbathsteps() / nheat_frequency; iter++) {
			gaugefield.perform_tasks(parameters.get_writefrequency(), parameters.get_overrelaxsteps());
			gaugefield.synchronize(0);
			gaugefield.print_gaugeobservables(iter);
			//    gaugefield.print_gaugeobservables_from_task(iter,0);
			//    gaugefield.print_gaugeobservables_from_task(iter,1);
			gaugefield.print_gaugeobservables(iter, gaugeout_name.str());
			gaugefield.print_kappa(iter, "kappa_clover.dat");
		}

		gaugefield.save("conf.save");
		logger.trace() << "... done";
		perform_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		total_timer.add();
		uint64_t totaltime = total_timer.getTime();
		general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);
		//print times from the devices...
		logger.info() << "## Device: Heatbath";
		(gaugefield.get_task_heatbath())->print_copy_times(totaltime);
		logger.info() << "## Device: Kappa";
		(gaugefield.get_task_kappa())->print_copy_times(totaltime);

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
	}

	return 0;

}
