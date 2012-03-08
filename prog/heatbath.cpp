#include "heatbath.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char* argv[])
{
	try {
		po::options_description desc("Allowed options");
		desc.add_options()
		("help,h", "Produce this help message")
		("input-file", po::value<std::string>()->required(), "File containing the input parameters")
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

		const char* inputfile = vm["input-file"].as<std::string>().c_str();
		inputparameters parameters;
		parameters.readfile(inputfile);
		parameters.print_info_heatbath(argv[0]);

		//name of file to store gauge observables
		stringstream gaugeout_name;
		gaugeout_name << "gaugeobservables_beta" << parameters.get_beta();

		fstream logfile;
		logfile.open("heatbath.log", std::ios::out | std::ios::app);
		if(logfile.is_open()) {
			parameters.print_info_heatbath(argv[0], &logfile);
			logfile.close();
		} else {
			throw File_Exception("heatbath.log");
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		init_timer.reset();

		Gaugefield_heatbath gaugefield;

		cl_device_type primary_device_type;
		//check whether GPU should be used
		if(parameters.get_use_gpu() == true) {
			primary_device_type = CL_DEVICE_TYPE_GPU;
		} else {
			primary_device_type = CL_DEVICE_TYPE_CPU;
		}
		gaugefield.init(1, primary_device_type, &parameters);
		logger.trace() << "initial gaugeobservables: ";
		gaugefield.print_gaugeobservables(0);

		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Heatbath
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();

		logger.trace() << "Start thermalization" ;
		int ntherm = parameters.get_thermalizationsteps();
		for(int iter = 0; iter < ntherm; iter++) gaugefield.perform_tasks(0);

		int nsteps = parameters.get_heatbathsteps();
		int overrelaxsteps = parameters.get_overrelaxsteps();
		int writefreq = parameters.get_writefrequency();
		int savefreq = parameters.get_savefrequency();

		logger.info() << "Start heatbath";

		for(int i = 0; i < nsteps; i++) {
			gaugefield.perform_tasks(overrelaxsteps);
			if( ( (i + 1) % writefreq ) == 0 ) {
				gaugefield.synchronize(0);
				gaugefield.print_gaugeobservables_from_task(i, 0, gaugeout_name.str());
			}
			if( parameters.get_saveconfigs() == true && ( (i + 1) % savefreq ) == 0 ) {
				gaugefield.synchronize(0);
				gaugefield.save(i);
			}
		}

		gaugefield.synchronize(0);
		gaugefield.save(nsteps);
		logger.trace() << "heatbath done";
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
