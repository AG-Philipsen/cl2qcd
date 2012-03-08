#include "inverter.h"

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
		parameters.print_info_inverter(argv[0]);

		ofstream ofile;
		ofile.open("inverter.log");
		if(ofile.is_open()) {
			parameters.print_info_inverter(argv[0], &ofile);
			ofile.close();
		} else {
			logger.warn() << "Could not open log file for inverter.";
		}

		//get name for file to which correlators are to be stored
		stringstream corr_fn;
		switch ( parameters.get_startcondition() ) {
			case START_FROM_SOURCE :
				corr_fn << parameters.sourcefile << "_correlators.dat" ;
				break;
			case HOT_START :
				corr_fn << "conf.hot_correlators.dat" ;
				break;
			case COLD_START :
				corr_fn << "conf.cold_correlators.dat" ;
				break;
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		init_timer.reset();
		Gaugefield_inverter gaugefield;

		//use 2 devices: one for solver, one for correlator
		int numtasks = 2;
		if(parameters.get_num_dev() != 2 )
			logger.warn() << "Only 1 device demanded by input file. All calculations performed on primary device.";

		cl_device_type primary_device;
		switch ( parameters.get_use_gpu() ) {
			case true :
				primary_device = CL_DEVICE_TYPE_GPU;
				break;
			case false :
				primary_device = CL_DEVICE_TYPE_CPU;
				break;
		}

		//check if correlator-device is a GPU and in that case exit because the kernels are not meant to be executed there
		if ( parameters.get_use_gpu() == false && parameters.get_num_dev() == 2) {
			throw Print_Error_Message("GPU cannot be used for correlator-calculation.", __FILE__, __LINE__);
		}

		logger.trace() << "Init gaugefield" ;
		gaugefield.init(numtasks, primary_device, &parameters);

		logger.info() << "Gaugeobservables:";
		gaugefield.print_gaugeobservables(0);
		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// inverter
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();
		/** @todo usage of solver_timer has to be checked. No output yet */
		usetimer solver_timer;

		logger.info() << "Perform inversion on device.." ;

		gaugefield.create_sources();
		gaugefield.perform_inversion(&solver_timer);

		//flavour_doublet_correlators does a sync at the beginning
		gaugefield.flavour_doublet_correlators(corr_fn.str());

		logger.trace() << "Inversion done" ;
		perform_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		total_timer.add();
		uint64_t totaltime = total_timer.getTime();
		general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);
		//print times from the devices...
		logger.info() << "## Device: Solver";
		(gaugefield.get_task_solver())->print_copy_times(totaltime);
		logger.info() << "## Device: Correlator";
		(gaugefield.get_task_correlator())->print_copy_times(totaltime);

		if(parameters.get_profile_solver() ) {
			stringstream profiling_out;
			profiling_out << argv[0] << "_profiling_data";
			fstream prof_file;
			prof_file.open(profiling_out.str(), std::ios::out | std::ios::app);
			if(prof_file.is_open()) {
				parameters.print_info_inverter(argv[0], &prof_file);
				prof_file.close();
			} else {
				logger.warn() << "Could not open " << profiling_out;
			}
			print_solver_profiling(profiling_out.str());
		}

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
