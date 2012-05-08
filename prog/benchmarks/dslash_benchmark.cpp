#include "../inverter.h"

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
		parameters.print_info_inverter(argv[0]);

		ofstream ofile;
		ofile.open("inverter.log");
		if(ofile.is_open()) {
			parameters.print_info_inverter(argv[0], &ofile);
			ofile.close();
		} else {
			logger.warn() << "Could not log file for inverter.";
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

		if(parameters.get_profile_solver() == false) {
			logger.warn() << "solver times will not be measured!";
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		init_timer.reset();
		Gaugefield_inverter gaugefield;

		//one needs 2 task here
		int numtasks = 2;
		if(parameters.get_num_dev() != numtasks )
			logger.warn() << "Exactly 2 devices demanded by benchmark executable. All calculations performed on primary device.";

		cl_device_type primary_device;
		switch ( parameters.get_use_gpu() ) {
			case true :
				primary_device = CL_DEVICE_TYPE_GPU;
				break;
			case false :
				primary_device = CL_DEVICE_TYPE_CPU;
				break;
		}

		logger.trace() << "Init gaugefield" ;
		gaugefield.init(numtasks, primary_device, &parameters);

		logger.info() << "Gaugeobservables:";
		gaugefield.print_gaugeobservables(0);

		//init needed buffers again
		//these are: 2 eoprec spinorfield, 1 gaugefield
		size_t eoprec_spinorfield_buffer_size = gaugefield.get_task_solver()->get_eoprec_spinorfield_buffer_size();
		cl_mem sf1, sf2, gf;

		gf = gaugefield.get_task_solver()->create_rw_buffer(gaugefield.get_task_solver()->getGaugefieldBufferSize());
		sf1 = gaugefield.get_task_solver()->create_rw_buffer(eoprec_spinorfield_buffer_size);
		sf2 = gaugefield.get_task_solver()->create_rw_buffer(eoprec_spinorfield_buffer_size);

		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// dslash-benchmark
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();

		int hmc_iter = parameters.get_hmcsteps();
		int iter;

		logger.trace() << "Perform " << hmc_iter << "of dslash benchmarking (EVEN + ODD) for each step";
		for(iter = 0; iter < hmc_iter; iter ++) {
			gaugefield.get_task_solver()->dslash_eo_device(sf1, sf2, gf, EVEN);
			gaugefield.get_task_solver()->dslash_eo_device(sf1, sf2, gf, ODD);
		}
		logger.trace() << "dslash benchmarking done" ;
		perform_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//release mem-obj
		int clerr;
		clerr = clReleaseMemObject(sf1);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(sf2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(gf);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);

		//CP: for the moment, init the buffers again in order that no segfault happens when freeing the already freed buffers
		gaugefield.get_task_solver()->fill_buffers();

		total_timer.add();
		uint64_t totaltime = total_timer.getTime();
		general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);
		//print times from the devices...
		logger.info() << "## Device: Solver";
		(gaugefield.get_task_solver())->print_copy_times(totaltime);

		//CP: this is just a fist version and will go into an own file later
		string profiling_out;
		profiling_out = string(argv[0]) + string("_profiling_data");

		fstream prof_file;
		prof_file.open(profiling_out.c_str(), std::ios::out | std::ios::app);
		if(prof_file.is_open()) {
			parameters.print_info_inverter(argv[0], &prof_file);
			prof_file.close();
		} else {
			logger.warn() << "Could not open " << profiling_out;
		}

		//print only dslash-infos
		const char * kernelName;
		kernelName = "dslash_eo";
		gaugefield.get_task_solver()->Opencl_Module::print_profiling(profiling_out, kernelName, (*gaugefield.get_task_solver()->get_timer(kernelName)).getTime(), (*gaugefield.get_task_solver()->get_timer(kernelName)).getNumMeas(), gaugefield.get_task_solver()->get_read_write_size(kernelName), gaugefield.get_task_solver()->get_flop_size(kernelName)) ;

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
}
