#include "../inverter.h"

int main(int argc, char* argv[])
{
	try {

		if(argc != 2) {
			logger.fatal() << "need file name for input parameters";
			throw File_Exception("No file given");
		}

		char* inputfile = argv[1];
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

		//one needs 1 task here
		int numtasks = 1;
		if(parameters.get_num_dev() != 2 )
			logger.warn() << "Only 1 device demanded by benchmark executable. All calculations performed on primary device.";

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
		size_t eoprec_spinorfield_size = sizeof(spinor) * gaugefield.get_task_solver()->get_parameters()->get_eoprec_spinorfieldsize();
		size_t gf_size = gaugefield.get_task_solver()->get_parameters()->get_gf_buf_size();
		cl_mem sf1, sf2, gf;

		gf = gaugefield.get_task_solver()->create_rw_buffer(gf_size);
		sf1 = gaugefield.get_task_solver()->create_rw_buffer(eoprec_spinorfield_size);
		sf2 = gaugefield.get_task_solver()->create_rw_buffer(eoprec_spinorfield_size);

		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// dslash-benchmark
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();

		int hmc_iter = parameters.get_hmcsteps();
		int iter;

		logger.trace() << "Perform " << hmc_iter << "of dslash benchmarking (EVEN + ODD) for each step";
		for(iter = 0; iter < hmc_iter; iter ++) {
			gaugefield.get_task_solver()->dslash_eoprec_device(sf1, sf2, gf, EVEN);
			gaugefield.get_task_solver()->dslash_eoprec_device(sf1, sf2, gf, ODD);
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

		//print only dslash-infos
		const char * kernelName;
		kernelName = "dslash_eoprec";
		gaugefield.get_task_solver()->Opencl_Module::print_profiling(profiling_out.str(), kernelName, (*gaugefield.get_task_solver()->get_timer(kernelName)).getTime(), (*gaugefield.get_task_solver()->get_timer(kernelName)).getNumMeas(), gaugefield.get_task_solver()->get_read_write_size(kernelName), gaugefield.get_task_solver()->get_flop_size(kernelName)) ;

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
