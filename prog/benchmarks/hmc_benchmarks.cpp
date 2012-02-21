#include "../hmc.h"

int main(int argc, char* argv[])
{

	if(argc != 2) {
		logger.fatal() << "need file name for input parameters";
		throw File_Exception("No file given");
	}

	char* inputfile = argv[1];
	inputparameters parameters;
	parameters.readfile(inputfile);
	parameters.print_info_hmc(argv[0]);

	//name of file to store gauge observables, print initial information
	/** @todo think about what is a senseful filename*/
	stringstream gaugeout_name;
	gaugeout_name << "hmc_output";

	fstream logfile;
	logfile.open("hmc.log", std::ios::out | std::ios::app);
	if(logfile.is_open()) {
		parameters.print_info_hmc(argv[0], &logfile);
		logfile.close();
	} else {
		logger.warn() << "Could not open hmc.log";
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//cl_int err;

	init_timer.reset();
	sourcefileparameters parameters_source;
	//hmc_observables obs;

	Gaugefield_hmc gaugefield;

	int numtasks = 1;
	cl_device_type primary_device;
	switch ( parameters.get_use_gpu() ) {
		case true :
			primary_device = CL_DEVICE_TYPE_GPU;
			break;
		case false :
			primary_device = CL_DEVICE_TYPE_CPU;
			break;
	}

	logger.trace() << "init gaugefield" ;
	gaugefield.init(numtasks, primary_device, &parameters);
	logger.trace() << "initial gaugeobservables:";
	gaugefield.print_gaugeobservables(0);
	init_timer.add();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// HMC-Benchmarks
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	perform_timer.reset();
	//CP: this is taken from gaugefield_hmc. However, I took out all parts out of the loop that have to do with changing the fields in order to always perform the same HMC step. 
	int hmc_iter = parameters.get_hmcsteps();
	int iter;
	//This is the random-number generator for the metropolis-step
  Random hmc_rnd_gen (parameters.get_host_seed());

	logger.debug() << "\tinit spinorfield and gaugemomentum" ;
	gaugefield.init_gaugemomentum_spinorfield();
	
	logger.debug() << "\tupdate gaugefield and gaugemomentum" ;
	size_t gfsize = get_parameters()->get_gf_buf_size();
	size_t gmsize = get_parameters()->get_gm_buf_size();	
	//copy u->u' p->p' for the integrator
	gaugefield.get_task_hmc(0)->copy_buffer_on_device(*(gaugefield.get_task_hmc(0)->get_gaugefield()), gaugefield.get_task_hmc(0)->get_clmem_new_u(), gfsize);
	gaugefield.get_task_hmc(0)->copy_buffer_on_device(gaugefield.get_task_hmc(0)->get_clmem_p(), gaugefield.get_task_hmc(0)->get_clmem_new_p(), gmsize);	
	
	logger.trace() << "Perform " << hmc_iter << "of benchmarking";
	for(iter = 0; iter < hmc_iter; iter ++) {
		//here, clmem_phi is inverted several times and stored in clmem_phi_inv
		gaugefield.integrator(solver_timer);
		//metropolis step: afterwards, the updated config is again in gaugefield and p
		logger.debug() << "\tperform Metropolis step: " ;
		//this call calculates also the HMC-Observables
		*obs = gaugefield.get_task_hmc(0)->metropolis(rnd_number, gaugefield.get_parameters()->get_beta());
		//CP: just reject the outcome of the metropolis step
		logger.trace() << "\tfinished HMC trajectory " << iter ;
	}
	logger.trace() << "HMC-benchmarking done";
	perform_timer.add();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	total_timer.add();
	general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);

	//CP: this is just a fist version and will go into an own file later
	stringstream profiling_out;
	profiling_out << argv[0] << "_profiling_data";

	fstream prof_file;
	prof_file.open(profiling_out.str(), std::ios::out | std::ios::app);
	if(prof_file.is_open()) {
		parameters.print_info_heatbath(argv[0], &prof_file);
		prof_file.close();
	} else {
		logger.warn() << "Could not open " << profiling_out;
	}
	gaugefield.print_profiling(profiling_out.str());

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	gaugefield.finalize();

	return 0;
}
