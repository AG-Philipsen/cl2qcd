#include "../hmc.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char* argv[])
{
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
	prng_init(parameters.get_host_seed());
	hmc_float rnd_number = prng_double();
	usetimer solver_timer;
	hmc_observables obs;

	logger.debug() << "\tinit spinorfield and gaugemomentum" ;
	gaugefield.init_gaugemomentum_spinorfield();

	logger.debug() << "\tupdate gaugefield and gaugemomentum" ;
	//copy u->u' p->p' for the integrator
	gaugefield.get_task_hmc(0)->copy_buffer_on_device(gaugefield.get_task_hmc(0)->get_gaugefield(), gaugefield.get_task_hmc(0)->get_clmem_new_u(), gaugefield.get_task_hmc(0)->getGaugefieldBufferSize());
	gaugefield.get_task_hmc(0)->copy_buffer_on_device(gaugefield.get_task_hmc(0)->get_clmem_p(), gaugefield.get_task_hmc(0)->get_clmem_new_p(), gaugefield.get_task_hmc(0)->get_gaugemomentum_buffer_size());
	logger.trace() << "Perform " << hmc_iter << "of benchmarking";
	for(iter = 0; iter < hmc_iter; iter ++) {
		//here, clmem_phi is inverted several times and stored in clmem_phi_inv
		gaugefield.integrator(&solver_timer);
		//metropolis step: afterwards, the updated config is again in gaugefield and p
		logger.debug() << "\tperform Metropolis step: " ;
		//this call calculates also the HMC-Observables
		obs = gaugefield.get_task_hmc(0)->metropolis(rnd_number, gaugefield.get_parameters()->get_beta());
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
	string profiling_out;
	profiling_out = string(argv[0]) + string("_profiling_data");

	fstream prof_file;
	prof_file.open(profiling_out.c_str(), std::ios::out | std::ios::app);
	if(prof_file.is_open()) {
		parameters.print_info_heatbath(argv[0], &prof_file);
		prof_file.close();
	} else {
		logger.warn() << "Could not open " << profiling_out;
	}
	gaugefield.print_profiling(profiling_out);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	gaugefield.finalize();

	return 0;
}
