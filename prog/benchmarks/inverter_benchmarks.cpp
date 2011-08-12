#include "../inverter.h"

int main(int argc, char* argv[])
{

	if(argc != 2) {
	  logger.fatal() << "need file name for input parameters";
	  return HMC_FILEERROR;
	}

	char* inputfile = argv[1];
	inputparameters parameters;
	parameters.readfile(inputfile);
	parameters.print_info_inverter(argv[0]);

	ofstream ofile;
	ofile.open("inverter.log");
	if(ofile.is_open()) {
	  parameters.print_info_inverter(argv[0],&ofile);
	  ofile.close();
	} else {
	  logger.warn() << "Could not log file for inverter.";
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	init_timer.reset();
	sourcefileparameters parameters_source;
	//CP: spinorfield on host for storage while copying between devices...
	spinorfield host_spinorfield [SPINORFIELDSIZE];

	Gaugefield_inversion gaugefield;
	cl_device_type devicetypes[parameters.get_num_dev()];
	gaugefield.init_devicetypes_array(devicetypes, &parameters);

	logger.trace() << "init gaugefield" ;
	gaugefield.init(parameters.get_num_dev(), devicetypes, &parameters);
	logger.trace()<< "initial gaugeobservables:";
	gaugefield.print_gaugeobservables(&poly_timer, &plaq_timer);
	logger.trace() << "copy gaugefield" ;
	gaugefield.copy_gaugefield_to_devices(&copy_to_from_dev_timer);
	init_timer.add();
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// inverter-benchmarks
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	perform_timer.reset();

	int hmc_iter = parameters.get_hmcsteps();
	int iter;	
	
	logger.trace() << "Perform " << hmc_iter << "of benchmarking";
	for(iter = 0; iter < hmc_iter; iter ++) {
		/** @todo Insert functions here */
	}
	logger.trace() << "inverter-benchmarking done" ;
	perform_timer.add();
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	total_timer.add();
	general_time_output(&total_timer, &init_timer, &perform_timer, &copy_to_from_dev_timer, &copy_on_dev_timer, &plaq_timer, &poly_timer);
	
	//CP: this is just a fist version and will go into an own file later
	stringstream profiling_out;
	profiling_out << argv[0] << "_profiling_data";

	fstream prof_file;
	prof_file.open(profiling_out.str(), std::ios::out | std::ios::app);
	if(prof_file.is_open()) {
	  parameters.print_info_heatbath(argv[0], &prof_file);
	  prof_file.close();
	} else {
	  logger.warn()<<"Could not open " << profiling_out;
	}
	gaugefield.get_devices()[0].print_profiling(profiling_out.str());
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int err = gaugefield.finalize();
	if (err!= HMC_SUCCESS) 
		logger.fatal() << "error in finalizing " << argv[0];
	return HMC_SUCCESS;
	
}
