#include "inverter.h"

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

	if(parameters.get_num_dev() == 1) {
#ifdef _USEGPU_
		devicetypes[0] = CL_DEVICE_TYPE_GPU;
#else
		devicetypes[0] = CL_DEVICE_TYPE_CPU;
#endif
	} else if(parameters.get_num_dev() == 2) {
		devicetypes[0] = CL_DEVICE_TYPE_GPU;
		devicetypes[1] = CL_DEVICE_TYPE_CPU;
	} else {
		logger.fatal() << "Number of devices too big, aborting..." ;
		return HMC_STDERR;
	}
	logger.trace() << "init gaugefield" ;
	gaugefield.init(parameters.get_num_dev(), devicetypes, &parameters);
	//cerr << "print initial gaugeobservables..." << endl;
	//  gaugefield.print_gaugeobservables(&polytime, &plaqtime);
	logger.trace() << "copy gaugefield" ;
	gaugefield.copy_gaugefield_to_devices(&copy_to_from_dev_timer);
	init_timer.add();
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// inverter
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	perform_timer.reset();
	logger.info() << "perform inversion on device.." ;
	if(parameters.get_num_dev() == 1) {
		gaugefield.perform_inversion_pointsource_ps_corr_devices(&copy_to_from_dev_timer, &copy_on_dev_timer,&solver_timer);
		/** @todo improve ls, gs, here*/
		gaugefield.get_devices_fermions()[0].ps_correlator_device(1,1);
	}
	else{	
		gaugefield.perform_inversion_pointsource_ps_corr_devices(&copy_to_from_dev_timer,&copy_on_dev_timer,&solver_timer);
		gaugefield.get_devices_fermions()[0].get_spinorfield_from_device(host_spinorfield, &copy_to_from_dev_timer);
		gaugefield.get_devices_fermions()[1].copy_spinorfield_to_device(host_spinorfield, &copy_to_from_dev_timer);
		/** @todo improve ls, gs, here*/
		gaugefield.get_devices_fermions()[1].ps_correlator_device(1, 1);
	}
	logger.trace() << "inversion done" ;
	perform_timer.add();
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	total_timer.add();
	inverter_time_output(&total_timer, &init_timer, &perform_timer, &copy_to_from_dev_timer, &copy_on_dev_timer, &solver_timer);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	gaugefield.finalize();

	return HMC_SUCCESS;
}
