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

	logger.trace() << "init gaugefield" ;
	gaugefield.init(parameters.get_num_dev(), devicetypes, &parameters);
	logger.trace()<< "initial gaugeobservables:";
	gaugefield.print_gaugeobservables(&poly_timer, &plaq_timer);
	logger.trace() << "copy gaugefield" ;
	gaugefield.copy_gaugefield_to_devices(&copy_to_from_dev_timer);
	init_timer.add();
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// inverter
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	perform_timer.reset();
	/** @todo usage of solver_timer has to be checked. No output yet */
	usetimer solver_timer;
	logger.info() << "perform inversion on device.." ;
	if(parameters.get_num_dev() == 1) {
		gaugefield.perform_inversion_pointsource_ps_corr_devices(&copy_to_from_dev_timer, &copy_on_dev_timer,&solver_timer);
		/** @todo improve ls, gs, here*/
		gaugefield.get_devices_fermions()[0].ps_correlator_device(gaugefield.get_devices_fermions()[0].get_clmem_corr());
	}
	else{	
		gaugefield.perform_inversion_pointsource_ps_corr_devices(&copy_to_from_dev_timer,&copy_on_dev_timer,&solver_timer);
		gaugefield.get_devices_fermions()[0].get_buffer_from_device(gaugefield.get_devices_fermions()[0].get_clmem_corr(), host_spinorfield, sizeof(spinor)*SPINORFIELDSIZE,  &copy_to_from_dev_timer);		
		gaugefield.get_devices_fermions()[1].copy_buffer_to_device(host_spinorfield, gaugefield.get_devices_fermions()[0].get_clmem_corr(), sizeof(spinor)*SPINORFIELDSIZE, &copy_to_from_dev_timer);
		/** @todo improve ls, gs, here*/
		gaugefield.get_devices_fermions()[1].ps_correlator_device(gaugefield.get_devices_fermions()[0].get_clmem_corr());
	}
	logger.trace() << "inversion done" ;
	perform_timer.add();
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	total_timer.add();
	general_time_output(&total_timer, &init_timer, &perform_timer, &copy_to_from_dev_timer, &copy_on_dev_timer, &plaq_timer, &poly_timer);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int err = gaugefield.finalize();
	if (err!= HMC_SUCCESS) 
		logger.fatal() << "error in finalizing " << argv[0];
	return HMC_SUCCESS;
	
}
