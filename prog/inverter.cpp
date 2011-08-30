#include "inverter.h"

int main(int argc, char* argv[])
{
  try{

    if(argc != 2) throw Print_Error_Message("Need file name for input parameters",__FILE__,__LINE__);

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
	cl_device_type* devicetypes = new cl_device_type[parameters.get_num_dev()];

	logger.trace() << "init gaugefield" ;
	gaugefield.init(parameters.get_num_dev(), devicetypes, &parameters);

	delete [] devicetypes;

	logger.trace()<< "initial gaugeobservables:";
	gaugefield.print_gaugeobservables(&poly_timer, &plaq_timer);
	logger.trace() << "copy gaugefield" ;
	gaugefield.copy_gaugefield_to_devices();
	init_timer.add();
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// inverter
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	perform_timer.reset();
	/** @todo usage of solver_timer has to be checked. No output yet */
	usetimer solver_timer;
	logger.info() << "perform inversion on device.." ;
	if(parameters.get_num_dev() == 1) {
		gaugefield.perform_inversion_pointsource_ps_corr_devices(&solver_timer);
		/** @todo improve ls, gs, here*/
		gaugefield.get_devices_fermions()[0].ps_correlator_device(gaugefield.get_devices_fermions()[0].get_clmem_corr());
	}
	else{	
		gaugefield.perform_inversion_pointsource_ps_corr_devices(&solver_timer);
		gaugefield.get_devices_fermions()[0].get_buffer_from_device(gaugefield.get_devices_fermions()[0].get_clmem_corr(), host_spinorfield, sizeof(spinor)*SPINORFIELDSIZE);		
		gaugefield.get_devices_fermions()[1].copy_buffer_to_device(host_spinorfield, gaugefield.get_devices_fermions()[0].get_clmem_corr(), sizeof(spinor)*SPINORFIELDSIZE);
		/** @todo improve ls, gs, here*/
		gaugefield.get_devices_fermions()[1].ps_correlator_device(gaugefield.get_devices_fermions()[0].get_clmem_corr());
	}
	logger.trace() << "inversion done" ;
	perform_timer.add();
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	total_timer.add();
	general_time_output(&total_timer, &init_timer, &perform_timer, gaugefield.get_devices_fermions()[0].get_copy_to(), gaugefield.get_devices_fermions()[0].get_copy_on(), &plaq_timer, &poly_timer);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	gaugefield.finalize();

  } //try
  //exceptions from Opencl classes
  catch (Opencl_Error& e) {
    logger.fatal()<<e.what();
    exit(1);
  }
  catch (File_Exception& fe) {
    logger.fatal()<<"Could not open file: "<<fe.get_filename();
    logger.fatal()<<"Aborting.";
    exit(1);
  }
  catch (Print_Error_Message& em) {
    logger.fatal()<<em.what();
    exit(1);
  }

	return 0;
	
}
