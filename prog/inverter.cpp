#include "inverter.h"

int main(int argc, char* argv[])
{

	if(argc != 2) {
		cerr << "need file name for input parameters" << endl;
		return HMC_FILEERROR;
	}

	char* progname = argv[0];
	print_hello(progname);

	char* inputfile = argv[1];
	inputparameters parameters;
	parameters.readfile(inputfile);
	print_info(&parameters,&cout);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	sourcefileparameters parameters_source;
	//CP: spinorfield on host for storage while copying between devices...
	spinorfield host_spinorfield [SPINORFIELDSIZE];

	Gaugefield_inversion gaugefield;
	cl_device_type devicetypes[parameters.get_num_dev()];

	if(parameters.get_num_dev() == 1){
	#ifdef _USEGPU_
		devicetypes[0] = CL_DEVICE_TYPE_GPU;
	#else
		devicetypes[0] = CL_DEVICE_TYPE_CPU;
	#endif
	}
	else if(parameters.get_num_dev() == 2){
		devicetypes[0] = CL_DEVICE_TYPE_GPU;
		devicetypes[1] = CL_DEVICE_TYPE_CPU;
	}
	else{
		logger.fatal() << "Number of devices too big, aborting..." ;
		return HMC_STDERR;
	}
	cerr << "init gaugefield" << endl;
	gaugefield.init(parameters.get_num_dev(), devicetypes, &parameters, &inittime);
	//cerr << "print initial gaugeobservables..." << endl;
	//	gaugefield.print_gaugeobservables(&polytime, &plaqtime);
	cerr << "copy gaugefield" << endl;
	gaugefield.copy_gaugefield_to_devices(&copytimer);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// inverter
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "perform inversion on device.." << endl;

	if(parameters.get_num_dev() == 1){
		gaugefield.perform_inversion_pointsource_ps_corr_devices(&copytimer,&singletimer,&Mtimer,&scalarprodtimer,&latimer,&dslashtimer,&Mdiagtimer,&solvertimer);
		/** @todo improve ls, gs, here*/
		gaugefield.get_devices_fermions()[0].ps_correlator_device(1,1, &latimer);
	}
	else{
		gaugefield.perform_inversion_pointsource_ps_corr_devices(&copytimer,&singletimer,&Mtimer,&scalarprodtimer,&latimer,&dslashtimer,&Mdiagtimer,&solvertimer);
		gaugefield.get_devices_fermions()[0].get_spinorfield_from_device(host_spinorfield, &copytimer);
		gaugefield.get_devices_fermions()[1].copy_spinorfield_to_device(host_spinorfield, &copytimer);
		/** @todo improve ls, gs, here*/
		gaugefield.get_devices_fermions()[1].ps_correlator_device(1,1, &latimer);
	}
	cout << "inversion done" << endl;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	gaugefield.finalize();

	return HMC_SUCCESS;
}
