#include "inverter.h"

void print_hello_bench(char* name)
{
	std::cout << "This is the Inverter-Benchmarking Programme " << name << endl;
	return;
}

int main(int argc, char* argv[])
{

	if(argc != 2) {
		cerr << "need file name for input parameters" << endl;
		return HMC_FILEERROR;
	}

	char* progname = argv[0];
	print_hello_bench(progname);

	char* inputfile = argv[1];
	inputparameters parameters;
	parameters.readfile(inputfile);
	print_info(&parameters,&cout);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	sourcefileparameters parameters_source;

	Gaugefield_inversion gaugefield;
	cl_device_type devicetypes[1];

#ifdef _USEGPU_
	devicetypes[0] = CL_DEVICE_TYPE_GPU;
#else
	devicetypes[0] = CL_DEVICE_TYPE_CPU;
#endif
	cerr << "init gaugefield" << endl;
	gaugefield.init(1, devicetypes, &parameters, &inittime);
	//cerr << "print initial gaugeobservables..." << endl;
	//	gaugefield.print_gaugeobservables(&polytime, &plaqtime);
	cerr << "copy gaugefield" << endl;
	gaugefield.copy_gaugefield_to_devices(&copytimer);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// inverter
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "perform Inverter-Benchmarking on device.." << endl;
	

	cout << "inverter-Benchmarking done" << endl;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	gaugefield.finalize();

	time_output_inverter(
   &totaltime,  &inittime,  &polytime,  &plaqtime,  &updatetime,  &overrelaxtime,  &copytime
  ,  &ferm_inittime, &singletimer, &Mtimer, &copytimer, &scalarprodtimer, &latimer,  &solvertimer,  &dslashtimer,  &Mdiagtimer);
	
	
	return HMC_SUCCESS;
}
