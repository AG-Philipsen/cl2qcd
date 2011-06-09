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

	Gaugefield_inversion gaugefield;
	cl_device_type devicetypes[1];

#ifdef _USEGPU_
	devicetypes[0] = CL_DEVICE_TYPE_GPU;
#else
	devicetypes[0] = CL_DEVICE_TYPE_CPU;
#endif
	cerr << "init gaugefield" << endl;
	gaugefield.init(1, devicetypes, &parameters, &inittime);
	cerr << "gaugefield has values: "<<endl;
	gaugefield.print_gaugeobservables (&copytimer, &copytimer);
	gaugefield.copy_gaugefield_to_devices(&copytimer);
	gaugefield.sync_gaugefield(&copytimer);
	cerr << "gaugefield has values: "<<endl;
	gaugefield.print_gaugeobservables (&copytimer, &copytimer);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// inverter
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	gaugefield.perform_inversion_pointsource_ps_corr_host();
	gaugefield.perform_inversion_pointsource_ps_corr_devices(&copytimer,&singletimer,&Mtimer,&scalarprodtimer,&latimer,&dslashtimer,&Mdiagtimer,&solvertimer);


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	gaugefield.finalize();

	return HMC_SUCCESS;
}
