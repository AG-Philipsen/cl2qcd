#include "hmc.h"

//this is the stand-alone HMC
//this most likely has to be changed according to the bigger one

int main(int argc, char* argv[])
{
	char* progname = argv[0];
	print_hello(progname);

	if( argc != 2 ) {
		std::cerr << "Please specify one input file." << std::endl;
		return HMC_FILEERROR;
	}

	char* inputfile = argv[1];
	inputparameters parameters;
	parameters.readfile(inputfile);
	print_info(&parameters);

	stringstream gaugeout_name;
	gaugeout_name << "gaugeobservables_beta" << parameters.get_beta();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	sourcefileparameters parameters_source;
	hmc_gaugefield * gaugefield;
	gaugefield = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
	hmc_rndarray rndarray;

	//init gaugefield according to the inputfile settings
	init_gaugefield(gaugefield, &parameters, &inittime);
	int err = init_random_seeds(rndarray, "rand_seeds", &inittime);
	if( err )
		return err;

	//TODO add a security function that ends the OpenCL-init if it takes too long (this is apparently necessary on the loewe)
#ifdef _USEGPU_
	opencl device(CL_DEVICE_TYPE_GPU, local_work_size, global_work_size, &inittime, &parameters);
#else /* _USEGPU_ */
	opencl device(CL_DEVICE_TYPE_CPU, local_work_size, global_work_size, &inittime, &parameters);
#endif /* _USEGPU_ */
	cout << endl << "OpenCL initialisaton time:\t" << inittime.getTime() << " [mus]" << endl;

	cout << "initial values of observables:\n\t" ;
	print_gaugeobservables(gaugefield, &polytime, &plaqtime);

	device.copy_gaugefield_to_device(gaugefield, &copytime);
	device.copy_rndarray_to_device(rndarray, &copytime);

	//TODO CP: port to OpenCL *g*
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Hybrid Monte Carlo
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//init section
	//CP: here, one has to get all the parameters necessary from the inputparameters
	// suppose that in parameters the number of hmc_iterations is given, store them in a variable here...
	int hmc_iter; //= ...
	int iter;
	int err = 0;
	//beta has to be saved here to give it to the metropolis step, all other parameters can be given via parameters
	hmc_float beta;
	// TODO give seed meaningful value, perhaps read it in from parameters
	int seed = 89343894543895;
	Random hmc_rnd_gen (seed);
	hmc_float rnd_number;

	hmc_spinor_field* phi = new hmc_spinor_field[SPINORFIELDSIZE];
	hmc_spinor_field* phi_inv = new hmc_spinor_field[SPINORFIELDSIZE];
	hmc_spinor_field* chi = new hmc_spinor_field[SPINORFIELDSIZE];
	hmc_gauge_momentum* p = new hmc_gauge_momentum[GAUGEMOMENTASIZE];
	hmc_gauge_momentum* new_p = new hmc_gauge_momentum[GAUGEMOMENTASIZE];
	//the "old" field is the "gaugefield introduced above
	hmc_gaugefield new_field;

	//TODO add implicit measurements(plaquette, deltaH, ...)
	//main hmc-loop
	for(iter = 0; iter < hmc_iter; iter ++) {
		//init gauge_momenta
		//TODO perhaps write a wrapper that automatically evaluates the err's
		err = generate_gaussian_gauge_momenta(p);

		//init/update spinorfield phi
		err = generate_gaussian_spinorfield(chi);
		err = md_update_spinorfield(chi, phi, &gaugefield, parameters);

		//update gaugefield and gauge_momenta via leapfrog
		//here, phi is inverted several times and stored in phi_inv each time
// 		err = leapfrog(parameters, &gaugefield, p, phi, &new_field, new_p, phi_inv);

		//metropolis step: afterwards, the updated config is again in gaugefield and p
		//generate new random-number
		rnd_number = hmc_rnd_gen.doub();
		err = metropolis(rnd_number, beta, phi, phi_inv, &gaugefield, p, &new_field, new_p);

		//TODO if wanted, measurements can be added here...
	}

	delete [] phi;
	delete [] phi_inv;
	delete [] chi;
	delete [] p;
	delete [] new_p;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	totaltime.add();
	save_gaugefield(gaugefield, &parameters, nsteps);
	time_output(
	  &totaltime, &inittime, &polytime, &plaqtime, &updatetime, &overrelaxtime, &copytime
#ifdef _FERMIONS_
	  , &inittimer, &singletimer, &Mtimer, &copytimer, &scalarprodtimer, &latimer, &solvertimer, &dslashtimer, &Mdiagtimer
#endif /* _FERMIONS_ */
	);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	free(gaugefield);
	device.finalize();

	return HMC_SUCCESS;
}
