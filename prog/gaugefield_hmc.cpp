#include "gaugefield_hmc.h"


hmc_error Gaugefield_hmc::init_devices(cl_device_type* devicetypes, usetimer* timer)
{
	if(get_num_ocl_devices() != 1) {
		//LZ: so far, we only use !!! 1 !!! device
		//this needs generalisation to several devices and subsets!!!!!
		cerr << "only 1 device possible..." << endl;
	}

	if(get_num_ocl_devices() > 0) {
		Opencl_hmc* dev_tmp = new Opencl_hmc[get_num_ocl_devices()];
		set_devices(dev_tmp);
	}


	for(int n = 0; n < get_num_ocl_devices(); n++) {
		cout << "init device #" << n << endl;
		get_devices_hmc()[n].init(devicetypes[n], timer, get_parameters());
	}
	return HMC_SUCCESS;
}

hmc_error Gaugefield_hmc::finalize()
{
	hmc_error err = HMC_SUCCESS;
	err |= Gaugefield_inversion::finalize();

	return err;
}

hmc_error Gaugefield_hmc::free_devices()
{
	if(get_num_ocl_devices() > 0)
		delete [] get_devices_hmc();
	return HMC_SUCCESS;
}

Opencl_hmc * Gaugefield_hmc::get_devices_hmc ()
{
	return  (Opencl_hmc*)get_devices();
}

hmc_error Gaugefield_hmc::perform_hmc_step(int iter, usetimer* copytimer, usetimer* singletimer, usetimer* Mtimer, usetimer* scalarprodtimer, usetimer* latimer, usetimer* dslashtimer, usetimer* Mdiagtimer, usetimer* solvertimer){
	
	//global and local work sizes;
	//LZ: should eventually be moved inside opencl_fermions class
#ifdef _USEGPU_
	const size_t ls = NUMTHREADS; /// @todo have local work size depend on kernel properties (and device? autotune?)
#else
	const size_t ls = 1; // nothing else makes sense on CPU
#endif

#ifdef _USEGPU_
	size_t gs = 4 * NUMTHREADS * get_devices()[0].max_compute_units; /// @todo autotune
#else
	size_t gs = get_devices()[0].max_compute_units;
#endif

	const cl_uint num_groups = (gs + ls - 1) / ls;
	gs = ls * num_groups;
	
	cout << "\tinit gauge momentum" << endl;
		//init gauge_momenta
// 		get_devices_hmc()[0].generate_gaussian_gaugemomenta_device();
		#ifdef _FERMIONS_
		//init/update spinorfield phi
		cout << "\tinit spinorfield " << endl;
// 		get_devices_hmc()[0].generate_gaussian_spinorfield_device();
		get_devices_hmc()[0].set_float_to_global_squarenorm_device(clmem_chi, clmem_energy_init, ls, gs, scalarprodtimer)
		cout << "\tperform md update of spinorfield" << endl;
// 		get_devices_hmc()[0].md_update_spinorfield_device();
		if(err!=HMC_SUCCESS) {cout << "\t\t\terror: " << err << endl; return HMC_STDERR; }
		#endif
	
		//update gaugefield and gauge_momenta via leapfrog
		//here, phi is inverted several times and stored in phi_inv
		cout << "\tperform leapfrog to update gaugefield and gaugemomentum" << endl;
	
// 		get_devices_hmc()[0].copy_gaugefield_device(gaugefield, new_field);
// 		get_devices_hmc()[0].copy_gaugemomenta_device(p, new_p);
		
// 		get_devices_hmc()[0].leapfrog_device(&parameters, 
// 									 #ifdef _FERMIONS_
// 									 phi, phi_inv,
// 									 #endif
// 									 new_field, new_p
// 									 );
		
		cout << "\tobservables of new config:\n\t" ;
// 		print_gaugeobservables(new_field, &polytime, &plaqtime);
		//metropolis step: afterwards, the updated config is again in gaugefield and p
		cout << "\tperform Metropolis step: " << endl;
		//generate new random-number
// 		rnd_number = hmc_rnd_gen.doub();
// 		get_devices_hmc()[0].metropolis_device(rnd_number, &parameters, 
// 										 #ifdef _FERMIONS_ 
// 										 phi, phi_inv, energy_init, 
// 										 #endif 
// 										 gaugefield, p, new_field, new_p);

		cout<< "\tfinished HMC trajectory " << iter << endl;
		/** @todo CP: measurements should be added here... */
		
// 		print_gaugeobservables(gaugefield, &plaqtime, &polytime, iter, gaugeout_name.str());
	
}

//needed variables for HMC:
/*
variables:
energy_init
chi (gaussian spinorfield)
p, p_new (ae fields)
u_new (gaugefield)
phi, phi_inv (spinorfields)

functions:
generate_gaussian_gaugemomenta_device
generate_gaussian_spinorfield_device
md_update_spinorfield_device()
leapfrog_device
copy_gaugefield_device
copy_gaugemomenta_device
force_device
hamiltonian_device

kernels:
generate_gaussian_spinorfield
generate_gaussian_gaugemomenta
md_update_gaugefield
md_update_gaugemomenta
gauge_force
fermion_force
s_gauge
s_fermion

*/

