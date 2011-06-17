#include "gaugefield_hmc.h"
#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

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
	for(int n = 0; n < get_num_ocl_devices(); n++)
		err |= get_devices_hmc()[n].finalize_hmc();
	
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

hmc_error Gaugefield_hmc::perform_hmc_step(int iter, hmc_float rnd_number, usetimer* copytimer, usetimer* singletimer, usetimer* Mtimer, usetimer* scalarprodtimer, usetimer* latimer, usetimer* dslashtimer, usetimer* Mdiagtimer, usetimer* solvertimer){
	
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
	
	/////////////////////////////////////////////////////////////////////
	//HMC-algorithm
	
	logger.trace() << "\tinit gauge momentum" ;
	//init gauge_momenta
	get_devices_hmc()[0].generate_gaussian_gaugemomenta_device(ls, gs, latimer);
	//init/update spinorfield phi
	logger.trace() << "\tinit spinorfield " ;
	
	//NOTE: one does not have to use phi as initial spinorfield in order to save one variable!!!
	//	original alg:
	//		generate_gaussian_spinorfield(chi)
	//		energy_init = |chi|^2
	//		md_update_spinorfield_device(chi, phi): phi = Qminus chi
	//	this can be changed to:
	//		generate_gaussian_spinorfield(phi_inv)
	//		energy_init = |phi_inv|^2
	//		md_update_spinorfield_device(phi_inv, phi): phi = Qminus phi_inv
	//	saving one variable in global mem!!
	
	get_devices_hmc()[0].generate_gaussian_spinorfield_device(ls, gs, latimer);
	get_devices_hmc()[0].calc_spinorfield_init_energy_device(ls, gs, scalarprodtimer);
	logger.trace() << "\tperform md update of spinorfield" ;
	get_devices_hmc()[0].md_update_spinorfield_device(ls, gs, Mtimer);
	
	//update gaugefield and gauge_momenta via leapfrog
	//here, clmem_phi is inverted several times and stored in clmem_phi_inv
	logger.trace() << "\tperform leapfrog to update gaugefield and gaugemomentum" ;
	
	/** @todo these have to be reconsidered! */
	get_devices_hmc()[0].copy_gaugefield_old_new_device(ls, gs, copytimer);
	get_devices_hmc()[0].copy_gaugemomenta_old_new_device(ls, gs, copytimer);
		
	get_devices_hmc()[0].leapfrog_device(ls, gs, latimer);
		
	logger.trace() << "\tobservables of new config:\t" ;
// 		print_gaugeobservables(new_field, &polytime, &plaqtime);
	//metropolis step: afterwards, the updated config is again in gaugefield and p
	logger.trace() << "\tperform Metropolis step: " ;
	
	//this call calculates deltaH on the device
	get_devices_hmc()[0].hamiltonian_device(ls, gs, latimer);
	hmc_float deltah;
	get_devices_hmc()[0].get_deltah_from_device(&deltah, ls, gs, copytimer);

	/** @todo CP:  export deltah */
	hmc_float compare_prob;
	if(deltah<0){
		compare_prob = exp(deltah);
	}else{
		compare_prob = 1.0;
	}
	
	// SL: the following can be tuned, whether it is more costly to draw always the rnd number even when compare_prob=1
	//     and whether the "if compare_prob==1" costs more or less than always evaluating the exp ...
	if(rnd_number <= compare_prob){
		// perform the change nonprimed->primed !
		get_devices_hmc()[0].copy_gaugefield_new_old_device(ls, gs, latimer);
		get_devices_hmc()[0].copy_gaugemomenta_new_old_device(ls, gs, latimer);
		// SL: this works as long as p and field are pointers to the *original* memory locations!
		logger.trace() << "\t\tnew configuration accepted" ;
	}
	else{
		logger.trace() << "\t\tnew configuration rejected" ;
	}
		logger.trace()<< "\tfinished HMC trajectory " << iter ;
		/** @todo CP: measurements should be added here... */
		
// 		print_gaugeobservables(gaugefield, &plaqtime, &polytime, iter, gaugeout_name.str());
	return HMC_SUCCESS;
}

//needed variables for HMC:
/*
variables:
energy_init
chi (gaussian spinorfield)
p, p_new (ae fields)
u_new (gaugefield)
phi, phi_inv (spinorfields)s

functions:
generate_gaussian_gaugemomenta_device
generate_gaussian_spinorfield_device
md_update_spinorfield_device()
leapfrog_device
copy_gaugefield_device
copy_gaugemomenta_device
force_device
hamiltonian_device
QplusQminus_device (this should go into opencl_fermions...)
Qplus_device (this should go into opencl_fermions...)
Qminus_device (this should go into opencl_fermions...)

kernels:
generate_gaussian_spinorfield
generate_gaussian_gaugemomenta
md_update_gaugefield
md_update_gaugemomenta
gauge_force
fermion_force
s_gauge
s_fermion
gamma5 (this should go into opencl_fermions...)

*/

