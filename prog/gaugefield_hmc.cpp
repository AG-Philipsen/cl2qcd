#include "gaugefield_hmc.h"
#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

hmc_error Gaugefield_hmc::init_devices(cl_device_type* devicetypes, usetimer* timer)
{
// 	if(get_num_ocl_devices() != 1) {
// 		//LZ: so far, we only use !!! 1 !!! device
// 		//this needs generalisation to several devices and subsets!!!!!
// 		cerr << "only 1 device possible..." << endl;
// 	}
	if(get_num_ocl_devices() > 0) {
		Opencl_hmc* dev_tmp = new Opencl_hmc[this->get_num_ocl_devices()];
		alloc_devicetypes();
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

hmc_error Gaugefield_hmc::perform_hmc_step(int dev, inputparameters *parameters, hmc_observables *obs, int iter, hmc_float rnd_number, const string outname, usetimer* copytimer, usetimer* singletimer, usetimer* Mtimer, usetimer* scalarprodtimer, usetimer* latimer, usetimer* dslashtimer, usetimer* Mdiagtimer, usetimer* solvertimer){
	
	//global and local work sizes;
	//LZ: should eventually be moved inside opencl_fermions class
#ifdef _USEGPU_
	const size_t ls = NUMTHREADS; /// @todo have local work size depend on kernel properties (and device? autotune?)
#else
	const size_t ls = 1; // nothing else makes sense on CPU
#endif

#ifdef _USEGPU_
	size_t gs = 4 * NUMTHREADS * get_devices()[dev].max_compute_units; /// @todo autotune
#else
	size_t gs = get_devices()[dev].max_compute_units;
#endif

	const cl_uint num_groups = (gs + ls - 1) / ls;
	gs = ls * num_groups;
	
	//variables for interesting observables
	hmc_float deltah;
	hmc_float plaq;
	hmc_float poly;
	
	/////////////////////////////////////////////////////////////////////
	//HMC-algorithm
	
	//init gauge_momenta, saved in clmem_p
	logger.debug() << "\tinit gauge momentum" ;
	get_devices_hmc()[dev].generate_gaussian_gaugemomenta_device(ls, gs, latimer);
	
	//init/update spinorfield phi
	logger.debug() << "\tinit spinorfield " ;
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
	get_devices_hmc()[dev].generate_gaussian_spinorfield_device(ls, gs, latimer);
	get_devices_hmc()[dev].calc_spinorfield_init_energy_device(ls, gs, scalarprodtimer);
	logger.debug() << "\tperform md update of spinorfield" ;
	get_devices_hmc()[dev].md_update_spinorfield_device(ls, gs, Mtimer);
	
	//update gaugefield and gauge_momenta via leapfrog
	//here, clmem_phi is inverted several times and stored in clmem_phi_inv
	logger.debug() << "\tperform leapfrog to update gaugefield and gaugemomentum" ;
	
	/** @todo these have to be reconsidered! */
	//copy u->u' p->p' for the leapfrog
	get_devices_hmc()[dev].copy_gaugefield_old_new_device(ls, gs, copytimer);
	get_devices_hmc()[dev].copy_gaugemomenta_old_new_device(ls, gs, copytimer);
		
	get_devices_hmc()[dev].leapfrog_device((*parameters).get_tau(), (*parameters).get_integrationsteps1(), (*parameters).get_integrationsteps2(), ls, gs, latimer);
		
// 	logger.debug() << "\tobservables of new config:\t" ;
// 		print_gaugeobservables(new_field, &polytime, &plaqtime);

	//metropolis step: afterwards, the updated config is again in gaugefield and p
	logger.debug() << "\tperform Metropolis step: " ;
	//this call calculates also the HMC-Observables
	*obs = get_devices_hmc()[dev].metropolis(rnd_number, (*parameters).get_beta(), outname, ls, gs, latimer);

	if((*obs).accept == 1){
		// perform the change nonprimed->primed !
		get_devices_hmc()[dev].copy_gaugefield_new_old_device(ls, gs, latimer);
		get_devices_hmc()[dev].copy_gaugemomenta_new_old_device(ls, gs, latimer);
		// SL: this works as long as p and field are pointers to the *original* memory locations!
		logger.debug() << "\t\tnew configuration accepted" ;
	}
	else{
		logger.debug() << "\t\tnew configuration rejected" ;
	}
	logger.trace()<< "\tfinished HMC trajectory " << iter ;
	
	return HMC_SUCCESS;
}

void Gaugefield_hmc::print_hmcobservables(hmc_observables obs, int iter, std::string filename)
{
	hmc_float exp_deltaH = exp(obs.deltaH);
	logger.trace() << "Observables: " << obs.plaq << "\t" << obs.tplaq << "\t" << obs.splaq << "\t" << obs.poly.re << "\t" << obs.poly.im <<  "\t" << obs.deltaH << "\t" << exp_deltaH << "\t" << obs.prob << "\t" << obs.accept ;
// 	printf("Observables:%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n",iter,obs.plaq,obs.tplaq,obs.splaq,obs.poly.re,obs.poly.im,obs.deltaH, exp_deltaH, obs.prob, obs.accept );
	std::fstream hmcout;
	hmcout.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!hmcout.is_open()) exit(HMC_FILEERROR);
	hmcout.width(8);
	hmcout << iter;
	hmcout << "\t";
	hmcout.precision(15);
	hmcout << obs.plaq << "\t" << obs.tplaq << "\t" << obs.splaq << "\t" << obs.poly.re << "\t" << obs.poly.im << "\t" << sqrt(obs.poly.re * obs.poly.re + obs.poly.im * obs.poly.im) <<  "\t" << obs.deltaH << "\t" << exp_deltaH << "\t" << obs.prob << "\t" << obs.accept << std::endl;
	hmcout.close();
	return;
}

