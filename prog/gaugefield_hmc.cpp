#include "gaugefield_hmc.h"
#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

void Gaugefield_hmc::init_devices(cl_device_type* devicetypes)
{
//  if(get_num_ocl_devices() != 1) {
//    //LZ: so far, we only use !!! 1 !!! device
//    //this needs generalisation to several devices and subsets!!!!!
//    cerr << "only 1 device possible..." << endl;
//  }
	if(get_num_ocl_devices() > 0) {
		Opencl_hmc* dev_tmp = new Opencl_hmc[this->get_num_ocl_devices()];
		alloc_devicetypes();
		set_devices(dev_tmp);
	}

	for(int n = 0; n < get_num_ocl_devices(); n++) {
		logger.debug() << "init device #" << n;
		get_devices_hmc()[n].init(devicetypes[n], get_parameters(), this->get_numrndstates());
	}
	return;
}

void Gaugefield_hmc::finalize()
{
	Gaugefield_inversion::finalize();
	for(int n = 0; n < get_num_ocl_devices(); n++)
		get_devices_hmc()[n].finalize_hmc();
	return;
}

void Gaugefield_hmc::free_devices()
{
	if(get_num_ocl_devices() > 0)
		delete [] get_devices_hmc();
	return;
}

Opencl_hmc * Gaugefield_hmc::get_devices_hmc ()
{
	return  (Opencl_hmc*)get_devices();
}

void Gaugefield_hmc::perform_hmc_step(int dev, inputparameters *parameters, hmc_observables *obs, int iter, hmc_float rnd_number)
{
	/////////////////////////////////////////////////////////////////////
	//HMC-algorithm

	//init gauge_momenta, saved in clmem_p
	logger.debug() << "\tinit gauge momentum" ;
	get_devices_hmc()[dev].generate_gaussian_gaugemomenta_device();

	//init/update spinorfield phi
	logger.debug() << "\tinit spinorfield " ;
	//NOTE: one does not have to use phi as initial spinorfield in order to save one variable!!!
	//  original alg:
	//    generate_gaussian_spinorfield(chi)
	//    energy_init = |chi|^2
	//    md_update_spinorfield_device(chi, phi): phi = Qminus chi
	//  this can be changed to:
	//    generate_gaussian_spinorfield(phi_inv)
	//    energy_init = |phi_inv|^2
	//    md_update_spinorfield_device(phi_inv, phi): phi = Qminus phi_inv
	//  saving one variable in global mem!!
	get_devices_hmc()[dev].generate_gaussian_spinorfield_device();
	get_devices_hmc()[dev].calc_spinorfield_init_energy();
	logger.debug() << "\tperform md update of spinorfield" ;
	get_devices_hmc()[dev].md_update_spinorfield();

	//update gaugefield and gauge_momenta via leapfrog
	//here, clmem_phi is inverted several times and stored in clmem_phi_inv
	logger.debug() << "\tperform leapfrog to update gaugefield and gaugemomentum" ;

	//copy u->u' p->p' for the leapfrog
	get_devices_hmc()[dev].copy_buffer_on_device(get_devices_hmc()[dev].get_clmem_gaugefield(), get_devices_hmc()[dev].get_clmem_new_u(), NDIM * parameters->get_volspace() * NTIME * sizeof(Matrixsu3));
	get_devices_hmc()[dev].copy_buffer_on_device(get_devices_hmc()[dev].get_clmem_p(), get_devices_hmc()[dev].get_clmem_new_p(), sizeof(ae) * get_parameters()->get_gaugemomentasize());
	///@todo this timer is not used at the moment, compare to inverter.cpp
	usetimer solvertimer;
	get_devices_hmc()[dev].leapfrog((*parameters).get_tau(), (*parameters).get_integrationsteps1(), (*parameters).get_integrationsteps2(), &solvertimer);

	//metropolis step: afterwards, the updated config is again in gaugefield and p
	logger.debug() << "\tperform Metropolis step: " ;
	//this call calculates also the HMC-Observables
	*obs = get_devices_hmc()[dev].metropolis(rnd_number, (*parameters).get_beta());

	if((*obs).accept == 1) {
		// perform the change nonprimed->primed !
		get_devices_hmc()[dev].copy_buffer_on_device(get_devices_hmc()[dev].get_clmem_new_u(), get_devices_hmc()[dev].get_clmem_gaugefield(), NDIM * parameters->get_volspace() * NTIME * sizeof(Matrixsu3));
		get_devices_hmc()[dev].copy_buffer_on_device(get_devices_hmc()[dev].get_clmem_new_p(), get_devices_hmc()[dev].get_clmem_p(), sizeof(ae) * get_parameters()->get_gaugemomentasize());
		logger.debug() << "\t\tnew configuration accepted" ;
	} else {
		logger.debug() << "\t\tnew configuration rejected" ;
	}
	logger.trace() << "\tfinished HMC trajectory " << iter ;

	return;
}

void Gaugefield_hmc::print_hmcobservables(hmc_observables obs, int iter, std::string filename)
{
	hmc_float exp_deltaH = exp(obs.deltaH);
	logger.trace() << "Observables: " << obs.plaq << "\t" << obs.tplaq << "\t" << obs.splaq << "\t" << obs.poly.re << "\t" << obs.poly.im <<  "\t" << obs.deltaH << "\t" << exp_deltaH << "\t" << obs.prob << "\t" << obs.accept ;
//  printf("Observables:%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n",iter,obs.plaq,obs.tplaq,obs.splaq,obs.poly.re,obs.poly.im,obs.deltaH, exp_deltaH, obs.prob, obs.accept );
	std::fstream hmcout;
	hmcout.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!hmcout.is_open()) throw File_Exception(filename);
	hmcout.width(8);
	hmcout << iter;
	hmcout << "\t";
	hmcout.precision(15);
	hmcout << obs.plaq << "\t" << obs.tplaq << "\t" << obs.splaq << "\t" << obs.poly.re << "\t" << obs.poly.im << "\t" << sqrt(obs.poly.re * obs.poly.re + obs.poly.im * obs.poly.im) <<  "\t" << obs.deltaH << "\t" << exp_deltaH << "\t" << obs.prob << "\t" << obs.accept << std::endl;
	hmcout.close();
	return;
}

