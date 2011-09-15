#include "gaugefield_hmc_tmp.h"


Opencl_Module_Hmc* Gaugefield_hmc::get_task_hmc(int dev)
{
	//if more than one device is used, the element dev from the array must be called here!!
	return (Opencl_Module_Hmc*)opencl_modules[task_hmc];
}

void Gaugefield_hmc::init_tasks()
{
	task_hmc = 0;

	opencl_modules = new Opencl_Module* [get_num_tasks()];

	//LZ: right now, each task carries exactly one opencl device -> thus the below allocation with [1]. Could be generalized in future
	opencl_modules[task_hmc] = new Opencl_Module_Hmc[1];
	get_task_hmc(0)->init(queue[task_hmc], get_clmem_gaugefield(), get_parameters(), get_max_compute_units(task_hmc), get_double_ext(task_hmc));

	return;
}

void Gaugefield_hmc::delete_variables()
{
	Gaugefield_hybrid::delete_variables();
	return;
}

void Gaugefield_hmc::finalize_opencl()
{

	Gaugefield_hybrid::finalize_opencl();
	return;
}

void Gaugefield_hmc::perform_hmc_step(hmc_observables *obs, int iter, hmc_float rnd_number, usetimer* solver_timer)
{
	size_t gfsize = get_parameters()->get_gf_buf_size();
	size_t gmsize = get_parameters()->get_gm_buf_size();

	//init gauge_momenta, saved in clmem_p
	logger.debug() << "\tinit gauge momentum" ;
	get_task_hmc(0)->generate_gaussian_gaugemomenta_device();

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
	get_task_hmc(0)->generate_gaussian_spinorfield_device();
	get_task_hmc(0)->calc_spinorfield_init_energy();
	logger.debug() << "\tperform md update of spinorfield" ;
	get_task_hmc(0)->md_update_spinorfield();

	//update gaugefield and gauge_momenta via leapfrog
	//here, clmem_phi is inverted several times and stored in clmem_phi_inv
	logger.debug() << "\tperform leapfrog to update gaugefield and gaugemomentum" ;

	//copy u->u' p->p' for the leapfrog
	get_task_hmc(0)->copy_buffer_on_device(*(get_task_hmc(0)->get_gaugefield()), get_task_hmc(0)->get_clmem_new_u(), gfsize);
	get_task_hmc(0)->copy_buffer_on_device(get_task_hmc(0)->get_clmem_p(), get_task_hmc(0)->get_clmem_new_p(), gmsize);

	get_task_hmc(0)->leapfrog(get_parameters()->get_tau(), get_parameters()->get_integrationsteps1(), get_parameters()->get_integrationsteps2(), solver_timer);

	//metropolis step: afterwards, the updated config is again in gaugefield and p
	logger.debug() << "\tperform Metropolis step: " ;
	//this call calculates also the HMC-Observables
	*obs = get_task_hmc(0)->metropolis(rnd_number, get_parameters()->get_beta());

	if((*obs).accept == 1) {
		// perform the change nonprimed->primed !
		get_task_hmc(0)->copy_buffer_on_device(get_task_hmc(0)->get_clmem_new_u(), *(get_task_hmc(0)->get_gaugefield()), gfsize);
		get_task_hmc(0)->copy_buffer_on_device(get_task_hmc(0)->get_clmem_new_p(), get_task_hmc(0)->get_clmem_p(), gmsize);
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