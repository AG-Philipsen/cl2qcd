#include "gaugefield_hmc.h"


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
   /// @todo this must be generalized if more than one device is used for one task
  for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
    opencl_modules[ntask]->finalize();
  }

	Gaugefield_hybrid::finalize_opencl();
	return;
}

//////////////////////////////////////////////////////////
// Methods for HMC-Algorithm

void Gaugefield_hmc::perform_hmc_step(hmc_observables *obs, int iter, hmc_float rnd_number, usetimer* solver_timer)
{
	size_t gfsize = get_parameters()->get_gf_buf_size();
	size_t gmsize = get_parameters()->get_gm_buf_size();
	
	logger.debug() << "\tinit spinorfield and gaugemomentum" ;
	this->init_gaugemomentum_spinorfield();
	
	logger.debug() << "\tupdate gaugefield and gaugemomentum" ;
	//copy u->u' p->p' for the integrator
	get_task_hmc(0)->copy_buffer_on_device(*(get_task_hmc(0)->get_gaugefield()), get_task_hmc(0)->get_clmem_new_u(), gfsize);
	get_task_hmc(0)->copy_buffer_on_device(get_task_hmc(0)->get_clmem_p(), get_task_hmc(0)->get_clmem_new_p(), gmsize);

	//here, clmem_phi is inverted several times and stored in clmem_phi_inv
	this->integrator(solver_timer);

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

void Gaugefield_hmc::print_hmcobservables(hmc_observables obs, int iter)
{
	hmc_float exp_deltaH = exp(obs.deltaH);
	//	logger.info() << setw(8) << setfill(' ') << iter << "\t" << setprecision(15) << obs.plaq << "\t" << obs.tplaq << "\t" << obs.splaq << "\t" << obs.poly.re << "\t" << obs.poly.im << "\t" << sqrt(obs.poly.re * obs.poly.re + obs.poly.im * obs.poly.im) <<  "\t" << obs.deltaH << "\t" << exp_deltaH << "\t" << obs.prob << "\t" << obs.accept;

	//short version of output, all obs are collected in the output file anyways...
	logger.info() << setw(8) << setfill(' ') << iter << "\t" << setprecision(15) << obs.plaq << "\t" << obs.poly.re << "\t" << obs.poly.im << "\t" <<  exp_deltaH;

	return;
}

void Gaugefield_hmc::calc_total_force(usetimer * solvertimer)
{
	//CP: make sure that the output field is set to zero
	get_task_hmc(0)->set_zero_clmem_force_device();
	if(get_parameters()->get_use_smearing() == true) {
		get_task_hmc(0)->smear_gaugefield(*(get_task_hmc(0)->get_gaugefield()));
	}
	get_task_hmc(0)->calc_fermion_force(solvertimer);
	if(get_parameters()->get_use_smearing() == true) {
		get_task_hmc(0)->stout_smeared_fermion_force_device();
		get_task_hmc(0)->unsmear_gaugefield(*(get_task_hmc(0)->get_gaugefield()));
	}
	get_task_hmc(0)->calc_gauge_force();
	/// @todo check again that the force-vector is updated the right way!!
}

void Gaugefield_hmc::md_update_gaugemomentum(hmc_float eps, usetimer * solvertimer){
	calc_total_force(solvertimer);
	get_task_hmc(0)->md_update_gaugemomentum_device(-1.*eps);
	return;
}

void Gaugefield_hmc::md_update_gaugemomentum_gauge(hmc_float eps){
	get_task_hmc(0)->set_zero_clmem_force_device();
	get_task_hmc(0)->calc_gauge_force();
	get_task_hmc(0)->md_update_gaugemomentum_device(-1.*eps);
	return;
}

void Gaugefield_hmc::md_update_gaugemomentum_fermion(hmc_float eps, usetimer * solvertimer){
	get_task_hmc(0)->set_zero_clmem_force_device();
	if(get_parameters()->get_use_smearing() == true) {
		get_task_hmc(0)->smear_gaugefield(*(get_task_hmc(0)->get_gaugefield()));
	}
	get_task_hmc(0)->calc_fermion_force(solvertimer);
	if(get_parameters()->get_use_smearing() == true) {
		get_task_hmc(0)->stout_smeared_fermion_force_device();
		get_task_hmc(0)->unsmear_gaugefield(*(get_task_hmc(0)->get_gaugefield()));
	}
	get_task_hmc(0)->md_update_gaugemomentum_device(-1.*eps);
	return;
}

void Gaugefield_hmc::md_update_gaugefield(hmc_float eps){
	get_task_hmc(0)->md_update_gaugefield_device(eps);
	return;
}

void Gaugefield_hmc::integrator(usetimer * solvertimer){
	if(get_parameters()->get_integrator() == LEAPFROG){
		this->leapfrog(solvertimer);
	}
	else if(get_parameters()->get_integrator() == TWOMN){
		this->twomn(solvertimer);
	}
	return;
}

void Gaugefield_hmc::leapfrog(usetimer * solvertimer)
{
	//it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already when this function is called the first time

	if(get_parameters()->get_num_timescales() == 1) {
		//this is the simplest case, just using 1 timescale
		int steps = get_parameters()->get_integrationsteps1();
		hmc_float stepsize = get_parameters()->get_tau() / ((hmc_float) steps);
		hmc_float stepsize_half = 0.5 * stepsize;

		logger.debug() << "\t\tinitial step:";
		md_update_gaugemomentum(stepsize_half, solvertimer);
		if(steps > 1) logger.debug() << "\t\tperform " << steps - 1 << " intermediate steps " ;
		for(int k = 1; k < steps; k++) {
			md_update_gaugefield(stepsize);
			md_update_gaugemomentum(stepsize, solvertimer);
		}
		logger.debug() << "\t\tfinal step" ;
		md_update_gaugefield(stepsize);
		md_update_gaugemomentum(stepsize_half, solvertimer);
		logger.debug() << "\t\tfinished leapfrog";
	} 
	else if (get_parameters()->get_num_timescales() == 2) {
		int steps1 = get_parameters()->get_integrationsteps1();
		int steps2 = get_parameters()->get_integrationsteps2();
		
		int mult = steps1 % steps2;
		if( mult != 0) Print_Error_Message("integrationsteps1 must be a multiple of integrationssteps2, nothing else is implemented yet. Aborting...");
		
		//this is the number of int.steps for the fermion during one gauge-int.step
		//this is done after hep-lat/0209037. See also hep-lat/050611v2 for a more advanced versions
		int m = steps1/steps2;
		
		//this uses 2 timescales (more is not implemented yet): timescale1 for the gauge-part, timescale2 for the fermion part
		hmc_float stepsize = get_parameters()->get_tau() / ((hmc_float) steps1);
		hmc_float stepsize2 = get_parameters()->get_tau() / ((hmc_float) steps2);
		hmc_float stepsize_half = 0.5 * stepsize;
		hmc_float stepsize2_half = 0.5 * stepsize2;

		logger.debug() << "\t\tinitial step:";
		//this corresponds to V_s2(deltaTau/2)
		md_update_gaugemomentum_fermion(stepsize2_half, solvertimer);
		//now, m steps "more" are performed for the gauge-part
		//this corresponds to [V_s1(deltaTau/2/m) V_t(deltaTau/m) V_s1(deltaTau/2/m) ]^m
		for(int l = 0; l < m; l++) {
			md_update_gaugemomentum_gauge(stepsize_half);
			md_update_gaugefield(stepsize2);
			md_update_gaugemomentum_gauge(stepsize_half);
		}
		if(steps2 > 1) logger.debug() << "\t\tperform " << steps2 - 1 << " intermediate steps " ;
		for(int k = 1; k < steps2; k++) {
			//this corresponds to V_s2(deltaTau)
			md_update_gaugemomentum_fermion(stepsize2, solvertimer);
			for(int l = 0; l < m; l++) {
				//this corresponds to [V_s1(deltaTau/2/m) V_t(deltaTau/m) V_s1(deltaTau/2/m) ]^m
				md_update_gaugemomentum_gauge(stepsize_half);
				md_update_gaugefield(stepsize2);
				md_update_gaugemomentum_gauge(stepsize_half);
			}
		}
		logger.debug() << "\t\tfinal step" ;
		//this corresponds to the missing V_s2(deltaTau/2)
		md_update_gaugemomentum_fermion(stepsize2_half, solvertimer);
		logger.debug() << "\t\tfinished leapfrog";
	}
	else 
		Print_Error_Message("More than 2 timescales is not implemented yet. Aborting...");
}

void Gaugefield_hmc::twomn(usetimer * solvertimer)
{
	//it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already

	if(get_parameters()->get_num_timescales() == 1) {
		//this is the simplest case, just using 1 timescale
		int steps = get_parameters()->get_integrationsteps1();
		hmc_float stepsize = get_parameters()->get_tau() / ((hmc_float) steps);
		hmc_float stepsize_half = 0.5 * stepsize;
		
	} 
	else if (get_parameters()->get_num_timescales() == 2) {
		int steps1 = get_parameters()->get_integrationsteps1();
		int steps2 = get_parameters()->get_integrationsteps2();
		
		int mult = steps1 % steps2;
		if( mult != 0) Print_Error_Message("integrationsteps1 must be a multiple of integrationssteps2, nothing else is implemented yet. Aborting...");
		
		//this is the number of int.steps for the fermion during one gauge-int.step
		//this is done after hep-lat/0209037. See also hep-lat/050611v2 for a more advanced versions
		int m = steps1/steps2;
		
		//this uses 2 timescales (more is not implemented yet): timescale1 for the gauge-part, timescale2 for the fermion part
		hmc_float stepsize = get_parameters()->get_tau() / ((hmc_float) steps1);
		hmc_float stepsize2 = get_parameters()->get_tau() / ((hmc_float) steps2);
		hmc_float stepsize_half = 0.5 * stepsize;
		hmc_float stepsize2_half = 0.5 * stepsize2;
		
	}
	else 
		Print_Error_Message("More than 2 timescales is not implemented yet. Aborting...");
}

void Gaugefield_hmc::init_gaugemomentum_spinorfield(){
	//init gauge_momenta, saved in clmem_p
	get_task_hmc(0)->generate_gaussian_gaugemomenta_device();
	//init/update spinorfield phi
	get_task_hmc(0)->generate_spinorfield_gaussian();
	get_task_hmc(0)->calc_spinorfield_init_energy();
	get_task_hmc(0)->md_update_spinorfield();

}
