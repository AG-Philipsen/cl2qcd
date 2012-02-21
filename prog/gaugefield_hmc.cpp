#include "gaugefield_hmc.h"


Opencl_Module_Hmc* Gaugefield_hmc::get_task_hmc(int dev)
{
	//@todo: if more than one device is used, the element dev from the array must be called here!!
	return (Opencl_Module_Hmc*)opencl_modules[dev];
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
	this->fermion_forces_call(solvertimer);
	get_task_hmc(0)->calc_gauge_force();
}

void Gaugefield_hmc::md_update_gaugemomentum(hmc_float eps, usetimer * solvertimer){
	get_task_hmc(0)->set_zero_clmem_force_device();
	calc_total_force(solvertimer);
	get_task_hmc(0)->md_update_gaugemomentum_device(-1.*eps);
	return;
}

void Gaugefield_hmc::md_update_gaugemomentum_gauge(hmc_float eps){
	logger.debug() << "update gauge with " << eps;
	get_task_hmc(0)->set_zero_clmem_force_device();
	get_task_hmc(0)->calc_gauge_force();
	get_task_hmc(0)->md_update_gaugemomentum_device(-1.*eps);
	return;
}

void Gaugefield_hmc::md_update_gaugemomentum_fermion(hmc_float eps, usetimer * solvertimer){
	logger.debug() << "update fermion with " << eps;
	get_task_hmc(0)->set_zero_clmem_force_device();
	this->fermion_forces_call(solvertimer);
	get_task_hmc(0)->md_update_gaugemomentum_device(-1.*eps);
	return;
}

void Gaugefield_hmc::fermion_forces_call(usetimer * solvertimer){
	//in case of stout-smearing we need every intermediate field for the force calculation
	//NOTE: if smearing is not used, this is just 0
	int rho_iter = get_parameters()->get_rho_iter();
	//array to save the intermediate fields
	//NOTE: One needs only rho_iter -1 here since the last iteration is saved in gf...
	//NOTE: If the original gf is also needed in the force calculation, one has to add it here
	//	or use the intermediate cl_mem obj gf_unsmeared. This is initialized in the smear_gaugefield function
	cl_mem * smeared_gfs;
	if(rho_iter > 0) smeared_gfs = new cl_mem [rho_iter -1];
	else smeared_gfs = NULL;
		
	if(get_parameters()->get_use_smearing() == true) {
		size_t gfsize = get_parameters()->get_gf_buf_size();
		for(int i = 0; i < rho_iter; i++)
			smeared_gfs[i] = get_task_hmc(0)->create_rw_buffer(gfsize);
		get_task_hmc(0)->smear_gaugefield(*(get_task_hmc(0)->get_gaugefield()), smeared_gfs);
	}
	get_task_hmc(0)->calc_fermion_force(solvertimer);
	if(get_parameters()->get_use_smearing() == true) {
		get_task_hmc(0)->stout_smeared_fermion_force_device(smeared_gfs);
		get_task_hmc(0)->unsmear_gaugefield(*(get_task_hmc(0)->get_gaugefield()));
		for(int i = 0; i < rho_iter; i++){
			cl_int clerr =clReleaseMemObject(smeared_gfs[i]);
			if(clerr != CL_SUCCESS) Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		}
	}
}

void Gaugefield_hmc::md_update_gaugefield(hmc_float eps){
    		logger.debug() << "update gf with " << eps;
	get_task_hmc(0)->md_update_gaugefield_device(eps);
	return;
}

void Gaugefield_hmc::integrator(usetimer * solvertimer){
	//CP: at the moment, one can only use the same type of integrator if one uses more then one timescale...
	if (get_parameters()->get_num_timescales() == 2) {
		if(( get_parameters()->get_integrator(0) != get_parameters()->get_integrator(1)  )){
			logger.fatal() << "Different timescales must use the same integrator up to now!\nAborting...";
			exit(1);
		}
	}
	if(get_parameters()->get_integrator(1) == LEAPFROG){
		this->leapfrog(solvertimer);
	}
	else if(get_parameters()->get_integrator(1) == TWOMN){
		this->twomn(solvertimer);
	}
	return;
}

void Gaugefield_hmc::leapfrog(usetimer * solvertimer)
{
	//it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already when this function is called the first time

	if(get_parameters()->get_num_timescales() == 1) {
		logger.debug() << "\t\tstarting leapfrog...";
		int n0 = get_parameters()->get_integrationsteps0();
		hmc_float deltaTau0 = get_parameters()->get_tau() / ((hmc_float) n0);
		hmc_float deltaTau0_half = 0.5 * deltaTau0;

		logger.debug() << "\t\tinitial step:";
		md_update_gaugemomentum(deltaTau0_half, solvertimer);
		if(n0 > 1) logger.debug() << "\t\tperform " << n0 - 1 << " intermediate steps " ;
		for(int k = 1; k < n0; k++) {
			md_update_gaugefield(deltaTau0);
			md_update_gaugemomentum(deltaTau0, solvertimer);
		}
		logger.debug() << "\t\tfinal step" ;
		md_update_gaugefield(deltaTau0);
		md_update_gaugemomentum(deltaTau0_half, solvertimer);
		logger.debug() << "\t\t...finished leapfrog";
	} 
	else if (get_parameters()->get_num_timescales() == 2) {
		logger.debug() << "start leapfrog with 2 timescales..";
		//this uses 2 timescales (more is not implemented yet): timescale0 for the gauge-part, timescale1 for the fermion part
 		//this is done after hep-lat/0209037. See also hep-lat/0506011v2 for a more advanced version		
		int n0 = get_parameters()->get_integrationsteps0();
		int n1 = get_parameters()->get_integrationsteps1();
		hmc_float deltaTau1 = get_parameters()->get_tau() / ((hmc_float) n1);
		hmc_float deltaTau0 = deltaTau1 / ( (hmc_float) n0 );
		hmc_float deltaTau0_half = 0.5 * deltaTau0;
		hmc_float deltaTau1_half = 0.5 * deltaTau1;

		logger.debug() << "\t\tinitial step:";
		//this corresponds to V_s2(deltaTau/2)
		md_update_gaugemomentum_fermion(deltaTau1_half, solvertimer);
		//now, m steps "more" are performed for the gauge-part
		//this corresponds to [V_s1(deltaTau/2/m) V_t(deltaTau/m) V_s1(deltaTau/2/m) ]^m
		for(int l = 0; l < n0; l++) {
			md_update_gaugemomentum_gauge(deltaTau0_half);
			md_update_gaugefield(deltaTau0);
			md_update_gaugemomentum_gauge(deltaTau0_half);
		}
		if(n1 > 1) logger.debug() << "\t\tperform " << n1 - 1 << " intermediate steps " ;
		for(int k = 1; k < n1; k++) {
			//this corresponds to V_s2(deltaTau)
			md_update_gaugemomentum_fermion(deltaTau1, solvertimer);
			for(int l = 0; l < n0; l++) {
				//this corresponds to [V_s1(deltaTau/2/m) V_t(deltaTau/m) V_s1(deltaTau/2/m) ]^m
				md_update_gaugemomentum_gauge(deltaTau0_half);
				md_update_gaugefield(deltaTau0);
				md_update_gaugemomentum_gauge(deltaTau0_half);
			}
		}
		logger.debug() << "\t\tfinal step" ;
		//this corresponds to the missing V_s2(deltaTau/2)
		md_update_gaugemomentum_fermion(deltaTau1_half, solvertimer);
		logger.debug() << "\t\tfinished leapfrog";
	}
	else 
		Print_Error_Message("More than 2 timescales is not implemented yet. Aborting...");
}

void Gaugefield_hmc::twomn(usetimer * solvertimer)
{
	//it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already
	if(get_parameters()->get_num_timescales() == 1) {
		logger.debug() << "\t\tstarting 2MN...";
		int n0 = get_parameters()->get_integrationsteps0();
		hmc_float deltaTau0 = get_parameters()->get_tau() / ((hmc_float) n0);
		hmc_float deltaTau0_half = 0.5 * deltaTau0;
		hmc_float lambda_times_deltaTau0 = deltaTau0*get_parameters()->get_lambda0();
		hmc_float one_minus_2_lambda = 1. - 2.*get_parameters()->get_lambda0();
		hmc_float one_minus_2_lambda_times_deltaTau0 = one_minus_2_lambda * deltaTau0;

		logger.debug() << "\t\tinitial step:";
		md_update_gaugemomentum(lambda_times_deltaTau0, solvertimer);
		md_update_gaugefield(deltaTau0_half);
		md_update_gaugemomentum(one_minus_2_lambda_times_deltaTau0, solvertimer);
		md_update_gaugefield(deltaTau0_half);
		
		if(n0 > 1) logger.debug() << "\t\tperform " << n0 - 1 << " intermediate steps " ;
		for(int k = 1; k < n0; k++) {
			md_update_gaugemomentum(2.*lambda_times_deltaTau0, solvertimer);
			md_update_gaugefield(deltaTau0_half);
			md_update_gaugemomentum(one_minus_2_lambda_times_deltaTau0, solvertimer);
			md_update_gaugefield(deltaTau0_half);
		}
		logger.debug() << "\t\tfinal step" ;
		md_update_gaugemomentum(lambda_times_deltaTau0, solvertimer);
		logger.debug() << "\t\tfinished 2MN";
	} 
	else if (get_parameters()->get_num_timescales() == 2) {
 		//this is done after hep-lat/0209037. See also hep-lat/0506011v2 for a more advanced version		
		int n0 = get_parameters()->get_integrationsteps0();
		int n1 = get_parameters()->get_integrationsteps1();
		
		//this uses 2 timescales (more is not implemented yet): timescale1 for the gauge-part, timescale2 for the fermion part
		hmc_float deltaTau1 = get_parameters()->get_tau() / ((hmc_float) n1);
		//NOTE: With 2MN, the stepsize for the lower integration step is deltaTau1/(2 N0)!!
		hmc_float deltaTau0= deltaTau1 / ( 2.* (hmc_float) n0 );
		hmc_float deltaTau0_half = 0.5 * deltaTau0;
		//hmc_float deltaTau1_half = 0.5 * deltaTau1;
		hmc_float lambda0_times_deltaTau0= deltaTau0*get_parameters()->get_lambda0();
		hmc_float lambda1_times_deltaTau1 = deltaTau1*get_parameters()->get_lambda1();
		hmc_float one_minus_2_lambda0 = 1. - 2.*get_parameters()->get_lambda0();
		hmc_float one_minus_2_lambda1 = 1. - 2.*get_parameters()->get_lambda1();
		hmc_float one_minus_2_lambda0_times_deltaTau0 = one_minus_2_lambda0 * deltaTau0;
		hmc_float one_minus_2_lambda1_times_deltaTau1 = one_minus_2_lambda1 * deltaTau1;

		logger.debug() << "\t\tinitial step:";
		//this corresponds to V_s2(lambda*deltaTau)
		md_update_gaugemomentum_fermion(lambda1_times_deltaTau1, solvertimer);
		//now, n0 steps "more" are performed for the gauge-part
		//this corresponds to [exp(lambda*eps T(V_gauge) ) exp( eps/2 V ) exp( (1 - 2lamdba) *eps T(V_gauge) ) exp( eps/2 V ) exp( lamdba*eps T(V_gauge) ) ]^m
		for(int l = 0; l < n0; l++) {
			md_update_gaugemomentum_gauge(lambda0_times_deltaTau0);
			md_update_gaugefield(deltaTau0_half);
			md_update_gaugemomentum_gauge(one_minus_2_lambda0_times_deltaTau0);
			md_update_gaugefield(deltaTau0_half);
			md_update_gaugemomentum_gauge(lambda0_times_deltaTau0);
		}
		//this corresponds to V_s2( ( 1 - 2lambda) *deltaTau)
		md_update_gaugemomentum_fermion(one_minus_2_lambda1_times_deltaTau1, solvertimer);
		//now, m steps "more" are performed for the gauge-part (again)
		for(int l = 0; l < n0; l++) {
			md_update_gaugemomentum_gauge(lambda0_times_deltaTau0);
			md_update_gaugefield(deltaTau0_half);
			md_update_gaugemomentum_gauge(one_minus_2_lambda0_times_deltaTau0);
			md_update_gaugefield(deltaTau0_half);
			md_update_gaugemomentum_gauge(lambda0_times_deltaTau0);
		}
		//the last V_s2(lambda*deltaTau) can be pulled into the intermediate steps
		if(n1 > 1) logger.debug() << "\t\tperform " << n1 - 1 << " intermediate steps " ;
		for(int k = 1; k < n1; k++) {
			//this corresponds to V_s2(deltaTau)
			md_update_gaugemomentum_fermion(2.*lambda1_times_deltaTau1, solvertimer);
			for(int l = 0; l < n0; l++) {
				md_update_gaugemomentum_gauge(lambda0_times_deltaTau0);
				md_update_gaugefield(deltaTau0_half);
				md_update_gaugemomentum_gauge(one_minus_2_lambda0_times_deltaTau0);
				md_update_gaugefield(deltaTau0_half);
				md_update_gaugemomentum_gauge(lambda0_times_deltaTau0);
			}
			md_update_gaugemomentum_fermion(one_minus_2_lambda1_times_deltaTau1, solvertimer);
			for(int l = 0; l < n0; l++) {
				md_update_gaugemomentum_gauge(lambda0_times_deltaTau0);
				md_update_gaugefield(deltaTau0_half);
				md_update_gaugemomentum_gauge(one_minus_2_lambda0_times_deltaTau0);
				md_update_gaugefield(deltaTau0_half);
				md_update_gaugemomentum_gauge(lambda0_times_deltaTau0);
			}
		}
		logger.debug() << "\t\tfinal step" ;
		//this corresponds to the missing V_s2(lambda*deltaTau)
		md_update_gaugemomentum_fermion(lambda1_times_deltaTau1, solvertimer);
		logger.debug() << "\t\tfinished 2MN";
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
