#include "opencl_hmc.h"
#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

hmc_error Opencl_hmc::fill_collect_options(stringstream* collect_options)
{

	Opencl_fermions::fill_collect_options(collect_options);
	*collect_options <<  " -DBETA=" << get_parameters()->get_beta() << " -DGAUGEMOMENTASIZE=" << get_parameters()->get_gaugemomentasize();
	return HMC_SUCCESS;	
}

hmc_error Opencl_hmc::fill_buffers()
{
	Opencl_fermions::fill_buffers();

	size_t local_work_size;
	size_t global_work_size;
	cl_uint num_groups;
	//CP: This has no effect yet!!
	char * kernelname = "dummy";
	get_work_sizes(&local_work_size, &global_work_size, &num_groups, Opencl::get_device_type(), kernelname);	


	int spinorfield_size = sizeof(spinor) * SPINORFIELDSIZE;
// 	int eoprec_spinorfield_size = sizeof(spinor) * EOPREC_SPINORFIELDSIZE;
	int gaugemomentum_size = sizeof(ae) * GAUGEMOMENTASIZE2;
// 	int complex_size = sizeof(hmc_complex);
	int float_size = sizeof(hmc_float);
// 	int global_buf_size = complex_size * num_groups;
// 	int global_buf_size_float = float_size * num_groups;
	hmc_complex one = hmc_complex_one;
	hmc_complex minusone = hmc_complex_minusone;
// 	hmc_float tmp;

	//init mem-objects

	logger.trace() << "Create buffer for HMC...";
	clmem_force = create_rw_buffer(gaugemomentum_size);
	clmem_phi_inv = create_rw_buffer(spinorfield_size);
	clmem_phi = create_rw_buffer(spinorfield_size);
	clmem_new_u = create_rw_buffer(sizeof(s_gaugefield));
	clmem_p = create_rw_buffer(gaugemomentum_size);
	clmem_new_p = create_rw_buffer(gaugemomentum_size);
	clmem_energy_init = create_rw_buffer(float_size);
	clmem_p2 = create_rw_buffer(float_size);
	clmem_new_p2 = create_rw_buffer(float_size);
	clmem_s_fermion = create_rw_buffer(float_size);
	
	return HMC_SUCCESS;
}

void Opencl_hmc::fill_kernels()
{
	//fill kernels of Mother classes
	Opencl_fermions::fill_kernels();

	basic_hmc_code = basic_fermion_code << "types_hmc.h";

	//init kernels for HMC
	set_zero_gaugemomentum = createKernel("set_zero_gaugemomentum") << basic_hmc_code << "gaugemomentum_zero.cl";
	generate_gaussian_spinorfield = createKernel("generate_gaussian_spinorfield") << basic_hmc_code << "random.cl" << "spinorfield_gaussian.cl";
	generate_gaussian_gaugemomenta = createKernel("generate_gaussian_gaugemomenta") << basic_hmc_code << "random.cl" << "gaugemomentum_gaussian.cl";
	md_update_gaugefield = createKernel("md_update_gaugefield") << basic_hmc_code << "md_update_gaugefield.cl";
	md_update_gaugemomenta = createKernel("md_update_gaugemomenta") << basic_hmc_code << "gaugemomentum.cl" << "md_update_gaugemomenta.cl";
	gauge_force = createKernel("gauge_force") << basic_hmc_code << "gaugemomentum.cl" << "force_gauge.cl";
	fermion_force = createKernel("fermion_force") << basic_hmc_code << "gaugemomentum.cl" << "fermionmatrix.cl" << "force_fermion.cl";
	if(get_parameters()->get_use_smearing() == true){
		stout_smear_fermion_force = createKernel("stout_smear_fermion_force") << basic_fermion_code << "stout_smear_fermion_force.cl";
	}
	gaugemomentum_squarenorm = createKernel("gaugemomentum_squarenorm") << basic_hmc_code << "gaugemomentum_squarenorm.cl";
}

hmc_error Opencl_hmc::init(cl_device_type wanted_device_type, inputparameters* parameters)
{
	hmc_error err = Opencl_fermions::init(wanted_device_type, parameters);
	
	return err;
}

hmc_error Opencl_hmc::finalize_hmc()
{

	logger.debug() << "release HMC-variables.." ;
	if(clReleaseMemObject(clmem_energy_init) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_p2) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_new_p2) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_p) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_new_p) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_new_u) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_phi_inv) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_force) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;

	logger.debug() << "release HMC-kernels.." ;
	if(clReleaseKernel(generate_gaussian_spinorfield) != CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(generate_gaussian_gaugemomenta) != CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(md_update_gaugefield) != CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(md_update_gaugemomenta) != CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(gauge_force) != CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(fermion_force) != CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(set_zero_gaugemomentum) != CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(get_parameters()->get_use_smearing() == true){
		if(clReleaseKernel(stout_smear_fermion_force) != CL_SUCCESS) return HMC_RELEASEKERNELERR;
	}
	return HMC_SUCCESS;
}

////////////////////////////////////////////////////
//Methods needed for the HMC-algorithm

hmc_error Opencl_hmc::generate_gaussian_gaugemomenta_device(const size_t ls, const size_t gs)
{
	int clerr;
	//this is always applied to clmem_new_p
	clerr = clSetKernelArg(generate_gaussian_gaugemomenta, 0, sizeof(cl_mem), &clmem_p);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(generate_gaussian_gaugemomenta, 1, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( generate_gaussian_gaugemomenta , gs, ls);

	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::generate_gaussian_spinorfield_device(const size_t ls, const size_t gs)
{
	int clerr;
	//this is always applied to clmem_phi_inv, which can be done since the gaussian field is only needed in the beginning
	clerr = clSetKernelArg(generate_gaussian_spinorfield, 0, sizeof(cl_mem), &clmem_phi_inv);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(generate_gaussian_spinorfield, 1, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}

	enqueueKernel(generate_gaussian_spinorfield  , gs, ls);

	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::md_update_spinorfield_device(const size_t local_work_size, const size_t global_work_size)
{
	//suppose the initial gaussian field is saved in clmem_phi_inv (see above).
	//	then the "phi" = Dpsi from the algorithm is stored in clmem_phi
	//	which then has to be the source of the inversion
	int err =  Opencl_fermions::Qplus(clmem_phi_inv, clmem_phi , get_clmem_gaugefield(), local_work_size, global_work_size);

	usetimer noop;
	hmc_float s_fermion;
	set_float_to_global_squarenorm_device(clmem_phi, clmem_s_fermion, local_work_size, global_work_size);
	get_buffer_from_device(clmem_s_fermion, &s_fermion, sizeof(hmc_float), &noop);
	logger.debug() << "\tsquarenorm init field after update = " << s_fermion;
		
	if(err != HMC_SUCCESS){
		logger.fatal() << "error occured in md_update_spinorfield_device.. ";
		exit (HMC_STDERR);
	}
	return HMC_SUCCESS;
}


//CP: steps2 is not used at the moment....
hmc_error Opencl_hmc::leapfrog_device(hmc_float tau, int steps1, int steps2, usetimer *copy_to, usetimer * copy_on, usetimer * solvertimer, const size_t ls, const size_t gs)
{
	//it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already
	//this is the number of int.steps for the fermion during one gauge-int.step
	//this is done after hep-lat/0209037. See also hep-lat/050611v2 for a more advanced versions
	int mult = steps1%steps2;
	if( mult != 0){
		logger.fatal() << "integrationsteps1 must be a multiple of integrationssteps2, nothing else is implemented yet. Aborting...";
		exit(HMC_STDERR);
	}
	int m = steps1/steps2;
	if(m == 1){
		//this is the simplest case, just using 1 timescale
		hmc_float stepsize = (tau) / ((hmc_float) steps1);
		hmc_float stepsize_half = 0.5 * stepsize;
		
		//initial step
		logger.debug() << "\t\tinitial step:";
		calc_total_force(copy_to, copy_on, solvertimer, ls, gs);
		md_update_gaugemomentum_device(-1.*stepsize_half, ls, gs);
		//intermediate steps
		if(steps1 > 1) logger.debug() << "\t\tperform " << steps1 - 1 << " intermediate steps " ;
		for(int k = 1; k < steps1; k++) {
			md_update_gaugefield_device(stepsize, ls, gs);
			calc_total_force(copy_to, copy_on, solvertimer, ls, gs);
			md_update_gaugemomentum_device(-1.*stepsize, ls, gs);
		}
		//final step
		logger.debug() << "\t\tfinal step" ;
		md_update_gaugefield_device(stepsize, ls, gs);
		calc_total_force(copy_to, copy_on, solvertimer, ls, gs);
		md_update_gaugemomentum_device(-1.*stepsize_half, ls, gs);
		logger.debug() << "\t\tfinished leapfrog";
	}
	else {
		//this uses 2 timescales (more is not implemented yet): timescale1 for the gauge-part, timescale2 for the fermion part
		hmc_float stepsize = (tau) / ((hmc_float) steps1);
		hmc_float stepsize2 = (tau) / ((hmc_float) steps2);
		hmc_float stepsize_half = 0.5 * stepsize;
		hmc_float stepsize2_half = 0.5 * stepsize2;
		
		//initial step
		logger.debug() << "\t\tinitial step:";
		//this corresponds to V_s2(deltaTau/2)
		set_zero_clmem_force_device(ls, gs);
		calc_fermion_force(copy_to, copy_on, solvertimer, ls, gs);
		md_update_gaugemomentum_device(-1.*stepsize2_half, ls, gs);
		//now, m steps "more" a performed for the gauge-part
		//this corresponds to [V_s1(deltaTau/2/m) V_t(deltaTau/m) V_s1(deltaTau/2/m) ]^m
		for(int l=0; l<m; l++){
			set_zero_clmem_force_device(ls, gs);
			calc_gauge_force(ls, gs);
			md_update_gaugemomentum_device(-1.*stepsize_half, ls, gs);
			md_update_gaugefield_device(stepsize2, ls, gs);
			set_zero_clmem_force_device(ls, gs);
			calc_gauge_force(ls, gs);
			md_update_gaugemomentum_device(-1.*stepsize_half, ls, gs);
		}
		//intermediate steps
		if(steps1 > 1) logger.debug() << "\t\tperform " << steps2 - 1 << " intermediate steps " ;
		for(int k = 1; k < steps2; k++) {
			//this corresponds to V_s2(deltaTau)
			set_zero_clmem_force_device(ls, gs);
			calc_fermion_force(copy_to, copy_on, solvertimer, ls, gs);
			md_update_gaugemomentum_device(-1.*stepsize2, ls, gs);
			for(int l=0; l<m; l++){
				//this corresponds to [V_s1(deltaTau/2/m) V_t(deltaTau/m) V_s1(deltaTau/2/m) ]^m
				set_zero_clmem_force_device(ls, gs);
				calc_gauge_force(ls, gs);
				md_update_gaugemomentum_device(-1.*stepsize_half, ls, gs);
				md_update_gaugefield_device(stepsize2, ls, gs);
				set_zero_clmem_force_device(ls, gs);
				calc_gauge_force(ls, gs);
				md_update_gaugemomentum_device(-1.*stepsize_half, ls, gs);
			}
		}
		//final step
		logger.debug() << "\t\tfinal step" ;
		//this corresponds to the missing V_s2(deltaTau/2)
		set_zero_clmem_force_device(ls, gs);
		calc_fermion_force(copy_to, copy_on, solvertimer, ls, gs);
		md_update_gaugemomentum_device(-1.*stepsize2_half, ls, gs);
		logger.debug() << "\t\tfinished leapfrog";
	}
	
	return HMC_SUCCESS;
}

void Opencl_hmc::calc_total_force(usetimer *copy_to, usetimer * copy_on, usetimer * solvertimer, const size_t ls, const size_t gs){
	//CP: make sure that the output field is set to zero
	set_zero_clmem_force_device(ls, gs);
	this->calc_fermion_force(copy_to, copy_on, solvertimer, ls, gs);
	this->calc_gauge_force(ls, gs);
}

void Opencl_hmc::calc_fermion_force(usetimer *copy_to, usetimer * copy_on, usetimer * solvertimer, const size_t ls, const size_t gs){
	int err = HMC_SUCCESS;
	//this is only used when smearing is activated
	cl_mem gf_tmp;
	
	if(get_parameters()->get_use_smearing() == true){
		logger.debug() << "\t\t\tsave unsmeared gaugefield...";
		gf_tmp = create_rw_buffer(sizeof(s_gaugefield));
		copy_buffer_on_device(get_clmem_gaugefield(), gf_tmp, sizeof(s_gaugefield), copy_on);
		logger.debug() << "\t\t\tsmear gaugefield...";
		stout_smear_device(ls, gs);
	}
	
	//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
	if(get_parameters()->get_use_eo() == true){
		if(get_parameters()->get_use_cg() == true){
			//this is broken right now since the CG doesnt work!!
			logger.debug() << "\t\tcalc fermion force ingredients using cg and eoprec";
			logger.fatal() << "this is not implemented yet. Aborting..";
			exit(HMC_STDERR); 
		}
		else{
			logger.debug() << "\t\tcalc fermion force ingredients using bicgstab and eoprec";
			logger.fatal() << "this is not implemented yet. Aborting..";
			exit(HMC_STDERR);
		}
	}
	else{
		if(get_parameters()->get_use_cg() == true){
			//this is broken right now since the CG doesnt work!!
			logger.debug() << "\t\tcalc fermion force ingredients using cg and no eoprec";
			logger.fatal() << "this is not implemented yet. Aborting..";
			exit(HMC_STDERR);
		}
		else{
			logger.debug() << "\t\tcalc fermion force ingredients using bicgstab and no eoprec";
			/**
			* The first inversion calculates 
			* Y = phi = (Qplus)^-1 psi
			* out of 
			* Qplus phi = psi
			* This is also the energy of the final field!
			*/
			logger.debug() << "\t\t\tstart solver";
			
			/** @todo at the moment, we can only put in a cold spinorfield 
				* or a point-source spinorfield as trial-solution 
				*/
			/**
			  * Trial solution for the spinorfield
	 		  */
			set_spinorfield_cold_device(get_clmem_inout(), ls, gs);			
			
			//debugging
			hmc_float s_fermion;
			set_float_to_global_squarenorm_device(get_clmem_inout(), clmem_s_fermion, local_work_size, global_work_size);
			get_buffer_from_device(clmem_s_fermion, &s_fermion, sizeof(hmc_float), copy_to);
			logger.debug() << "\tsquarenorm of inv.field before = " << s_fermion;
			 
			err = Opencl_fermions::solver_device(Qplus_call, get_clmem_inout(), get_clmem_phi(), clmem_new_u, copy_to, copy_on, solvertimer, ls, gs, get_parameters()->get_cgmax());
			if (err != HMC_SUCCESS) logger.debug() << "\t\t\tsolver did not solve!!";
			else logger.debug() << "\t\t\tsolver solved!";
			
			//debugging
			set_float_to_global_squarenorm_device(get_clmem_inout(), clmem_s_fermion, local_work_size, global_work_size);
			get_buffer_from_device(clmem_s_fermion, &s_fermion, sizeof(hmc_float), copy_to);
			logger.debug() << "\tsquarenorm of inv.field after = " << s_fermion;
			
			//store this result in clmem_phi_inv
			copy_buffer_on_device(get_clmem_inout(), clmem_phi_inv, sizeof(spinor) * SPINORFIELDSIZE, copy_on);
			
			/**
			 * Now, one has to calculate 
			 * X = (Qminus)^-1 Y = (Qminus)^-1 (Qplus)^-1 psi = (QplusQminus)^-1 psi ??
			 * out of 
			 * Qminus clmem_inout = clmem_phi_inv
			 * Therefore, when calculating the final energy of the spinorfield,
			 * 	one can just take clmem_phi_inv (see also above)!!
			 */
			
			//copy former solution to clmem_source
			copy_buffer_on_device(get_clmem_inout(), get_clmem_source(), sizeof(spinor) * SPINORFIELDSIZE, copy_on);
			logger.debug() << "\t\t\tstart solver";
			
			//debugging
			set_float_to_global_squarenorm_device(get_clmem_inout(), clmem_s_fermion, local_work_size, global_work_size);
			get_buffer_from_device(clmem_s_fermion, &s_fermion, sizeof(hmc_float), copy_to);
			logger.debug() << "\tsquarenorm of inv.field before = " << s_fermion;
			 
			//this sets clmem_inout cold as trial-solution
			set_spinorfield_cold_device(get_clmem_inout(), ls, gs);
			
			err = Opencl_fermions::solver_device(Qminus_call, get_clmem_inout(), get_clmem_source(), clmem_new_u, copy_to, copy_on, solvertimer, ls, gs, get_parameters()->get_cgmax());
			if (err != HMC_SUCCESS) logger.debug() << "\t\t\tsolver did not solve!!";
			else logger.debug() << "\t\t\tsolver solved!";
			
			//debugging
			set_float_to_global_squarenorm_device(get_clmem_inout(), clmem_s_fermion, local_work_size, global_work_size);
			get_buffer_from_device(clmem_s_fermion, &s_fermion, sizeof(hmc_float), copy_to);
			logger.debug() << "\tsquarenorm of inv.field after = " << s_fermion;
		}	
	}
	logger.debug() << "\t\tcalc fermion_force...";
	//CP: this always calls fermion_force(Y,X) with Y = clmem_phi_inv, X = clmem_inout
	fermion_force_device(ls, gs);
	
	if(get_parameters()->get_use_smearing() == true){
		//CP: The fermion_force call above used the smeared gaugefield if wished,
		//	now the force by the thin link has to be determined...
		logger.debug() << "\t\t\tcalc stout-smeared fermion_force...";
		stout_smeared_fermion_force_device(ls, gs);
		logger.debug() << "\t\t\trestore unsmeared gaugefield...";
		copy_buffer_on_device(gf_tmp, get_clmem_gaugefield(), sizeof(s_gaugefield), copy_on);
		if(clReleaseMemObject(gf_tmp) != CL_SUCCESS) exit(HMC_OCLERROR);
	}	
}

void Opencl_hmc::calc_gauge_force(const size_t ls, const size_t gs){
	logger.debug() << "\t\tcalc gauge_force...";
	gauge_force_device(ls, gs);
}

hmc_observables Opencl_hmc::metropolis(hmc_float rnd, hmc_float beta, const size_t local_work_size, const size_t global_work_size, usetimer * timer)
{
	//Calc Hamiltonian
	logger.debug() << "Calculate Hamiltonian";
	
	//Gauge-Part
	hmc_float tplaq, splaq, plaq;
	hmc_float tplaq_new, splaq_new, plaq_new;
	hmc_complex poly;
	hmc_complex poly_new;
	//In this call, the observables are calculated already with appropiate Weighting factor of 2.0/(VOL4D*NDIM*(NDIM-1)*NC)
	Opencl::gaugeobservables(get_clmem_gaugefield(), &plaq,  &tplaq, &splaq, &poly);
	Opencl::gaugeobservables(clmem_new_u, &plaq_new,  &tplaq_new, &splaq_new, &poly_new);
	//plaq has to be divided by the norm-factor and multiplied by NC 
	//	(because this is in the defintion of the gauge action and not in the normalization) to get s_gauge
	hmc_float factor = 2.0 / static_cast<hmc_float>(VOL4D * NDIM * (NDIM - 1) );
	/** NOTE: the minus here is introduced to fit tmlqcd!!! */
	hmc_float deltaH = -(plaq - plaq_new) * beta / factor;
	
	logger.debug() << "\tS_gauge(old field) = " << setprecision(10) <<plaq << "\t" << plaq* beta  /factor;
	logger.debug() << "\tS_gauge(new field) = " << setprecision(10) <<plaq_new << "\t" << plaq_new* beta /factor;
	logger.debug() << "\tdeltaS_gauge = " << setprecision(10) << deltaH;
	
	//Gaugemomentum-Part
	hmc_float p2, new_p2;
	set_float_to_gaugemomentum_squarenorm_device(clmem_p, clmem_p2, local_work_size, global_work_size);
	set_float_to_gaugemomentum_squarenorm_device(clmem_new_p, clmem_new_p2, local_work_size, global_work_size);
	Opencl_fermions::get_buffer_from_device(clmem_p2, &p2, sizeof(hmc_float), timer);

	Opencl_fermions::get_buffer_from_device(clmem_new_p2, &new_p2, sizeof(hmc_float), timer);
	//the energy is half the squarenorm
	deltaH += 0.5 * (p2 - new_p2);
	
	logger.debug() << "\tS_gaugemom(old field) = " << setprecision(10) << 0.5*p2;
	logger.debug() << "\tS_gaugemom(new field) = " << setprecision(10) <<0.5*new_p2;
	logger.debug() << "\tdeltaS_gaugemom = " << setprecision(10) <<0.5 * (p2 - new_p2);

	//Fermion-Part:
	hmc_float spinor_energy_init, s_fermion;
	//initial energy has been computed in the beginning...
	Opencl_fermions::get_buffer_from_device(clmem_energy_init, &spinor_energy_init, sizeof(hmc_float), timer);
	// sum_links phi*_i (M^+M)_ij^-1 phi_j
	// the inversion with respect to the input-gaussian field and the new gaugefield has taken place during the leapfrog
	// 	and was saved in clmem_phi_inv during that step.
	set_float_to_global_squarenorm_device(clmem_phi_inv, clmem_s_fermion, local_work_size, global_work_size);
	get_buffer_from_device(clmem_s_fermion, &s_fermion, sizeof(hmc_float), timer);
	deltaH += spinor_energy_init - s_fermion;
 
 	logger.debug() << "\tS_ferm(old field) = " << setprecision(10) <<  spinor_energy_init;
	logger.debug() << "\tS_ferm(new field) = " << setprecision(10)<< s_fermion;
	logger.debug() << "\tdeltaS_ferm = " << spinor_energy_init - s_fermion;

	//Metropolis-Part
	hmc_float compare_prob;
	if(deltaH < 0) {
		compare_prob = exp(deltaH);
	} else {
		compare_prob = 1.0;
	}
	logger.debug() << "\tdeltaH = " << deltaH << "\tAcc-Prop = " << compare_prob;
	hmc_observables tmp;
	if(rnd <= compare_prob) {
		tmp.accept = 1;
		tmp.plaq = plaq_new;
		tmp.tplaq = tplaq_new;
		tmp.splaq = splaq_new;
		tmp.poly = poly_new;
		tmp.deltaH = deltaH;
		tmp.prob = compare_prob;
	} else {
		tmp.accept = 0;
		tmp.plaq = plaq;
		tmp.tplaq = tplaq;
		tmp.splaq = splaq;
		tmp.poly = poly;
		tmp.deltaH = deltaH;
		tmp.prob = compare_prob;
	}

	return tmp;
}

hmc_error Opencl_hmc::calc_spinorfield_init_energy_device(const size_t local_work_size, const size_t global_work_size)
{
	//Suppose the initial spinorfield is saved in phi_inv
	Opencl_fermions::set_float_to_global_squarenorm_device(clmem_phi_inv, clmem_energy_init, local_work_size, global_work_size);

	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::md_update_gaugemomentum_device(hmc_float eps, const size_t ls, const size_t gs)
{
	//__kernel void md_update_gaugemomenta(hmc_float eps, __global ae * p_inout, __global ae* force_in){

	hmc_float tmp = eps;

	int clerr;
	clerr = clSetKernelArg(md_update_gaugemomenta, 0, sizeof(hmc_float), &tmp);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(md_update_gaugemomenta, 1, sizeof(cl_mem), &clmem_new_p);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(md_update_gaugemomenta, 2, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	
	enqueueKernel(md_update_gaugemomenta  , gs, ls);
	
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::md_update_gaugefield_device(hmc_float eps, const size_t ls, const size_t gs)
{
	// __kernel void md_update_gaugefield(hmc_float eps, __global ae * p_in, __global ocl_s_gaugefield * u_inout){
	hmc_float tmp = eps;
	int clerr;
	//this is always applied to clmem_force
	clerr = clSetKernelArg(md_update_gaugefield, 0, sizeof(hmc_float), &tmp);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(md_update_gaugefield, 1, sizeof(cl_mem), &clmem_new_p);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(md_update_gaugefield, 2, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	
	enqueueKernel( md_update_gaugefield , gs, ls);
	
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::set_zero_clmem_force_device(const size_t ls, const size_t gs)
{
	int clerr;
	//this is always applied to clmem_force
	clerr = clSetKernelArg(set_zero_gaugemomentum, 0, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	
	enqueueKernel( set_zero_gaugemomentum , gs, ls);

	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::gauge_force_device(const size_t ls, const size_t gs)
{
	int clerr;
	clerr = clSetKernelArg(gauge_force, 0, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(gauge_force, 1, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	
	enqueueKernel( gauge_force , gs, ls);
	
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::fermion_force_device(const size_t ls, const size_t gs)
{
	//fermion_force(field, Y, X, out);
	cl_mem tmp = get_clmem_inout();
	int clerr;
	clerr = clSetKernelArg(fermion_force, 0, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(fermion_force, 1, sizeof(cl_mem), &clmem_phi_inv);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(fermion_force, 2, sizeof(cl_mem), &tmp);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(fermion_force, 3, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	
	enqueueKernel( fermion_force , gs, ls);
	
	return HMC_SUCCESS;

}

hmc_error Opencl_hmc::stout_smeared_fermion_force_device(const size_t ls, const size_t gs)
{
		
	return HMC_SUCCESS;

}

hmc_error Opencl_hmc::set_float_to_gaugemomentum_squarenorm_device(cl_mem clmem_in, cl_mem clmem_out, const size_t ls, const size_t gs)
{
	int clerr = CL_SUCCESS;
	//__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out){
	clerr = clSetKernelArg(gaugemomentum_squarenorm, 0, sizeof(cl_mem), &clmem_in);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
//  /** @todo add reduction */
	clerr = clSetKernelArg(gaugemomentum_squarenorm, 1, sizeof(cl_mem), &clmem_out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	
	enqueueKernel(gaugemomentum_squarenorm  , gs, ls);
	
	return HMC_SUCCESS;
}

////////////////////////////////////////////////////
//Access to members

cl_mem Opencl_hmc::get_clmem_p(){
	return clmem_p;
}

cl_mem Opencl_hmc::get_clmem_new_p(){
	return clmem_new_p;
}

cl_mem Opencl_hmc::get_clmem_new_u(){
	return clmem_new_u;
}

cl_mem Opencl_hmc::get_clmem_phi(){
	return clmem_phi;
}

#ifdef _PROFILING_
usetimer* Opencl_hmc::get_timer(char * in){
	usetimer *noop = NULL;
	noop = Opencl_fermions::get_timer(in);
	if(noop != NULL) return noop;

	if (strcmp(in, "generate_gaussian_spinorfield") == 0){
    return &this->timer_generate_gaussian_spinorfield;
	}	
	if (strcmp(in, "generate_gaussian_gaugemomenta") == 0){
    return &this->timer_generate_gaussian_gaugemomenta;
	}	
	if (strcmp(in, "md_update_gaugefield") == 0){
    return &this->timer_md_update_gaugefield;
	}
	if (strcmp(in, "md_update_gaugemomenta") == 0){
    return &this->timer_md_update_gaugemomenta;
	}
	if (strcmp(in, "gauge_force") == 0){
    return &this->timer_gauge_force;
	}
	if (strcmp(in, "fermion_force") == 0){
    return &this->timer_fermion_force;
	}
	if (strcmp(in, "set_zero_gaugemomentum") == 0){
    return &this->timer_set_zero_gaugemomentum;
	}	
	if (strcmp(in, "gaugemomentum_squarenorm") == 0){
    return &this->timer_gaugemomentum_squarenorm;
	}
	if (strcmp(in, "stout_smear_fermion_force") == 0){
		return &this->timer_stout_smear_fermion_force;
	}
	//if the kernelname has not matched, return NULL
	else{
		return NULL;
	}
}

int Opencl_hmc::get_read_write_size(char * in, inputparameters * parameters){
	Opencl::get_read_write_size(in, parameters);
		//Depending on the compile-options, one has different sizes...
	int D = (*parameters).get_float_size();
	int R = (*parameters).get_mat_size();
	int S;
	if((*parameters).get_use_eo() == 1)
	  S = EOPREC_SPINORFIELDSIZE;
	else
	  S = SPINORFIELDSIZE;
	//this is the same as in the function above
	if (strcmp(in, "generate_gaussian_spinorfield") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "generate_gaussian_gaugemomenta") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "md_update_gaugefield") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "md_update_gaugemomenta") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "gauge_force") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "fermion_force") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "set_zero_gaugemomentum;") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "gaugemomentum_squarenorm") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "stout_smear_fermion_force") == 0){
    return 10000000000000000000;
	}	
	return 0;	
}

void Opencl_hmc::print_profiling(std::string filename){
	Opencl::print_profiling(filename);
	char * kernelName;
	kernelName = "generate_gaussian_spinorfield";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "generate_gaussian_gaugemomenta";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "md_update_gaugefield";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "md_update_gaugemomenta";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "gauge_force";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "fermion_force";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "set_zero_gaugemomentum";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "gaugemomentum_squarenorm";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "stout_smear_fermion_force";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	
}
#endif
