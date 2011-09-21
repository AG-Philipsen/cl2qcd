#include "opencl_module_hmc.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

void Opencl_Module_Hmc::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module_Fermions::fill_collect_options(collect_options);
	*collect_options <<  " -DBETA=" << get_parameters()->get_beta() << " -DGAUGEMOMENTASIZE=" << get_parameters()->get_gaugemomentasize();
	return;
}


void Opencl_Module_Hmc::fill_buffers()
{

	Opencl_Module_Fermions::fill_buffers();

	int spinorfield_size = get_parameters()->get_sf_buf_size();
	int gaugemomentum_size = get_parameters()->get_gm_buf_size();
	int gaugefield_size = get_parameters()->get_gf_buf_size();
	int float_size = sizeof(hmc_float);
	hmc_complex one = hmc_complex_one;
	hmc_complex minusone = hmc_complex_minusone;

	//init mem-objects

	logger.trace() << "Create buffer for HMC...";
	clmem_force = create_rw_buffer(gaugemomentum_size);
	clmem_phi_inv = create_rw_buffer(spinorfield_size);
	clmem_phi = create_rw_buffer(spinorfield_size);
	clmem_new_u = create_rw_buffer(gaugefield_size);
	clmem_p = create_rw_buffer(gaugemomentum_size);
	clmem_new_p = create_rw_buffer(gaugemomentum_size);
	clmem_energy_init = create_rw_buffer(float_size);
	clmem_p2 = create_rw_buffer(float_size);
	clmem_new_p2 = create_rw_buffer(float_size);
	clmem_s_fermion = create_rw_buffer(float_size);

	return;
}

void Opencl_Module_Hmc::fill_kernels()
{
	Opencl_Module_Fermions::fill_kernels();

	basic_hmc_code = basic_fermion_code << "types_hmc.h";

	//init kernels for HMC
	set_zero_gaugemomentum = createKernel("set_zero_gaugemomentum") << basic_hmc_code << "gaugemomentum_zero.cl";
	generate_gaussian_spinorfield = createKernel("generate_gaussian_spinorfield") << basic_hmc_code << "random.cl" << "spinorfield_gaussian.cl";
	generate_gaussian_gaugemomenta = createKernel("generate_gaussian_gaugemomenta") << basic_hmc_code << "random.cl" << "gaugemomentum_gaussian.cl";
	md_update_gaugefield = createKernel("md_update_gaugefield") << basic_hmc_code << "md_update_gaugefield.cl";
	md_update_gaugemomenta = createKernel("md_update_gaugemomenta") << basic_hmc_code << "operations_gaugemomentum.cl" << "md_update_gaugemomenta.cl";
	gauge_force = createKernel("gauge_force") << basic_hmc_code << "operations_gaugemomentum.cl" << "force_gauge.cl";
	fermion_force = createKernel("fermion_force") << basic_hmc_code << "operations_gaugemomentum.cl" << "fermionmatrix.cl" << "force_fermion.cl";
	if(get_parameters()->get_use_smearing() == true) {
		stout_smear_fermion_force = createKernel("stout_smear_fermion_force") << basic_fermion_code << "stout_smear_fermion_force.cl";
	}
	gaugemomentum_squarenorm = createKernel("gaugemomentum_squarenorm") << basic_hmc_code << "gaugemomentum_squarenorm.cl";

	return;
}

void Opencl_Module_Hmc::clear_kernels()
{
	Opencl_Module_Fermions::clear_kernels();

	cl_uint clerr = CL_SUCCESS;

	logger.debug() << "release HMC-kernels.." ;
	clerr = clReleaseKernel(generate_gaussian_spinorfield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(generate_gaussian_gaugemomenta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(md_update_gaugefield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(md_update_gaugemomenta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(gauge_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(fermion_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(set_zero_gaugemomentum);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(get_parameters()->get_use_smearing() == true) {
		clerr = clReleaseKernel(stout_smear_fermion_force);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}

	return;
}

void Opencl_Module_Hmc::clear_buffers()
{
	Opencl_Module_Fermions::clear_buffers();

	cl_uint clerr = CL_SUCCESS;

	logger.debug() << "release HMC-variables.." ;
	clerr = clReleaseMemObject(clmem_energy_init);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_p2);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_new_p2);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_p);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_new_p);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_new_u);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_phi_inv);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);

	return;
}

void Opencl_Module_Hmc::get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups)
{
	Opencl_Module_Fermions::get_work_sizes(kernel, dev_type, ls, gs, num_groups);

	return;
}

////////////////////////////////////////////////////
//Access to members

cl_mem Opencl_Module_Hmc::get_clmem_p()
{
	return clmem_p;
}

cl_mem Opencl_Module_Hmc::get_clmem_new_p()
{
	return clmem_new_p;
}

cl_mem Opencl_Module_Hmc::get_clmem_new_u()
{
	return clmem_new_u;
}

cl_mem Opencl_Module_Hmc::get_clmem_phi()
{
	return clmem_phi;
}

#ifdef _PROFILING_
usetimer* Opencl_Module_Hmc::get_timer(char * in)
{
	usetimer *noop = NULL;
	noop = Opencl_Module_Fermions::get_timer(in);
	if(noop != NULL) return noop;

	if (strcmp(in, "generate_gaussian_spinorfield") == 0) {
		return &this->timer_generate_gaussian_spinorfield;
	}
	if (strcmp(in, "generate_gaussian_gaugemomenta") == 0) {
		return &this->timer_generate_gaussian_gaugemomenta;
	}
	if (strcmp(in, "md_update_gaugefield") == 0) {
		return &this->timer_md_update_gaugefield;
	}
	if (strcmp(in, "md_update_gaugemomenta") == 0) {
		return &this->timer_md_update_gaugemomenta;
	}
	if (strcmp(in, "gauge_force") == 0) {
		return &this->timer_gauge_force;
	}
	if (strcmp(in, "fermion_force") == 0) {
		return &this->timer_fermion_force;
	}
	if (strcmp(in, "set_zero_gaugemomentum") == 0) {
		return &this->timer_set_zero_gaugemomentum;
	}
	if (strcmp(in, "gaugemomentum_squarenorm") == 0) {
		return &this->timer_gaugemomentum_squarenorm;
	}
	if (strcmp(in, "stout_smear_fermion_force") == 0) {
		return &this->timer_stout_smear_fermion_force;
	}
	//if the kernelname has not matched, return NULL
	else {
		return NULL;
	}
}

int Opencl_Module_Hmc::get_read_write_size(char * in, inputparameters * parameters)
{
	Opencl_Module_Fermions::get_read_write_size(in, parameters);
	//Depending on the compile-options, one has different sizes...
	int D = (*parameters).get_float_size();
	int R = (*parameters).get_mat_size();
	int S;
	if((*parameters).get_use_eo() == 1)
		S = get_parameters()->get_eoprec_spinorfieldsize();
	else
		S = get_parameters()->get_spinorfieldsize();
	//this is the same as in the function above
	if (strcmp(in, "generate_gaussian_spinorfield") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "generate_gaussian_gaugemomenta") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "md_update_gaugefield") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "md_update_gaugemomenta") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "gauge_force") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "fermion_force") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "set_zero_gaugemomentum;") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "gaugemomentum_squarenorm") == 0) {
		return 10000000000000000000;
	}
	if (strcmp(in, "stout_smear_fermion_force") == 0) {
		return 10000000000000000000;
	}
	return 0;
}

void Opencl_Module_Hmc::print_profiling(std::string filename)
{
	Opencl_Module_Fermions::print_profiling(filename);
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

////////////////////////////////////////////////////
//Methods needed for the HMC-algorithm

void Opencl_Module_Hmc::generate_gaussian_gaugemomenta_device()
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(generate_gaussian_gaugemomenta, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_new_p
	int clerr = clSetKernelArg(generate_gaussian_gaugemomenta, 0, sizeof(cl_mem), &clmem_p);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(generate_gaussian_gaugemomenta, 1, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( generate_gaussian_gaugemomenta , gs2, ls2);
}

void Opencl_Module_Hmc::generate_gaussian_spinorfield_device()
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(generate_gaussian_spinorfield, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_phi_inv, which can be done since the gaussian field is only needed in the beginning
	int clerr = clSetKernelArg(generate_gaussian_spinorfield, 0, sizeof(cl_mem), &clmem_phi_inv);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(generate_gaussian_spinorfield, 1, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(generate_gaussian_spinorfield  , gs2, ls2);
}

void Opencl_Module_Hmc::md_update_spinorfield()
{
	//suppose the initial gaussian field is saved in clmem_phi_inv (see above).
	//  then the "phi" = Dpsi from the algorithm is stored in clmem_phi
	//  which then has to be the source of the inversion
	Opencl_Module_Fermions::Qplus(clmem_phi_inv, clmem_phi , *get_gaugefield());

	//debugging
	hmc_float s_fermion;
	set_float_to_global_squarenorm_device(clmem_phi, clmem_s_fermion);
	get_buffer_from_device(clmem_s_fermion, &s_fermion, sizeof(hmc_float));
	logger.debug() << "\tsquarenorm init field after update = " << s_fermion;
}


void Opencl_Module_Hmc::calc_fermion_force(usetimer * solvertimer)
{
	if(get_parameters()->get_use_eo() == true){
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		if(get_parameters()->get_use_cg() == true) { 
		//this is broken right now since the CG doesnt work!!
			throw Print_Error_Message("\t\tcalc fermion force ingredients using cg is not implemented yet. Aborting..");
		} else  {
			logger.debug() << "\t\tcalc fermion force ingredients using bicgstab";
	}
	else{
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield 
		if(get_parameters()->get_use_cg() == true) { 
			//this is broken right now since the CG doesnt work!!
			throw Print_Error_Message("\t\tcalc fermion force ingredients using cg is not implemented yet. Aborting..");
		} else  {
			logger.debug() << "\t\tcalc fermion force ingredients using bicgstab";

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
			set_spinorfield_cold_device(get_clmem_inout());

			//debugging
			hmc_float s_fermion;
			set_float_to_global_squarenorm_device(get_clmem_inout(), clmem_s_fermion);
			get_buffer_from_device(clmem_s_fermion, &s_fermion, sizeof(hmc_float));
			logger.debug() << "\tsquarenorm of inv.field before = " << s_fermion;

			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			Opencl_Module_Fermions::solver(Qplus_call, this->get_clmem_inout(), this->get_clmem_phi(), this->clmem_new_u, solvertimer);

			//debugging
			set_float_to_global_squarenorm_device(get_clmem_inout(), clmem_s_fermion);
			get_buffer_from_device(clmem_s_fermion, &s_fermion, sizeof(hmc_float));
			logger.debug() << "\tsquarenorm of inv.field after = " << s_fermion;

			//store this result in clmem_phi_inv
			copy_buffer_on_device(get_clmem_inout(), clmem_phi_inv, sizeof(spinor) * get_parameters()->get_spinorfieldsize());

			/**
			 * Now, one has to calculate
			 * X = (Qminus)^-1 Y = (Qminus)^-1 (Qplus)^-1 psi = (QplusQminus)^-1 psi ??
			 * out of
			 * Qminus clmem_inout = clmem_phi_inv
			 * Therefore, when calculating the final energy of the spinorfield,
			 *  one can just take clmem_phi_inv (see also above)!!
			 */

			//copy former solution to clmem_source
			copy_buffer_on_device(get_clmem_inout(), get_clmem_source(), sizeof(spinor) * get_parameters()->get_spinorfieldsize());
			logger.debug() << "\t\t\tstart solver";

			//debugging
			set_float_to_global_squarenorm_device(get_clmem_inout(), clmem_s_fermion);
			get_buffer_from_device(clmem_s_fermion, &s_fermion, sizeof(hmc_float));
			logger.debug() << "\tsquarenorm of inv.field before = " << s_fermion;

			//this sets clmem_inout cold as trial-solution
			set_spinorfield_cold_device(get_clmem_inout());

			Opencl_Module_Fermions::solver(Qminus_call, get_clmem_inout(), get_clmem_source(), clmem_new_u, solvertimer);

			//debugging
			set_float_to_global_squarenorm_device(get_clmem_inout(), clmem_s_fermion);
			get_buffer_from_device(clmem_s_fermion, &s_fermion, sizeof(hmc_float));
			logger.debug() << "\tsquarenorm of inv.field after = " << s_fermion;
		}
		logger.debug() << "\t\tcalc fermion_force...";
		//CP: this always calls fermion_force(Y,X) with Y = clmem_phi_inv, X = clmem_inout
		fermion_force_device();
	}
}

void Opencl_Module_Hmc::calc_gauge_force()
{
	logger.debug() << "\t\tcalc gauge_force...";
	gauge_force_device();
}

hmc_observables Opencl_Module_Hmc::metropolis(hmc_float rnd, hmc_float beta)
{
	//Calc Hamiltonian
	logger.debug() << "Calculate Hamiltonian";

	//Gauge-Part
	hmc_float tplaq, splaq, plaq;
	hmc_float tplaq_new, splaq_new, plaq_new;
	hmc_complex poly;
	hmc_complex poly_new;
	//In this call, the observables are calculated already with appropiate Weighting factor of 2.0/(VOL4D*NDIM*(NDIM-1)*NC)
	Opencl_Module::gaugeobservables(*get_gaugefield(), &plaq,  &tplaq, &splaq, &poly);
	Opencl_Module::gaugeobservables(clmem_new_u, &plaq_new,  &tplaq_new, &splaq_new, &poly_new);
	//plaq has to be divided by the norm-factor and multiplied by NC
	//  (because this is in the defintion of the gauge action and not in the normalization) to get s_gauge
	hmc_float factor = 2.0 / static_cast<hmc_float>(parameters->get_vol4d() * NDIM * (NDIM - 1) );
	/** NOTE: the minus here is introduced to fit tmlqcd!!! */
	hmc_float deltaH = -(plaq - plaq_new) * beta / factor;

	logger.debug() << "\tS_gauge(old field) = " << setprecision(10) << plaq << "\t" << plaq* beta  / factor;
	logger.debug() << "\tS_gauge(new field) = " << setprecision(10) << plaq_new << "\t" << plaq_new* beta / factor;
	logger.debug() << "\tdeltaS_gauge = " << setprecision(10) << deltaH;

	//Gaugemomentum-Part
	hmc_float p2, new_p2;
	set_float_to_gaugemomentum_squarenorm_device(clmem_p, clmem_p2);
	set_float_to_gaugemomentum_squarenorm_device(clmem_new_p, clmem_new_p2);
	Opencl_Module_Hmc::get_buffer_from_device(clmem_p2, &p2, sizeof(hmc_float));
	Opencl_Module_Hmc::get_buffer_from_device(clmem_new_p2, &new_p2, sizeof(hmc_float));
	//the energy is half the squarenorm
	deltaH += 0.5 * (p2 - new_p2);

	logger.debug() << "\tS_gaugemom(old field) = " << setprecision(10) << 0.5 * p2;
	logger.debug() << "\tS_gaugemom(new field) = " << setprecision(10) << 0.5 * new_p2;
	logger.debug() << "\tdeltaS_gaugemom = " << setprecision(10) << 0.5 * (p2 - new_p2);

	//Fermion-Part:
	hmc_float spinor_energy_init, s_fermion;
	//initial energy has been computed in the beginning...
	Opencl_Module_Hmc::get_buffer_from_device(clmem_energy_init, &spinor_energy_init, sizeof(hmc_float));
	// sum_links phi*_i (M^+M)_ij^-1 phi_j
	// the inversion with respect to the input-gaussian field and the new gaugefield has taken place during the leapfrog
	//  and was saved in clmem_phi_inv during that step.
	set_float_to_global_squarenorm_device(clmem_phi_inv, clmem_s_fermion);
	get_buffer_from_device(clmem_s_fermion, &s_fermion, sizeof(hmc_float));
	deltaH += spinor_energy_init - s_fermion;

	logger.debug() << "\tS_ferm(old field) = " << setprecision(10) <<  spinor_energy_init;
	logger.debug() << "\tS_ferm(new field) = " << setprecision(10) << s_fermion;
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

void Opencl_Module_Hmc::calc_spinorfield_init_energy()
{
	//Suppose the initial spinorfield is saved in phi_inv
	//  it is created in generate_gaussian_spinorfield_device
	Opencl_Module_Fermions::set_float_to_global_squarenorm_device(clmem_phi_inv, clmem_energy_init);
}

void Opencl_Module_Hmc::md_update_gaugemomentum_device(hmc_float eps)
{
	//__kernel void md_update_gaugemomenta(hmc_float eps, __global ae * p_inout, __global ae* force_in){
	hmc_float tmp = eps;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(md_update_gaugemomenta, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(md_update_gaugemomenta, 0, sizeof(hmc_float), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(md_update_gaugemomenta, 1, sizeof(cl_mem), &clmem_new_p);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(md_update_gaugemomenta, 2, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(md_update_gaugemomenta  , gs2, ls2);
}

void Opencl_Module_Hmc::md_update_gaugefield_device(hmc_float eps)
{
	// __kernel void md_update_gaugefield(hmc_float eps, __global ae * p_in, __global ocl_s_gaugefield * u_inout){
	hmc_float tmp = eps;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(md_update_gaugefield, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_force
	int clerr = clSetKernelArg(md_update_gaugefield, 0, sizeof(hmc_float), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(md_update_gaugefield, 1, sizeof(cl_mem), &clmem_new_p);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(md_update_gaugefield, 2, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( md_update_gaugefield , gs2, ls2);
}

void Opencl_Module_Hmc::set_zero_clmem_force_device()
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(md_update_gaugemomenta, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_force
	int clerr = clSetKernelArg(set_zero_gaugemomentum, 0, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( set_zero_gaugemomentum , gs2, ls2);
}

void Opencl_Module_Hmc::gauge_force_device()
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gauge_force, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gauge_force, 0, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(gauge_force, 1, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( gauge_force , gs2, ls2);
}

void Opencl_Module_Hmc::fermion_force_device()
{
	//fermion_force(field, Y, X, out);
	cl_mem tmp = get_clmem_inout();
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(fermion_force, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(fermion_force, 0, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 1, sizeof(cl_mem), &clmem_phi_inv);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 2, sizeof(cl_mem), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 3, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( fermion_force , gs2, ls2);
}

void Opencl_Module_Hmc::stout_smeared_fermion_force_device()
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(stout_smear_fermion_force, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
}

void Opencl_Module_Hmc::set_float_to_gaugemomentum_squarenorm_device(cl_mem clmem_in, cl_mem clmem_out)
{
	//__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gaugemomentum_squarenorm, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out){
	int clerr = clSetKernelArg(gaugemomentum_squarenorm, 0, sizeof(cl_mem), &clmem_in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

//  /** @todo add reduction */
	clerr = clSetKernelArg(gaugemomentum_squarenorm, 1, sizeof(cl_mem), &clmem_out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(gaugemomentum_squarenorm  , gs2, ls2);
}
