#include "opencl_module_hmc.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"
#include "meta/util.hpp"

using namespace std;

void Opencl_Module_Hmc::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module_Fermions::fill_collect_options(collect_options);
	*collect_options <<  " -DBETA=" << get_parameters().get_beta() << " -DGAUGEMOMENTASIZE=" << meta::get_vol4d(get_parameters()) * NDIM;
	//in case of tlsym gauge action
	if(meta::get_use_rectangles(get_parameters()) == true) {
		*collect_options <<  " -DC0=" << meta::get_c0(get_parameters()) << " -DC1=" << meta::get_c1(get_parameters());
	}
	if(use_soa) {
		*collect_options << " -DGAUGEMOMENTA_STRIDE=" << calculateStride(meta::get_vol4d(get_parameters()) * NDIM, sizeof(hmc_float));
	}
	return;
}


void Opencl_Module_Hmc::fill_buffers()
{

	Opencl_Module_Fermions::fill_buffers();
	///@todo CP: some of the above buffers are not used and can be deleted again!! especially in the eo-case

	int spinorfield_size = meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor);
	int eoprec_spinorfield_size = get_eoprec_spinorfield_buffer_size();
	int gaugemomentum_size = get_gaugemomentum_buffer_size();
	int gaugefield_size = getGaugefieldBufferSize();
	int float_size = sizeof(hmc_float);
	hmc_complex one = hmc_complex_one;
	hmc_complex minusone = hmc_complex_minusone;

	//init mem-objects

	logger.trace() << "Create buffer for HMC...";
	clmem_force = create_rw_buffer(gaugemomentum_size);
	if(get_parameters().get_use_eo() == true) {
		///@TODO in this case, the objects cl_mem_inout, source, tmp from the fermions module can be released again!!
		clmem_phi_inv_eo = create_rw_buffer(eoprec_spinorfield_size);
		clmem_phi_eo = create_rw_buffer(eoprec_spinorfield_size);
		//in debug-mode, this field is currently used temporarily...
		if(logger.beDebug()) clmem_phi_inv = create_rw_buffer(spinorfield_size);
		if(get_parameters().get_use_mp() ) {
			clmem_phi_mp_eo = create_rw_buffer(eoprec_spinorfield_size);
		}
	} else {
		///@TODO in this case, the object cl_mem_source from the fermions module can be released again!!
		clmem_phi = create_rw_buffer(spinorfield_size);
		clmem_phi_inv = create_rw_buffer(spinorfield_size);
		if(get_parameters().get_use_mp() ) {
			clmem_phi_mp = create_rw_buffer(spinorfield_size);
		}
	}

	clmem_new_u = create_rw_buffer(gaugefield_size);
	clmem_p = create_rw_buffer(gaugemomentum_size);
	clmem_new_p = create_rw_buffer(gaugemomentum_size);
	clmem_s_fermion_init = create_rw_buffer(float_size);
	if(get_parameters().get_use_mp() ) {
		clmem_s_fermion_mp_init = create_rw_buffer(float_size);
	}
	clmem_p2 = create_rw_buffer(float_size);
	clmem_new_p2 = create_rw_buffer(float_size);
	clmem_s_fermion = create_rw_buffer(float_size);

	return;
}

void Opencl_Module_Hmc::fill_kernels()
{
	Opencl_Module_Fermions::fill_kernels();

	basic_hmc_code = basic_fermion_code << "types_hmc.h" << "operations_gaugemomentum.cl";

	//init kernels for HMC
	if(get_parameters().get_use_eo() == true) {
		generate_gaussian_spinorfield_eo = createKernel("generate_gaussian_spinorfield_eo") << basic_hmc_code << prng_code << "spinorfield_eo_gaussian.cl";
		fermion_force_eo = createKernel("fermion_force_eo") << basic_hmc_code << "fermionmatrix.cl" << "force_fermion_eo.cl";
	} else {
		generate_gaussian_spinorfield = createKernel("generate_gaussian_spinorfield") << basic_hmc_code << prng_code << "spinorfield_gaussian.cl";
	}
	fermion_force = createKernel("fermion_force") << basic_hmc_code << "fermionmatrix.cl" << "force_fermion.cl";
	set_zero_gaugemomentum = createKernel("set_zero_gaugemomentum") << basic_hmc_code <<  "gaugemomentum_zero.cl";
	generate_gaussian_gaugemomenta = createKernel("generate_gaussian_gaugemomenta") << basic_hmc_code << prng_code << "gaugemomentum_gaussian.cl";
	md_update_gaugefield = createKernel("md_update_gaugefield") << basic_hmc_code << "md_update_gaugefield.cl";
	md_update_gaugemomenta = createKernel("md_update_gaugemomenta") << basic_hmc_code  << "md_update_gaugemomenta.cl";
	gauge_force = createKernel("gauge_force") << basic_hmc_code  << "force_gauge.cl";
	if(meta::get_use_rectangles(get_parameters()) == true) {
		//at the time of writing this kernel, the OpenCL compiler crashed the kernel using optimizations
		gauge_force_tlsym = createKernel("gauge_force_tlsym") << basic_hmc_code << "force_gauge_tlsym.cl";
	}
	if(get_parameters().get_use_smearing() == true) {
		stout_smear_fermion_force = createKernel("stout_smear_fermion_force") << basic_hmc_code << "force_fermion_stout_smear.cl";
	}
	gaugemomentum_squarenorm = createKernel("gaugemomentum_squarenorm") << basic_hmc_code << "gaugemomentum_squarenorm.cl";
	if(use_soa) {
		gaugemomentum_convert_to_soa = createKernel("gaugemomentum_convert_to_soa") << basic_hmc_code << "gaugemomentum_convert.cl";
		gaugemomentum_convert_from_soa = createKernel("gaugemomentum_convert_from_soa") << basic_hmc_code << "gaugemomentum_convert.cl";
	} else {
		gaugemomentum_convert_to_soa = 0;
		gaugemomentum_convert_from_soa = 0;
	}
}

void Opencl_Module_Hmc::clear_kernels()
{
	Opencl_Module_Fermions::clear_kernels();

	cl_uint clerr = CL_SUCCESS;

	logger.debug() << "release HMC-kernels.." ;
	if(get_parameters().get_use_eo() == true) {
		clerr = clReleaseKernel(generate_gaussian_spinorfield_eo);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(fermion_force_eo);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	} else {
		clerr = clReleaseKernel(generate_gaussian_spinorfield);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(fermion_force);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	clerr = clReleaseKernel(generate_gaussian_gaugemomenta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(md_update_gaugefield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(md_update_gaugemomenta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(gauge_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(meta::get_use_rectangles(get_parameters()) == true) {
		clerr = clReleaseKernel(gauge_force_tlsym);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	clerr = clReleaseKernel(set_zero_gaugemomentum);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(get_parameters().get_use_smearing() == true) {
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
	clerr = clReleaseMemObject(clmem_s_fermion_init);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	if(get_parameters().get_use_mp() ) {
		clerr = clReleaseMemObject(clmem_s_fermion_mp_init);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}
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
	clerr = clReleaseMemObject(clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	if(get_parameters().get_use_eo() == true) {
		clerr = clReleaseMemObject(clmem_phi_inv_eo);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(clmem_phi_eo);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		if(get_parameters().get_use_mp() ) {
			clerr = clReleaseMemObject(clmem_phi_mp_eo);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		}
	} else {
		clerr = clReleaseMemObject(clmem_phi_inv);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(clmem_phi);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		if(get_parameters().get_use_mp() ) {
			clerr = clReleaseMemObject(clmem_phi_mp);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		}
	}
	return;
}

void Opencl_Module_Hmc::get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups)
{
	Opencl_Module_Fermions::get_work_sizes(kernel, dev_type, ls, gs, num_groups);

	// kernels that use random numbers must not exceed the size of the random state array
	if(kernel == generate_gaussian_gaugemomenta
	   || kernel == generate_gaussian_spinorfield
	   || kernel == generate_gaussian_spinorfield_eo) {
		if(*gs > get_num_rndstates()) {
			*gs = get_num_rndstates();
		}
	}

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

cl_mem Opencl_Module_Hmc::get_clmem_phi_eo()
{
	return clmem_phi_eo;
}

cl_mem Opencl_Module_Hmc::get_clmem_phi_mp()
{
	return clmem_phi_mp;
}

cl_mem Opencl_Module_Hmc::get_clmem_phi_mp_eo()
{
	return clmem_phi_mp_eo;
}

cl_mem Opencl_Module_Hmc::get_clmem_s_fermion_init()
{
	return clmem_s_fermion_init;
}

cl_mem Opencl_Module_Hmc::get_clmem_s_fermion_mp_init()
{
	return clmem_s_fermion_mp_init;
}

#ifdef _PROFILING_
usetimer* Opencl_Module_Hmc::get_timer(const char * in)
{
	logger.trace() << "Opencl_Module_Hmc::get_timer(char*)";
	usetimer *noop = NULL;
	noop = Opencl_Module_Fermions::get_timer(in);
	if(noop != NULL) return noop;

	if (strcmp(in, "generate_gaussian_spinorfield") == 0) {
		return &this->timer_generate_gaussian_spinorfield;
	}
	if (strcmp(in, "generate_gaussian_spinorfield_eo") == 0) {
		return &this->timer_generate_gaussian_spinorfield_eo;
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
	if (strcmp(in, "gauge_force_tlsym") == 0) {
		return &this->timer_gauge_force_tlsym;
	}
	if (strcmp(in, "fermion_force") == 0) {
		return &this->timer_fermion_force;
	}
	if (strcmp(in, "fermion_force_eo") == 0) {
		return &this->timer_fermion_force_eo;
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

#endif

size_t Opencl_Module_Hmc::get_read_write_size(char * in)
{
	size_t result = Opencl_Module_Fermions::get_read_write_size(in);
	if (result != 0) return result;
//Depending on the compile-options, one has different sizes...
	size_t D = meta::get_float_size(parameters);
	//this returns the number of entries in an su3-matrix
	size_t R = meta::get_mat_size(parameters);
	//this is the number of spinors in the system (or number of sites)
	size_t S = meta::get_spinorfieldsize(get_parameters());
	size_t Seo = meta::get_eoprec_spinorfieldsize(get_parameters());
	//this is the number of links in the system (and of gaugemomenta)
	size_t G = meta::get_vol4d(get_parameters()) * NDIM;
	//factor for complex numbers
	int C = 2;
	//this is the same as in the function above
	//NOTE: 1 spinor has NC*NDIM = 12 complex entries
	//NOTE: 1 ae has NC*NC-1 = 8 real entries
	int A = meta::get_su3algebrasize();
	if (strcmp(in, "generate_gaussian_spinorfield") == 0) {
		//this kernel writes 1 spinor
		return ( 12 * C ) * D * S;
	}
	if (strcmp(in, "generate_gaussian_spinorfield_eo") == 0) {
		//this kernel writes 1 spinor
		return ( 12 * C ) * D * Seo;
	}
	if (strcmp(in, "generate_gaussian_gaugemomenta") == 0) {
		//this kernel writes 1 ae
		return (A) * D * G;
	}
	if (strcmp(in, "md_update_gaugefield") == 0) {
		//this kernel reads 1 ae and 1 su3 matrix and writes 1 su3 matrix for every link
		return (A + (1 + 1) * R * C) * D * G;
	}
	if (strcmp(in, "md_update_gaugemomenta") == 0) {
		//this kernel reads 2 ae and writes 1 ae per link
		return ((2 + 1) * A) * D * G;
	}
	if (strcmp(in, "gauge_force") == 0) {
		//this kernel reads ingredients for 1 staple plus 1 su3matrix and writes 1 ae for every link
		return G * D * (R * C * ( 6 * (NDIM - 1) + 1 ) + A );
	}
	if (strcmp(in, "gauge_force_tlsym") == 0) {
		//this kernel reads ingredients for 1 rect-staple plus 1 su3matrix and writes 1 ae for every link
		//the rect staple is the same as the normal staple, but with 3 add. contributions and 2 add. matrices in each contribution, so instead of 2 * 3 = 6, one reads 6 * 5 matrices in each direction
		return G * D * (R * C * ( 6 * 5 * (NDIM - 1 ) + 1 ) + A );
	}
	if (strcmp(in, "fermion_force") == 0) {
		//this kernel reads 16 spinors, 8 su3matrices and writes 8 ae per site
		//NOTE: the kernel now runs over all ae instead of all sites, but this must be equivalent!
		return (C * 12 * (16) + C * 8 * R + 8 * A) * D * S;
	}
	if (strcmp(in, "fermion_force_eo") == 0) {
		//this kernel reads 16 spinors, 8 su3matrices and writes 1 ae per site
		return (C * 12 * (16) + C * 8 * R + 8 * A) * D * Seo;
	}
	if (strcmp(in, "set_zero_gaugemomentum;") == 0) {
		//this kernel writes 1 ae per link
		return G * D * A;
	}
	if (strcmp(in, "gaugemomentum_squarenorm") == 0) {
		//this kernel reads 1 ae and writes 1 float per link
		return G * D * ( A + 1);
	}
	if (strcmp(in, "stout_smear_fermion_force") == 0) {
		return 10000000000000000000;
	}
	return 0;
}

uint64_t Opencl_Module_Hmc::get_flop_size(const char * in)
{
	uint64_t result = Opencl_Module_Fermions::get_flop_size(in);
	if (result != 0) return result;
	//this is the number of spinors in the system (or number of sites)
	size_t S = meta::get_spinorfieldsize(get_parameters());
	size_t Seo = meta::get_eoprec_spinorfieldsize(get_parameters());
	//this is the number of links in the system (and of gaugemomenta)
	uint64_t G = meta::get_vol4d(get_parameters()) * NDIM;
	//NOTE: 1 ae has NC*NC-1 = 8 real entries
	uint64_t A = meta::get_su3algebrasize();
	//this returns the number of entries in an su3-matrix
	uint64_t R = meta::get_mat_size(parameters);
	//this is the same as in the function above
	if (strcmp(in, "generate_gaussian_spinorfield") == 0) {
		//this kernel performs 12 multiplications per site
		///@todo ? I did not count the gaussian normal pair production, which is very complicated...
		return 12 * S;
	}
	if (strcmp(in, "generate_gaussian_spinorfield_eo") == 0) {
		//this kernel performs 12 multiplications per site
		///@todo ? I did not count the gaussian normal pair production, which is very complicated...
		return 12 * Seo;
	}
	if (strcmp(in, "generate_gaussian_gaugemomenta") == 0) {
		//this kernel performs 0 multiplications per site
		///@todo ? I did not count the gaussian normal pair production, which is very complicated...
		return 0;
	}
	if (strcmp(in, "md_update_gaugefield") == 0) {
		//this kernel performs one exp(i ae) ( = 327 flops + 1 su3 mult ) and 1 su3 mult per link
		return (meta::get_flop_su3_su3() * ( 1 + 1)  + 327 ) * G;
	}
	if (strcmp(in, "md_update_gaugemomenta") == 0) {
		//this kernel performs 1 real mult and 1 real add per ae
		return (1 + 1) * A * G;
	}
	if (strcmp(in, "gauge_force") == 0) {
		//this kernel calculates 1 staple (= 4*ND-1 su3_su3 + 2*ND-1 su3_add), 1 su3*su3, 1 tr_lambda_u (19 flops) plus 8 add and 8 mult per ae
		return ( 4 * (NDIM - 1) * meta::get_flop_su3_su3() + 2 * (NDIM - 1) * 18 + 1 * meta::get_flop_su3_su3() + 19  + A * ( 1 + 1 )
		       ) * G;
	}
	if (strcmp(in, "gauge_force_tlsym") == 0) {
		//this kernel calculates 1 rect-staple (= 24*ND-1 su3_su3 + 6*ND-1 su3_add), 1 su3*su3, 1 tr_lambda_u (19 flops) plus 8 add and 8 mult per ae
		//24 = 6 contr. per dir, 4 mat_mat per contr.
		return ( 24 * (NDIM - 1) * meta::get_flop_su3_su3() + 6 * (NDIM - 1) * 18 + 1 * meta::get_flop_su3_su3() + 19  + A * ( 1 + 1 )
		       ) * G;
	}
	if (strcmp(in, "fermion_force") == 0) {
		//this kernel performs NDIM * ( 4 * su3vec_acc (6 flops) + tr(v*u) (126 flops) + tr_lambda_u(19 flops) + update_ae(8*2 flops) + su3*su3 + su3*complex (flop_complex_mult * R ) ) per site
		//NOTE: the kernel now runs over all ae instead of all sites, but this must be equivalent!
		return Seo * NDIM * ( 4 * 6 + 126 + 19 + 8 * 2 + meta::get_flop_su3_su3() + meta::get_flop_complex_mult() * R );
	}
	if (strcmp(in, "fermion_force_eo") == 0) {
		//this kernel performs NDIM * ( 4 * su3vec_acc (6 flops) + tr(v*u) (126 flops) + tr_lambda_u(19 flops) + update_ae(8*2 flops) + su3*su3 + su3*complex (flop_complex_mult * R ) ) per site
		return Seo * NDIM * ( 4 * 6 + 126 + 19 + 8 * 2 + meta::get_flop_su3_su3() + meta::get_flop_complex_mult() * R );
	}
	if (strcmp(in, "set_zero_gaugemomentum;") == 0) {
		//this kernel performs 0 mults
		return 0;
	}
	if (strcmp(in, "gaugemomentum_squarenorm") == 0) {
		//this kernel performs 8 real mults and 8-1 real adds per ae
		return (8 + 7) * A * G;
	}
	if (strcmp(in, "stout_smear_fermion_force") == 0) {
		return 10000000000000000000;
	}
	return 0;
}

#ifdef _PROFILING_

void Opencl_Module_Hmc::print_profiling(std::string filename, int number)
{
	Opencl_Module_Fermions::print_profiling(filename, number);
	char * kernelName;
	kernelName = "generate_gaussian_spinorfield";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "generate_gaussian_spinorfield_eo";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "generate_gaussian_gaugemomenta";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "md_update_gaugefield";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "md_update_gaugemomenta";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "gauge_force";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "gauge_force_tlsym";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "fermion_force";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "fermion_force_eo";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "set_zero_gaugemomentum";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "gaugemomentum_squarenorm";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "stout_smear_fermion_force";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
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
	int clerr = clSetKernelArg(generate_gaussian_gaugemomenta, 0, sizeof(cl_mem), &clmem_p);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(generate_gaussian_gaugemomenta, 1, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( generate_gaussian_gaugemomenta , gs2, ls2);

	if(logger.beDebug()) {
		cl_mem force_tmp = create_rw_buffer(sizeof(hmc_float));
		hmc_float resid;
		this->set_float_to_gaugemomentum_squarenorm_device(clmem_p, force_tmp);
		get_buffer_from_device(force_tmp, &resid, sizeof(hmc_float));
		logger.debug() <<  "\tgaussian gaugemomenta:\t" << resid;
		int clerr = clReleaseMemObject(force_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		if(resid != resid) {
			throw Print_Error_Message("calculation of gaussian gm gave nan! Aborting...", __FILE__, __LINE__);
		}
		if(resid == INFINITY) {
			bool writeout = false;
			if(writeout) {
				//create buffer to store ae-field
				int ae_num = meta::get_vol4d(parameters) * NDIM;

				ae * ae_tmp = new ae[ae_num];

				//get buffer from device
				cout << "copy buffer to host" << endl;
				exportGaugemomentumBuffer(ae_tmp, clmem_p);

				//write out to file
				ofstream out("clmem_p_at_inf");
				if(!out) {
					cout << "Cannot open file.\n";
				}
				for(int i = 0; i < ae_num; i++) {
					out << i << "\t" << ae_tmp[i].e0 << endl;
					out << i << "\t" << ae_tmp[i].e1 << endl;
					out << i << "\t" << ae_tmp[i].e2 << endl;
					out << i << "\t" << ae_tmp[i].e3 << endl;
					out << i << "\t" << ae_tmp[i].e4 << endl;
					out << i << "\t" << ae_tmp[i].e5 << endl;
					out << i << "\t" << ae_tmp[i].e6 << endl;
					out << i << "\t" << ae_tmp[i].e7 << endl;
				}
				out.close();

				//calc sqnorm of ae_tmp
				hmc_float sqnorm = 0.;
				for(int i = 0; i < ae_num; i++) {
					sqnorm += ae_tmp[i].e0 * ae_tmp[i].e0;
					sqnorm += ae_tmp[i].e1 * ae_tmp[i].e1;
					sqnorm += ae_tmp[i].e2 * ae_tmp[i].e2;
					sqnorm += ae_tmp[i].e3 * ae_tmp[i].e3;
					sqnorm += ae_tmp[i].e4 * ae_tmp[i].e4;
					sqnorm += ae_tmp[i].e5 * ae_tmp[i].e5;
					sqnorm += ae_tmp[i].e6 * ae_tmp[i].e6;
					sqnorm += ae_tmp[i].e7 * ae_tmp[i].e7;
				}
				cout << "sqnrom: " << sqnorm << endl;
				free(ae_tmp);
			}
			throw Print_Error_Message("calculation of gaussian gm gave inf! Aborting...", __FILE__, __LINE__);
		}

	}

}

void Opencl_Module_Hmc::generate_spinorfield_gaussian()
{
	if(get_parameters().get_use_eo() == true) {
		this->generate_gaussian_spinorfield_eo_device();
	} else {
		this->generate_gaussian_spinorfield_device();
	}
	return;
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

	if(logger.beDebug()) {
		cl_mem force_tmp = create_rw_buffer(sizeof(hmc_float));
		hmc_float resid;
		this->set_float_to_global_squarenorm_device(clmem_phi_inv, force_tmp);
		get_buffer_from_device(force_tmp, &resid, sizeof(hmc_float));
		logger.debug() <<  "\tinit gaussian spinorfield:\t" << resid;
		int clerr = clReleaseMemObject(force_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		if(resid != resid) {
			throw Print_Error_Message("calculation of gaussian spinorfield gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}

void Opencl_Module_Hmc::generate_gaussian_spinorfield_eo_device()
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(generate_gaussian_spinorfield_eo, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_phi_inv_eoporec, which can be done since the gaussian field is only needed in the beginning
	int clerr = clSetKernelArg(generate_gaussian_spinorfield_eo, 0, sizeof(cl_mem), &clmem_phi_inv_eo);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(generate_gaussian_spinorfield_eo, 1, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(generate_gaussian_spinorfield_eo, gs2, ls2);

	if(logger.beDebug()) {
		cl_mem force_tmp = create_rw_buffer(sizeof(hmc_float));
		hmc_float resid;
		this->set_float_to_global_squarenorm_eoprec_device(clmem_phi_inv_eo, force_tmp);
		get_buffer_from_device(force_tmp, &resid, sizeof(hmc_float));
		logger.debug() <<  "\tinit gaussian spinorfield energy:\t" << resid;
		int clerr = clReleaseMemObject(force_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		if(resid != resid) {
			throw Print_Error_Message("calculation of gaussian spinorfield gave nan! Aborting...", __FILE__, __LINE__);
		}
		logger.debug() << "\t\tforce calculated from gaussian spinorfield:";
		this->set_zero_clmem_force_device();
		fermion_force_eo_device(clmem_phi_inv_eo,  clmem_phi_inv_eo, EVEN);
		fermion_force_eo_device(clmem_phi_inv_eo,  clmem_phi_inv_eo, ODD);
		this->set_zero_clmem_force_device();
	}

}

void Opencl_Module_Hmc::md_update_spinorfield(hmc_float kappa, hmc_float mubar)
{
	//suppose the initial gaussian field is saved in clmem_phi_inv (see above).
	//  then the "phi" = Dpsi from the algorithm is stored in clmem_phi
	//  which then has to be the source of the inversion
	if(get_parameters().get_use_eo() == true) {
		Opencl_Module_Fermions::Qplus_eo (clmem_phi_inv_eo, clmem_phi_eo , get_gaugefield(), kappa, mubar);
		if(logger.beDebug()) print_info_inv_field(clmem_phi_eo, true, "\tinit field after update ");
	} else {
		Opencl_Module_Fermions::Qplus(clmem_phi_inv, clmem_phi , get_gaugefield(), kappa, mubar);
		if(logger.beDebug()) print_info_inv_field(clmem_phi, false, "\tinit field after update ");
	}
}

void Opencl_Module_Hmc::md_update_spinorfield_mp(usetimer * solvertimer)
{
	///@todo solvertimer is not used here yet...
	//suppose the initial gaussian field is saved in clmem_phi_inv (see above).
	//  then the "phi" = Dpsi from the algorithm is stored in clmem_phi
	//  which then has to be the source of the inversion
	//in the mass preconditioning case, this is a bit more complicated and involves an inversion
	if(get_parameters().get_use_eo() == true) {
		//CP: Init tmp spinorfield
		int spinorfield_size = sizeof(spinor) * meta::get_eoprec_spinorfieldsize(get_parameters());
		cl_mem sf_eo_tmp;
		sf_eo_tmp = create_rw_buffer(spinorfield_size);

		//sf_eo_tmp = Qplus_eo(light_mass) phi_inv_eo
		Opencl_Module_Fermions::Qplus_eo (clmem_phi_inv_eo, sf_eo_tmp , get_gaugefield());

		//Now one needs ( Qplus_eo )^-1 (heavy_mass) using sf_eo_tmp as source to get phi_mp_eo
		//use always bicgstab here
		logger.debug() << "\t\t\tstart solver";

		/** @todo at the moment, we can only put in a cold spinorfield
		 * or a point-source spinorfield as trial-solution
		 */
		set_zero_spinorfield_eoprec_device(get_clmem_phi_mp_eo());

		gamma5_eo_device(get_clmem_phi_mp_eo());

		int converged = -1;
		if(logger.beDebug()) print_info_inv_field(get_clmem_phi_mp_eo(), true, "\tinv. field before inversion ");
		if(logger.beDebug()) print_info_inv_field(sf_eo_tmp, true, "\tsource before inversion ");
		converged = Opencl_Module_Fermions::bicgstab_eo(::Qplus_eo(this), this->get_clmem_phi_mp_eo(), sf_eo_tmp, this->clmem_new_u, get_parameters().get_solver_prec(), get_parameters().get_kappa_mp(), meta::get_mubar_mp(parameters));
		if (converged < 0) {
			if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
			else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
		} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
		if(logger.beDebug()) print_info_inv_field(get_clmem_phi_mp_eo(), true, "\tinv. field after inversion ");
		//add number of inverter iterations to counter
		*inversions0 += abs(converged);
		if(logger.beDebug()) print_info_inv_field(clmem_phi_mp_eo, true, "\tinit field after update ");

		int clerr = clReleaseMemObject(sf_eo_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);

	} else {
		//CP: Init tmp spinorfield
		int spinorfield_size = sizeof(spinor) * get_spinorfieldsize(get_parameters());
		cl_mem sf_tmp;
		sf_tmp = create_rw_buffer(spinorfield_size);

		//sf_tmp = Qplus(light_mass) phi_inv
		Opencl_Module_Fermions::Qplus (clmem_phi_inv, sf_tmp , get_gaugefield());

		//Now one needs ( Qplus )^-1 (heavy_mass) using sf_tmp as source to get phi_mp
		//use always bicgstab here
		logger.debug() << "\t\t\tstart solver";

		/** @todo at the moment, we can only put in a cold spinorfield
		 * or a point-source spinorfield as trial-solution
		 */
		set_zero_spinorfield_device(get_clmem_phi_mp());
		gamma5_device(get_clmem_phi_mp());

		int converged = -1;
		if(logger.beDebug()) print_info_inv_field(get_clmem_phi_mp(), false, "\tinv. field before inversion ");
		if(logger.beDebug()) print_info_inv_field(sf_tmp, false, "\tsource before inversion ");
		converged = Opencl_Module_Fermions::bicgstab(::Qplus(this), this->get_clmem_phi_mp(), sf_tmp, this->clmem_new_u, get_parameters().get_solver_prec(), get_parameters().get_kappa_mp(), get_mubar_mp(parameters));
		if (converged < 0) {
			if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
			else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
		} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
		if(logger.beDebug()) print_info_inv_field(get_clmem_phi_mp(), false, "\tinv. field after inversion ");
		//add number of inverter iterations to counter
		*inversions0 += abs(converged);
		if(logger.beDebug()) print_info_inv_field(clmem_phi_mp, false, "\tinit field after update ");

		int clerr = clReleaseMemObject(sf_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}
}

//this function takes to args kappa and mubar because one has to use it with different masses when mass-prec is used and when not
void Opencl_Module_Hmc::calc_fermion_force(usetimer * solvertimer, hmc_float kappa, hmc_float mubar)
{
	int converged = -1;
	if(get_parameters().get_use_eo() == true) {
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			/**
			 * The first inversion calculates
			 * X_even = phi = (Qplusminus_eo)^-1 psi
			 * out of
			 * Qplusminus_eo phi_even = psi
			 */
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			/**
			 * Trial solution for the spinorfield
			 */
			set_eoprec_spinorfield_cold_device(get_clmem_inout_eo());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\t\t\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::cg_eo(::QplusQminus_eo(this), this->get_clmem_inout_eo(), this->get_clmem_phi_eo(), this->clmem_new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			//add number of inverter iterations to counter
			*inversions1 += abs(converged);
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\t\t\tinv. field after inversion ");

			/**
			 * Y_even is now just
			 *  Y_even = (Qminus_eo) X_even = (Qminus_eo) (Qplusminus_eo)^-1 psi =
			 *    = (Qplus_eo)^-1 psi
			 */
			Opencl_Module_Fermions::Qminus_eo(this->get_clmem_inout_eo(), clmem_phi_inv_eo, this->clmem_new_u, kappa, mubar);
		} else {
			///@todo if wanted, solvertimer has to be used here..
			//logger.debug() << "\t\tcalc fermion force ingredients using bicgstab with eo.";
			/**
			 * The first inversion calculates
			 * Y_even = phi = (Qplus_eo)^-1 psi
			 * out of
			 * Qplus_eo phi = psi
			 * This is also the energy of the final field!
			*/
			//logger.debug() << "\t\tcalc Y_even...";
			//logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			/**
			 * Trial solution for the spinorfield
			 */
			set_zero_spinorfield_eoprec_device(get_clmem_inout_eo());
			gamma5_eo_device(get_clmem_inout_eo());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\t\t\tinv. field before inversion ");
			if(logger.beDebug()) print_info_inv_field(get_clmem_phi_eo(), true, "\t\t\tsource before inversion ");
			converged = Opencl_Module_Fermions::bicgstab_eo(::Qplus_eo(this), this->get_clmem_inout_eo(), this->get_clmem_phi_eo(), this->clmem_new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\t\t\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions1 += abs(converged);

			//store this result in clmem_phi_inv
			copy_buffer_on_device(get_clmem_inout_eo(), clmem_phi_inv_eo, get_eoprec_spinorfield_buffer_size());

			/**
			 * Now, one has to calculate
			 * X = (Qminus_eo)^-1 Y = (Qminus_eo)^-1 (Qplus_eo)^-1 psi = (QplusQminus_eo)^-1 psi ??
			 * out of
			 * Qminus_eo clmem_inout_eo = clmem_phi_inv_eo
			 * Therefore, when calculating the final energy of the spinorfield,
			 *  one can just take clmem_phi_inv_eo (see also above)!!
			 */

			//logger.debug() << "\t\tcalc X_even...";
			//copy former solution to clmem_source
			copy_buffer_on_device(get_clmem_inout_eo(), get_clmem_source_even(), get_eoprec_spinorfield_buffer_size());

			//this sets clmem_inout cold as trial-solution
			set_eoprec_spinorfield_cold_device(get_clmem_inout_eo());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\t\t\tinv. field before inversion ");
			if(logger.beDebug()) print_info_inv_field(get_clmem_source_even(), true, "\t\t\tsource before inversion ");
			converged = Opencl_Module_Fermions::bicgstab_eo(::Qminus_eo(this), get_clmem_inout_eo(), get_clmem_source_even(), clmem_new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\t\t\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions1 += abs(converged);
		}
		/**
		 * At this point, one has calculated X_odd and Y_odd.
		 * If one has a fermionmatrix
		 *  M = R + D
		 * these are:
		 *  X_odd = -R(-mu)_inv D X_even
		 *  Y_odd = -R(mu)_inv D Y_even
		 */

		/** @fixme below usages of dslash should work, but only because we use bicgstab above
		    proper implementation needs to make sure this is always the case */

		///@NOTE the following calculations could also go in a new function for convenience
		//calculate X_odd
		//therefore, clmem_tmp_eo_1 is used as intermediate state. The result is saved in clmem_inout, since
		//  this is used as a default in the force-function.
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			dslash_eo_device(get_clmem_inout_eo(), get_clmem_tmp_eo_1(), clmem_new_u, ODD, kappa);
			sax_eoprec_device(get_clmem_tmp_eo_1(), get_clmem_minusone(), get_clmem_tmp_eo_1());
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			dslash_eo_device(get_clmem_inout_eo(), get_clmem_tmp_eo_1(), clmem_new_u, ODD, kappa);
			M_tm_inverse_sitediagonal_minus_device(get_clmem_tmp_eo_1(), get_clmem_tmp_eo_2(), mubar);
			sax_eoprec_device(get_clmem_tmp_eo_2(), get_clmem_minusone(), get_clmem_tmp_eo_1());
		}

		/**
		 * If needed, this gives additional debug informations
		 */
		bool debug_hard = false;
		cl_mem x_tmp, y_tmp;
		cl_mem x_eo_tmp, x_eo_tmp2, y_eo_tmp, y_eo_tmp2;

		if(logger.beDebug() && debug_hard) {
			int spinorfield_size = sizeof(spinor) * get_spinorfieldsize(get_parameters());
			x_tmp = create_rw_buffer(spinorfield_size);
			int eo_spinorfield_size = sizeof(spinor) * get_eoprec_spinorfieldsize(get_parameters());
			x_eo_tmp = create_rw_buffer(eo_spinorfield_size);
			x_eo_tmp2 = create_rw_buffer(eo_spinorfield_size);

			this->convert_from_eoprec_device(get_clmem_inout_eo(), get_clmem_tmp_eo_1(), x_tmp);
			print_info_inv_field(get_clmem_inout_eo(), true, "\t\t\t\tX_even ");
			print_info_inv_field(get_clmem_tmp_eo_1(), true, "\t\t\t\tX_odd ");
			print_info_inv_field(x_tmp, false, "\t\t\t\tX = (X_even, X_odd) ");

			//save x_even and x_odd temporarily
			copy_buffer_on_device(get_clmem_inout_eo(), x_eo_tmp, eo_spinorfield_size);
			copy_buffer_on_device(get_clmem_tmp_eo_1(), x_eo_tmp2, eo_spinorfield_size);
			print_info_inv_field(x_eo_tmp, true, "\t\t\tx_even:\t");
			print_info_inv_field(x_eo_tmp2, true, "\t\t\tx_odd:\t");
		}

		//logger.debug() << "\t\tcalc eo fermion_force F(Y_even, X_odd)...";
		//Calc F(Y_even, X_odd) = F(clmem_phi_inv_eo, clmem_tmp_eo_1)
		fermion_force_eo_device(clmem_phi_inv_eo,  get_clmem_tmp_eo_1(), EVEN, kappa);

		//calculate Y_odd
		//therefore, clmem_tmp_eo_1 is used as intermediate state. The result is saved in clmem_phi_inv, since
		//  this is used as a default in the force-function.
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			dslash_eo_device(clmem_phi_inv_eo, get_clmem_tmp_eo_1(), clmem_new_u, ODD, kappa);
			sax_eoprec_device(get_clmem_tmp_eo_1(), get_clmem_minusone(), get_clmem_tmp_eo_1());
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			dslash_eo_device(clmem_phi_inv_eo, get_clmem_tmp_eo_1(), clmem_new_u, ODD, kappa);
			M_tm_inverse_sitediagonal_device(get_clmem_tmp_eo_1(), get_clmem_tmp_eo_2(), mubar);
			sax_eoprec_device(get_clmem_tmp_eo_2(), get_clmem_minusone(), get_clmem_tmp_eo_1());
		}

		//logger.debug() << "\t\tcalc eoprec fermion_force F(Y_odd, X_even)...";
		//Calc F(Y_odd, X_even) = F(clmem_tmp_eo_1, clmem_inout_eo)
		fermion_force_eo_device(get_clmem_tmp_eo_1(), get_clmem_inout_eo(), ODD, kappa);

		if(logger.beDebug() && debug_hard) {
			int spinorfield_size = sizeof(spinor) * get_spinorfieldsize(get_parameters());
			y_tmp = create_rw_buffer(spinorfield_size);
			int eo_spinorfield_size = sizeof(spinor) * get_eoprec_spinorfieldsize(get_parameters());
			y_eo_tmp = create_rw_buffer(eo_spinorfield_size);
			y_eo_tmp2 = create_rw_buffer(eo_spinorfield_size);
			//tmp field for differences
			cl_mem sf_diff = create_rw_buffer(eo_spinorfield_size);

			this->convert_from_eoprec_device(clmem_phi_inv_eo, get_clmem_tmp_eo_1(), y_tmp);
			print_info_inv_field(clmem_phi_inv_eo, true, "\t\t\t\tY_even ");
			print_info_inv_field(get_clmem_tmp_eo_1(), true, "\t\t\t\tY_odd ");
			print_info_inv_field(y_tmp, false, "\t\t\t\tY = (Y_even, Yodd) ");

			//save y_even and y_odd temporarily
			copy_buffer_on_device(clmem_phi_inv_eo, y_eo_tmp, eo_spinorfield_size);
			copy_buffer_on_device(get_clmem_tmp_eo_1(), y_eo_tmp2, eo_spinorfield_size);
			print_info_inv_field(y_eo_tmp, true, "\t\t\ty_even:\t");
			print_info_inv_field(y_eo_tmp2, true, "\t\t\ty_odd:\t");

			hmc_float compare_to_non_eo_diff;
			hmc_complex compare_to_non_eo_scal;
			hmc_complex scalprod;
			cl_mem scal_tmp = create_rw_buffer(sizeof(hmc_complex));

			logger.debug() << "\t\t\tperform some tests on eo force ingredients:";
			//calculate scalar product
			this->set_complex_to_scalar_product_eoprec_device(x_eo_tmp, x_eo_tmp2, scal_tmp);
			get_buffer_from_device(scal_tmp, &scalprod, sizeof(hmc_complex));
			logger.debug() << std::scientific << "\t\t\t(x_even, x_odd):\t" << scalprod.re << "\t" << scalprod.im;
			//calculate difference
			this->saxpy_eoprec_device(x_eo_tmp2, x_eo_tmp, get_clmem_one(), sf_diff);
			print_info_inv_field(sf_diff, true, "\t\t\t(x_even - x_odd):\t");

			//calculate scalar product
			this->set_complex_to_scalar_product_eoprec_device(y_eo_tmp, y_eo_tmp2, scal_tmp);
			get_buffer_from_device(scal_tmp, &scalprod, sizeof(hmc_complex));
			logger.debug() << std::scientific << "\t\t\t(y_even, y_odd):\t" << scalprod.re << "\t" << scalprod.im;
			//calculate difference
			this->saxpy_eoprec_device(y_eo_tmp2, y_eo_tmp, get_clmem_one(), sf_diff);
			print_info_inv_field(sf_diff, true, "\t\t\t(y_even - y_odd):\t");

			//calculate scalar product
			this->set_complex_to_scalar_product_eoprec_device(x_eo_tmp, y_eo_tmp2, scal_tmp);
			get_buffer_from_device(scal_tmp, &scalprod, sizeof(hmc_complex));
			logger.debug() << std::scientific << "\t\t\t(x_even, y_odd):\t" << scalprod.re << "\t" << scalprod.im;
			//calculate difference
			this->saxpy_eoprec_device(y_eo_tmp2, x_eo_tmp, get_clmem_one(), sf_diff);
			print_info_inv_field(sf_diff, true, "\t\t\t(x_even - y_odd):\t");

			//calculate scalar product
			this->set_complex_to_scalar_product_eoprec_device(y_eo_tmp, x_eo_tmp2, scal_tmp);
			get_buffer_from_device(scal_tmp, &scalprod, sizeof(hmc_complex));
			logger.debug() << std::scientific << "\t\t\t(y_even, x_odd):\t" << scalprod.re << "\t" << scalprod.im;
			//calculate difference
			this->saxpy_eoprec_device(x_eo_tmp2, y_eo_tmp, get_clmem_one(), sf_diff);
			print_info_inv_field(sf_diff, true, "\t\t\t(y_even - x_odd):\t");

			//calculate scalar product
			this->set_complex_to_scalar_product_eoprec_device(x_eo_tmp, y_eo_tmp, scal_tmp);
			get_buffer_from_device(scal_tmp, &scalprod, sizeof(hmc_complex));
			logger.debug() << std::scientific << "\t\t\t(x_even, y_even):\t" << scalprod.re << "\t" << scalprod.im;

			compare_to_non_eo_scal.re = scalprod.re;
			compare_to_non_eo_scal.im = scalprod.im;

			//calculate difference
			this->saxpy_eoprec_device(x_eo_tmp, y_eo_tmp, get_clmem_one(), sf_diff);
			compare_to_non_eo_diff = print_info_inv_field(sf_diff, true, "\t\t\t(y_even - x_even):\t");

			//calculate scalar product
			this->set_complex_to_scalar_product_eoprec_device(x_eo_tmp2, y_eo_tmp2, scal_tmp);
			get_buffer_from_device(scal_tmp, &scalprod, sizeof(hmc_complex));
			logger.debug() << std::scientific << "\t\t\t(x_even, y_odd):\t" << scalprod.re << "\t" << scalprod.im;
			//calculate difference
			this->saxpy_eoprec_device(x_eo_tmp2, y_eo_tmp2, get_clmem_one(), sf_diff);
			compare_to_non_eo_diff += print_info_inv_field(sf_diff, true, "\t\t\t(x_odd - y_odd):\t");

			compare_to_non_eo_scal.re += scalprod.re;
			compare_to_non_eo_scal.im += scalprod.im;

			logger.debug() << "\t\t\tnon-eo result should be:";
			logger.debug() << std::scientific <<  "\t\t\t(x,y):\t" << compare_to_non_eo_scal.re << " " << compare_to_non_eo_scal.im;
			logger.debug() << std::scientific << "\t\t\t(x-y):\t" << compare_to_non_eo_diff ;

			int clerr = clReleaseMemObject(scal_tmp);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clerr = clReleaseMemObject(sf_diff);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);

		}

		if(logger.beDebug() && debug_hard) {
			logger.debug() << "\t\t\tperform checks on eo force:";
			logger.debug() << "\t\t\tF(x_e, x_o):";
			this->set_zero_clmem_force_device();
			fermion_force_eo_device(x_eo_tmp, x_eo_tmp2, EVEN, kappa);
			logger.debug() << "\t\t\tF(x_e, y_o):";
			this->set_zero_clmem_force_device();
			fermion_force_eo_device(x_eo_tmp, y_eo_tmp2, EVEN, kappa);
			logger.debug() << "\t\t\tF(x_e, y_e):";
			this->set_zero_clmem_force_device();
			fermion_force_eo_device(x_eo_tmp, y_eo_tmp, ODD, kappa);
			logger.debug() << "\t\t\tF(y_e, y_o):";
			this->set_zero_clmem_force_device();
			fermion_force_eo_device(y_eo_tmp, y_eo_tmp2, EVEN, kappa);
			logger.debug() << "\t\t\tF(y_e, x_o):";
			this->set_zero_clmem_force_device();
			fermion_force_eo_device(y_eo_tmp, x_eo_tmp2, EVEN, kappa);
			logger.debug() << "\t\t\tF(y_e, x_e):";
			this->set_zero_clmem_force_device();
			fermion_force_eo_device(y_eo_tmp, x_eo_tmp, ODD, kappa);
			logger.debug() << "\t\t\tF(x_o, x_e):";
			this->set_zero_clmem_force_device();
			fermion_force_eo_device(x_eo_tmp2, x_eo_tmp, ODD, kappa);
			logger.debug() << "\t\t\tF(x_o, y_e):";
			this->set_zero_clmem_force_device();
			fermion_force_eo_device(x_eo_tmp2, y_eo_tmp, ODD, kappa);
			logger.debug() << "\t\t\tF(x_o, y_o):";
			this->set_zero_clmem_force_device();
			fermion_force_eo_device(x_eo_tmp2, y_eo_tmp2, EVEN, kappa);
			logger.debug() << "\t\t\tF(y_o, x_e):";
			this->set_zero_clmem_force_device();
			fermion_force_eo_device(y_eo_tmp2, x_eo_tmp, ODD, kappa);
			logger.debug() << "\t\t\tF(y_o, x_o):";
			this->set_zero_clmem_force_device();
			fermion_force_eo_device(y_eo_tmp2, x_eo_tmp2, EVEN, kappa);
			logger.debug() << "\t\t\tF(y_o, y_e):";
			this->set_zero_clmem_force_device();
			fermion_force_eo_device(y_eo_tmp2, y_eo_tmp, ODD, kappa);
			this->set_zero_clmem_force_device();
		}

		if(logger.beDebug() && debug_hard) {
			//free fields from hard debugging
			int clerr = clReleaseMemObject(y_eo_tmp);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clerr = clReleaseMemObject(y_eo_tmp2);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clerr = clReleaseMemObject(x_eo_tmp);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			//clerr = clReleaseMemObject(x_eo_tmp2);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		}

		if(logger.beDebug() && debug_hard) {
			int clerr = CL_SUCCESS;
			hmc_complex scalprod;
			logger.debug() << "\t\t\tperform some tests on non-eo force ingredients:";
			print_info_inv_field(x_tmp, false, "\t\t\tx:\t");
			print_info_inv_field(y_tmp, false, "\t\t\ty:\t");

			//calculate scalar product
			cl_mem scal_tmp = create_rw_buffer(sizeof(hmc_complex));
			this->set_complex_to_scalar_product_device(x_tmp, y_tmp, scal_tmp);
			get_buffer_from_device(scal_tmp, &scalprod, sizeof(hmc_complex));
			logger.debug() << std::scientific << "\t\t\t(x,y):\t" << scalprod.re << "\t" << scalprod.im;
			clerr = clReleaseMemObject(scal_tmp);

			//calculate difference
			this->saxpy_device(y_tmp, x_tmp, get_clmem_one(), get_clmem_inout());
			print_info_inv_field(get_clmem_inout(), false, "\t\t\t(y-x):\t");

			//calculate non-eo force
			this->set_zero_clmem_force_device();
			copy_buffer_on_device(x_tmp, clmem_phi_inv, meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor));
			copy_buffer_on_device(y_tmp, get_clmem_inout(), meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor));

			print_info_inv_field(clmem_phi_inv, false, "\t\t\tx before:\t");
			print_info_inv_field(get_clmem_inout(), false, "\t\t\ty before:\t");

			fermion_force_device(clmem_phi_inv, get_clmem_inout(), kappa);

			print_info_inv_field(clmem_phi_inv, false, "\t\t\tx after:\t");
			print_info_inv_field(get_clmem_inout(), false, "\t\t\ty after:\t");

			//free fields
			clerr = clReleaseMemObject(x_tmp);
			clerr = clReleaseMemObject(y_tmp);
		}
	} else {
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			/**
			 * The first inversion calculates
			 * X = phi = (Qplusminus)^-1 psi
			 * out of
			 * Qplusminus phi = psi
			 */
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			/**
			 * Trial solution for the spinorfield
			 */
			set_spinorfield_cold_device(get_clmem_inout());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			converged = Opencl_Module_Fermions::cg(::QplusQminus(this), this->get_clmem_inout(), this->get_clmem_phi(), this->clmem_new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions1 += abs(converged);

			/**
			 * Y is now just
			 *  Y = (Qminus) X = (Qminus) (Qplusminus)^-1 psi =
			 *    = (Qplus)^-1 psi
			 */
			Opencl_Module_Fermions::Qminus(this->get_clmem_inout(), clmem_phi_inv, this->clmem_new_u, kappa, mubar);

		} else  {
			logger.debug() << "\t\tcalc fermion force ingredients using bicgstab without eo";

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

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			converged = Opencl_Module_Fermions::bicgstab(::Qplus(this), this->get_clmem_inout(), this->get_clmem_phi(), this->clmem_new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions1 += abs(converged);

			//store this result in clmem_phi_inv
			copy_buffer_on_device(get_clmem_inout(), clmem_phi_inv, sizeof(spinor) * get_spinorfieldsize(get_parameters()));

			/**
			 * Now, one has to calculate
			 * X = (Qminus)^-1 Y = (Qminus)^-1 (Qplus)^-1 psi = (QplusQminus)^-1 psi ??
			 * out of
			 * Qminus clmem_inout = clmem_phi_inv
			 * Therefore, when calculating the final energy of the spinorfield,
			 *  one can just take clmem_phi_inv (see also above)!!
			 */

			//copy former solution to clmem_source
			copy_buffer_on_device(get_clmem_inout(), get_clmem_source(), sizeof(spinor) * get_spinorfieldsize(get_parameters()));
			logger.debug() << "\t\t\tstart solver";

			//this sets clmem_inout cold as trial-solution
			set_spinorfield_cold_device(get_clmem_inout());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::bicgstab(::Qminus(this), get_clmem_inout(), get_clmem_source(), clmem_new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions1 += abs(converged);

		}
		if(logger.beDebug()) {
			print_info_inv_field(clmem_phi_inv, false, "\tY ");
			print_info_inv_field(get_clmem_inout(), false, "\tX ");
		}
		logger.debug() << "\t\tcalc fermion_force...";
		fermion_force_device(clmem_phi_inv, get_clmem_inout(), kappa);
	}
}

void Opencl_Module_Hmc::calc_fermion_force_detratio(usetimer * solvertimer)
{
	/**
	 * For detratio = det(kappa, mubar) / det(kappa2, mubar2) = det(Q_1^+Q_1^-) / det(Q_2^+Q_2^-)
	 * the force has almost the same ingredients as in the above case:
	 *   F(detratio) = - ( - phi^+ deriv(Q_2) X + Y^+ deriv(Q_1) X  ) + h.c.;
	 * where deriv(Q_i) is the same fct. as above with different parameters,
	 * X = (Q_1^+ Q_1^-)^-1 Q_2^+ phi
	 * Y = (Q_1^+)^-1 Q_2^+ phi
	 * (In the case of Q_2 = 1 = const., one recovers the expression for the "normal" force
	 * The main differences are:
	 *   - invert Q_2^+ phi, not phi for X and Y
	 *   - one additional force term with different mass-parameters and without Y
	 */
	int converged = -1;
	hmc_float kappa = get_parameters().get_kappa();
	hmc_float mubar = get_mubar(get_parameters());
	hmc_float kappa2 = get_parameters().get_kappa_mp();
	hmc_float mubar2 = meta::get_mubar_mp(get_parameters());
	if(get_parameters().get_use_eo() == true) {
		//CP: Init tmp spinorfield
		int spinorfield_size = sizeof(spinor) * meta::get_eoprec_spinorfieldsize(get_parameters());
		cl_mem sf_eo_tmp;
		sf_eo_tmp = create_rw_buffer(spinorfield_size);
		//the source is now Q_2^+ phi = sf_eo_tmp
		Opencl_Module_Fermions::Qplus_eo (this->get_clmem_phi_eo(), sf_eo_tmp , get_gaugefield(), kappa2, mubar2);
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			/**
			 * The first inversion calculates
			 * X_even = phi = (Qplusminus_eo)^-1 sf_eo_tmp = (Qplusminus_eo)^-1 Q_2^+ phi
			 * out of
			 * Qplusminus_eo phi = sf_eo_tmp
			 */
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			/**
			 * Trial solution for the spinorfield
			 */
			set_eoprec_spinorfield_cold_device(get_clmem_inout_eo());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\t\t\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::cg_eo(::QplusQminus_eo(this), this->get_clmem_inout_eo(), sf_eo_tmp, this->clmem_new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			//add number of inverter iterations to counter
			*inversions1 += abs(converged);
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\t\t\tinv. field after inversion ");

			/**
			 * Y_even is now just
			 *  Y_even = (Qminus_eo) X_even = (Qminus_eo) (Qplusminus_eo)^-1 sf_eo_tmp =
			 *    = (Qplus_eo)^-1 Q_2^+ psi
			 */
			Opencl_Module_Fermions::Qminus_eo(this->get_clmem_inout_eo(), clmem_phi_inv_eo, this->clmem_new_u, kappa, mubar);
		} else {
			///@todo if wanted, solvertimer has to be used here..
			//logger.debug() << "\t\tcalc fermion force ingredients using bicgstab with eo.";
			/**
			 * The first inversion calculates
			 * Y_even = phi = (Qplus_eo)^-1 psi
			 * out of
			 * Qplus_eo phi = psi
			 */
			//logger.debug() << "\t\tcalc Y_even...";
			//logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			/**
			 * Trial solution for the spinorfield
			 */
			set_zero_spinorfield_eoprec_device(get_clmem_inout_eo());
			gamma5_eo_device(get_clmem_inout_eo());

			//if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\t\t\tinv. field before inversion ");
			//if(logger.beDebug()) print_info_inv_field(get_clmem_phi_eo(), true, "\t\t\tsource before inversion ");
			converged = Opencl_Module_Fermions::bicgstab_eo(::Qplus_eo(this), this->get_clmem_inout_eo(), sf_eo_tmp, this->clmem_new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			//if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\t\t\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions1 += abs(converged);

			//store this result in clmem_phi_inv
			copy_buffer_on_device(get_clmem_inout_eo(), clmem_phi_inv_eo, get_eoprec_spinorfield_buffer_size());

			/**
			 * Now, one has to calculate
			 * X = (Qminus_eo)^-1 Y = (Qminus_eo)^-1 (Qplus_eo)^-1 sf_eo_tmp = (QplusQminus_eo)^-1 sf_eo_tmp = (QplusQminus_eo)^-1 Q_2^+ psi
			 * out of
			 * Qminus_eo clmem_inout_eo = clmem_phi_inv_eo
			 */

			//logger.debug() << "\t\tcalc X_even...";
			//copy former solution to clmem_source
			copy_buffer_on_device(get_clmem_inout_eo(), get_clmem_source_even(), get_eoprec_spinorfield_buffer_size());

			//this sets clmem_inout cold as trial-solution
			set_eoprec_spinorfield_cold_device(get_clmem_inout_eo());

			//if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\t\t\tinv. field before inversion ");
			//if(logger.beDebug()) print_info_inv_field(get_clmem_source_even(), true, "\t\t\tsource before inversion ");
			converged = Opencl_Module_Fermions::bicgstab_eo(::Qminus_eo(this), get_clmem_inout_eo(), get_clmem_source_even(), clmem_new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			//if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\t\t\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions1 += abs(converged);
		}
		/**
		 * At this point, one has to calculate X_odd and Y_odd.
		 * If one has a fermionmatrix
		 *  M = R + D
		 * these are:
		 *  X_odd = -R(-mu)_inv D X_even
		 *  Y_odd = -R(mu)_inv D Y_even
		 */

		/** @fixme below usages of dslash should work, but only because we use bicgstab above
		proper implementation needs to make sure this is always the case */

		///@NOTE the following calculations could also go in a new function for convenience
		//calculate X_odd
		//therefore, clmem_tmp_eo_1 is used as intermediate state. The result is saved in clmem_inout, since
		//  this is used as a default in the force-function.
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			dslash_eo_device(get_clmem_inout_eo(), get_clmem_tmp_eo_1(), clmem_new_u, ODD);
			sax_eoprec_device(get_clmem_tmp_eo_1(), get_clmem_minusone(), get_clmem_tmp_eo_1());
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			dslash_eo_device(get_clmem_inout_eo(), get_clmem_tmp_eo_1(), clmem_new_u, ODD);
			M_tm_inverse_sitediagonal_minus_device(get_clmem_tmp_eo_1(), get_clmem_tmp_eo_2());
			sax_eoprec_device(get_clmem_tmp_eo_2(), get_clmem_minusone(), get_clmem_tmp_eo_1());
		}

		//logger.debug() << "\t\tcalc eo fermion_force F(Y_even, X_odd)...";
		//Calc F(Y_even, X_odd) = F(clmem_phi_inv_eo, clmem_tmp_eo_1)
		fermion_force_eo_device(clmem_phi_inv_eo,  get_clmem_tmp_eo_1(), EVEN, kappa);

		//calculate Y_odd
		//therefore, clmem_tmp_eo_1 is used as intermediate state. The result is saved in clmem_phi_inv, since
		//  this is used as a default in the force-function.
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			dslash_eo_device(clmem_phi_inv_eo, get_clmem_tmp_eo_1(), clmem_new_u, ODD);
			sax_eoprec_device(get_clmem_tmp_eo_1(), get_clmem_minusone(), get_clmem_tmp_eo_1());
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			dslash_eo_device(clmem_phi_inv_eo, get_clmem_tmp_eo_1(), clmem_new_u, ODD);
			M_tm_inverse_sitediagonal_device(get_clmem_tmp_eo_1(), get_clmem_tmp_eo_2());
			sax_eoprec_device(get_clmem_tmp_eo_2(), get_clmem_minusone(), get_clmem_tmp_eo_1());
		}

		//logger.debug() << "\t\tcalc eoprec fermion_force F(Y_odd, X_even)...";
		//Calc F(Y_odd, X_even) = F(clmem_tmp_eo_1, clmem_inout_eo)
		fermion_force_eo_device(get_clmem_tmp_eo_1(), get_clmem_inout_eo(), ODD, kappa);

		/**
		 *Now, one has the additional term - phi^+ deriv(Q_2) X
		 *This works exactly the same as above, except replacing Y -> -phi and different mass parameters
		 */

		//Y is not needed anymore, therefore use clmem_phi_inv_eo to store -phi
		sax_eoprec_device(get_clmem_phi_eo(), get_clmem_minusone(), clmem_phi_inv_eo);

		//logger.debug() << "\t\tcalc eo fermion_force F(Y_even, X_odd)...";
		//Calc F(Y_even, X_odd) = F(clmem_phi_inv_eo, clmem_tmp_eo_1)
		fermion_force_eo_device(clmem_phi_inv_eo,  get_clmem_tmp_eo_1(), EVEN, kappa2);

		//calculate phi_odd
		//this works in the same way as with Y above, since -phi_even is saved in the same buffer as Y_even
		//therefore, clmem_tmp_eo_1 is used as intermediate state. The result is saved in clmem_phi_inv, since
		//  this is used as a default in the force-function.
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			dslash_eo_device(clmem_phi_inv_eo, get_clmem_tmp_eo_1(), clmem_new_u, ODD);
			sax_eoprec_device(get_clmem_tmp_eo_1(), get_clmem_minusone(), get_clmem_tmp_eo_1());
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			dslash_eo_device(clmem_phi_inv_eo, get_clmem_tmp_eo_1(), clmem_new_u, ODD);
			M_tm_inverse_sitediagonal_device(get_clmem_tmp_eo_1(), get_clmem_tmp_eo_2());
			sax_eoprec_device(get_clmem_tmp_eo_2(), get_clmem_minusone(), get_clmem_tmp_eo_1());
		}

		//logger.debug() << "\t\tcalc eoprec fermion_force F(Y_odd, X_even)...";
		//Calc F(Y_odd, X_even) = F(clmem_tmp_eo_1, clmem_inout_eo)
		fermion_force_eo_device(get_clmem_tmp_eo_1(), get_clmem_inout_eo(), ODD, kappa2);


		int clerr = clReleaseMemObject(sf_eo_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	} else {
		//CP: Init tmp spinorfield
		int spinorfield_size = sizeof(spinor) * get_spinorfieldsize(get_parameters());
		cl_mem sf_tmp;
		sf_tmp = create_rw_buffer(spinorfield_size);
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		//the source is now Q_2^+ phi = sf_tmp
		Opencl_Module_Fermions::Qplus (this->get_clmem_phi(), sf_tmp , get_gaugefield(), kappa2, mubar2);
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			/**
			 * The first inversion calculates
			 * X = phi = (Qplusminus)^-1 sf_tmp
			 * out of
			 * Qplusminus phi = sf_tmp
			 */
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			/**
			 * Trial solution for the spinorfield
			 */
			set_spinorfield_cold_device(get_clmem_inout());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			converged = Opencl_Module_Fermions::cg(::QplusQminus(this), this->get_clmem_inout(), sf_tmp, this->clmem_new_u, get_parameters().get_force_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions1 += abs(converged);

			/**
			 * Y is now just
			 *  Y = (Qminus) X = (Qminus) (Qplusminus)^-1 sf_tmp =
			 *    = (Qplus)^-1 sf_tmp
			 */
			Opencl_Module_Fermions::Qminus(this->get_clmem_inout(), clmem_phi_inv, this->clmem_new_u);

		} else  {
			logger.debug() << "\t\tcalc fermion force ingredients using bicgstab without eo";

			/**
			 * The first inversion calculates
			 * Y = phi = (Qplus)^-1 sf_tmp
			 * out of
			 * Qplus phi = sf_tmp
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

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			converged = Opencl_Module_Fermions::bicgstab(::Qplus(this), this->get_clmem_inout(), sf_tmp, this->clmem_new_u, get_parameters().get_force_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions1 += abs(converged);

			//store this result in clmem_phi_inv
			copy_buffer_on_device(get_clmem_inout(), clmem_phi_inv, sizeof(spinor) * get_spinorfieldsize(get_parameters()));

			/**
			 * Now, one has to calculate
			 * X = (Qminus)^-1 Y = (Qminus)^-1 (Qplus)^-1 psi = (QplusQminus)^-1 psi ??
			 * out of
			 * Qminus clmem_inout = clmem_phi_inv
			 * Therefore, when calculating the final energy of the spinorfield,
			 *  one can just take clmem_phi_inv (see also above)!!
			 */

			//copy former solution to clmem_source
			copy_buffer_on_device(get_clmem_inout(), get_clmem_source(), sizeof(spinor) * get_spinorfieldsize(get_parameters()));
			logger.debug() << "\t\t\tstart solver";

			//this sets clmem_inout cold as trial-solution
			set_spinorfield_cold_device(get_clmem_inout());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::bicgstab(::Qminus(this), get_clmem_inout(), get_clmem_source(), clmem_new_u, get_parameters().get_force_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions1 += abs(converged);

		}
		if(logger.beDebug()) {
			print_info_inv_field(clmem_phi_inv, false, "\tY ");
			print_info_inv_field(get_clmem_inout(), false, "\tX ");
		}
		logger.debug() << "\t\tcalc fermion_force...";
		fermion_force_device(clmem_phi_inv, get_clmem_inout(), kappa);

		/**
		 *Now, one has the additional term - phi^+ deriv(Q_2) X
		 *This works exactly the same as above, except replacing Y -> -phi and different mass parameters
		 */

		//Y is not needed anymore, therefore use clmem_phi_inv_eo to store -phi
		sax_device(get_clmem_phi(), get_clmem_minusone(), clmem_phi_inv);

		fermion_force_device(clmem_phi_inv, get_clmem_inout(), kappa2);

		int clerr = clReleaseMemObject(sf_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}
}

void Opencl_Module_Hmc::calc_gauge_force()
{
	logger.debug() << "\t\tcalc gauge_force...";
	gauge_force_device();
	if(meta::get_use_rectangles(get_parameters()) == true) {
		logger.debug() << "\t\tcalc rect gauge_force...";
		gauge_force_tlsym_device();
	}
}

hmc_float Opencl_Module_Hmc::calc_s_fermion()
{
	logger.debug() << "calc final fermion energy...";
	//this function essentially performs the same steps as in the force-calculation, but with higher precision.
	//  therefore, comments are deleted here...
	//  Furthermore, in the bicgstab-case, the second inversions are not needed
	int converged = -1;
	if(get_parameters().get_use_eo() == true) {
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			logger.debug() << "\t\t\tstart solver";

			set_eoprec_spinorfield_cold_device(get_clmem_inout_eo());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::cg_eo(::QplusQminus_eo(this), this->get_clmem_inout_eo(), this->get_clmem_phi_eo(), this->clmem_new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions0 += abs(converged);

			Opencl_Module_Fermions::Qminus_eo(this->get_clmem_inout_eo(), clmem_phi_inv_eo, this->clmem_new_u);
		} else {
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			set_zero_spinorfield_eoprec_device(get_clmem_inout_eo());
			gamma5_eo_device(get_clmem_inout_eo());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\tinv. field before inversion ");
			if(logger.beDebug()) print_info_inv_field(get_clmem_phi_eo(), true, "\tsource before inversion ");
			converged = Opencl_Module_Fermions::bicgstab_eo(::Qplus_eo(this), this->get_clmem_inout_eo(), this->get_clmem_phi_eo(), this->clmem_new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions0 += abs(converged);

			//store this result in clmem_phi_inv
			copy_buffer_on_device(get_clmem_inout_eo(), clmem_phi_inv_eo, get_eoprec_spinorfield_buffer_size());

		}
	} else {
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			set_spinorfield_cold_device(get_clmem_inout());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::cg(::QplusQminus(this), this->get_clmem_inout(), this->get_clmem_phi(), this->clmem_new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions0 += abs(converged);

			Opencl_Module_Fermions::Qminus(this->get_clmem_inout(), clmem_phi_inv, this->clmem_new_u);

		} else  {

			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			set_spinorfield_cold_device(get_clmem_inout());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::bicgstab(::Qplus(this), this->get_clmem_inout(), this->get_clmem_phi(), this->clmem_new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions0 += abs(converged);

			//store this result in clmem_phi_inv
			copy_buffer_on_device(get_clmem_inout(), clmem_phi_inv, sizeof(spinor) * get_spinorfieldsize(get_parameters()));
		}
	}
	///@todo: this can be moved in the ifs above!!
	if(get_parameters().get_use_eo() == true) {
		set_float_to_global_squarenorm_eoprec_device(clmem_phi_inv_eo, clmem_s_fermion);
	} else {
		set_float_to_global_squarenorm_device(clmem_phi_inv, clmem_s_fermion);
	}
	hmc_float tmp;
	get_buffer_from_device(clmem_s_fermion, &tmp, sizeof(hmc_float));
	return tmp;
}

hmc_float Opencl_Module_Hmc::calc_s_fermion_mp()
{
	logger.debug() << "calc final fermion energy...";
	//this function essentially performs the same steps as in the non mass-prec case, however, one has to apply one more matrix multiplication
	//  therefore, comments are deleted here...
	//  Furthermore, in the bicgstab-case, the second inversions are not needed
	int converged = -1;
	if(get_parameters().get_use_eo() == true) {
		//CP: Init tmp spinorfield
		int spinorfield_size = sizeof(spinor) * meta::get_eoprec_spinorfieldsize(get_parameters());
		cl_mem sf_eo_tmp;
		sf_eo_tmp = create_rw_buffer(spinorfield_size);

		//sf_eo_tmp = Qplus_eo(heavy_mass) phi_mp_eo
		Opencl_Module_Fermions::Qplus_eo (this->get_clmem_phi_mp_eo(), sf_eo_tmp , get_gaugefield(), get_parameters().get_kappa_mp(), meta::get_mubar_mp(get_parameters()));
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			logger.debug() << "\t\t\tstart solver";
			set_eoprec_spinorfield_cold_device(get_clmem_inout_eo());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::cg_eo(::QplusQminus_eo(this), this->get_clmem_inout_eo(), sf_eo_tmp, this->clmem_new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions0 += abs(converged);

			Opencl_Module_Fermions::Qminus_eo(this->get_clmem_inout_eo(), clmem_phi_inv_eo, this->clmem_new_u);
		} else {
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			set_zero_spinorfield_eoprec_device(get_clmem_inout_eo());
			gamma5_eo_device(get_clmem_inout_eo());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\tinv. field before inversion ");
			if(logger.beDebug()) print_info_inv_field(sf_eo_tmp, true, "\tsource before inversion ");
			converged = Opencl_Module_Fermions::bicgstab_eo(::Qplus_eo(this), this->get_clmem_inout_eo(), sf_eo_tmp, this->clmem_new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout_eo(), true, "\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions0 += abs(converged);

			//store this result in clmem_phi_inv
			copy_buffer_on_device(get_clmem_inout_eo(), clmem_phi_inv_eo, get_eoprec_spinorfield_buffer_size());
		}
		int clerr = clReleaseMemObject(sf_eo_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	} else {
		//CP: Init tmp spinorfield
		int spinorfield_size = sizeof(spinor) * get_spinorfieldsize(get_parameters());
		cl_mem sf_tmp;
		sf_tmp = create_rw_buffer(spinorfield_size);

		//sf_eo_tmp = Qplus(light_mass) phi_mp
		Opencl_Module_Fermions::Qplus (this->get_clmem_phi_mp(), sf_tmp , get_gaugefield(), get_parameters().get_kappa_mp(), meta::get_mubar_mp(get_parameters()));
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			set_spinorfield_cold_device(get_clmem_inout());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::cg(::QplusQminus(this), this->get_clmem_inout(), sf_tmp, this->clmem_new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions0 += abs(converged);

			Opencl_Module_Fermions::Qminus(this->get_clmem_inout(), clmem_phi_inv, this->clmem_new_u);

		} else  {

			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			set_spinorfield_cold_device(get_clmem_inout());

			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field before inversion ");
			converged = Opencl_Module_Fermions::bicgstab(::Qplus(this), this->get_clmem_inout(), sf_tmp, this->clmem_new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) print_info_inv_field(get_clmem_inout(), false, "\tinv. field after inversion ");
			//add number of inverter iterations to counter
			*inversions0 += abs(converged);

			//store this result in clmem_phi_inv
			copy_buffer_on_device(get_clmem_inout(), clmem_phi_inv, sizeof(spinor) * get_spinorfieldsize(get_parameters()));
		}
		int clerr = clReleaseMemObject(sf_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}
	///@todo: this can be moved in the ifs above!!
	if(get_parameters().get_use_eo() == true) {
		set_float_to_global_squarenorm_eoprec_device(clmem_phi_inv_eo, clmem_s_fermion);
	} else {
		set_float_to_global_squarenorm_device(clmem_phi_inv, clmem_s_fermion);
	}
	hmc_float tmp;
	get_buffer_from_device(clmem_s_fermion, &tmp, sizeof(hmc_float));
	return tmp;
}

hmc_observables Opencl_Module_Hmc::metropolis(hmc_float rnd, hmc_float beta)
{
	//Calc Hamiltonian
	logger.debug() << "Calculate Hamiltonian";
	hmc_float deltaH = 0.;
	hmc_float s_old = 0.;
	hmc_float s_new = 0.;

	//Gauge-Part
	hmc_float tplaq, splaq, plaq;
	hmc_float tplaq_new, splaq_new, plaq_new;
	hmc_float rect_new = 0.;
	hmc_float rect = 0.;
	hmc_complex poly;
	hmc_complex poly_new;
	//In this call, the observables are calculated already with appropiate Weighting factor of 2.0/(VOL4D*NDIM*(NDIM-1)*NC)
	Opencl_Module::gaugeobservables(get_gaugefield(), &plaq,  &tplaq, &splaq, &poly);
	Opencl_Module::gaugeobservables(clmem_new_u, &plaq_new,  &tplaq_new, &splaq_new, &poly_new);
	//plaq has to be divided by the norm-factor to get s_gauge
	hmc_float factor = 1. / (meta::get_plaq_norm(get_parameters()));
	if(meta::get_use_rectangles(get_parameters()) == true) {
		Opencl_Module::gaugeobservables_rectangles(get_gaugefield(), &rect);
		Opencl_Module::gaugeobservables_rectangles(clmem_new_u, &rect_new);
		hmc_float c0 = meta::get_c0(get_parameters());
		hmc_float c1 = meta::get_c1(get_parameters());
		deltaH = - beta * ( c0 * (plaq - plaq_new) / factor + c1 * ( rect - rect_new )  );
		s_old = - beta * ( c0 * (plaq) / factor + c1 * ( rect )  );
		s_new = - beta * ( c0 * (plaq_new) / factor + c1 * ( rect_new )  );

	} else {
		/** NOTE: the minus here is introduced to fit tmlqcd!!! */
		deltaH = -(plaq - plaq_new) * beta / factor;
		s_old = -(plaq ) * beta / factor;
		s_new = -(plaq_new) * beta / factor;
	}

	logger.debug() << "\tS_gauge(old field) = " << setprecision(10) << s_old;
	logger.debug() << "\tS_gauge(new field) = " << setprecision(10) << s_new;
	logger.info() << "\tdeltaS_gauge = " << setprecision(10) << deltaH;

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
	logger.info() << "\tdeltaS_gaugemom = " << setprecision(10) << 0.5 * (p2 - new_p2);

	//Fermion-Part:
	if(! get_parameters().get_use_gauge_only() ) {
		hmc_float spinor_energy_init, s_fermion_final;
		//initial energy has been computed in the beginning...
		Opencl_Module_Hmc::get_buffer_from_device(clmem_s_fermion_init, &spinor_energy_init, sizeof(hmc_float));
		// sum_links phi*_i (M^+M)_ij^-1 phi_j
		s_fermion_final = calc_s_fermion();
		deltaH += spinor_energy_init - s_fermion_final;

		logger.debug() << "\tS_ferm(old field) = " << setprecision(10) <<  spinor_energy_init;
		logger.debug() << "\tS_ferm(new field) = " << setprecision(10) << s_fermion_final;
		logger.info() << "\tdeltaS_ferm = " << spinor_energy_init - s_fermion_final;
		if( get_parameters().get_use_mp() ) {
			hmc_float spinor_energy_mp_init, s_fermion_mp_final;
			//initial energy has been computed in the beginning...
			Opencl_Module_Hmc::get_buffer_from_device(clmem_s_fermion_mp_init, &spinor_energy_mp_init, sizeof(hmc_float));
			// sum_links phi*_i (M^+M)_ij^-1 phi_j
			s_fermion_mp_final = calc_s_fermion_mp();
			deltaH += spinor_energy_mp_init - s_fermion_mp_final;

			logger.debug() << "\tS_ferm_mp(old field) = " << setprecision(10) <<  spinor_energy_mp_init;
			logger.debug() << "\tS_ferm_mp(new field) = " << setprecision(10) << s_fermion_mp_final;
			logger.info() << "\tdeltaS_ferm_mp = " << spinor_energy_init - s_fermion_mp_final;
		}
	}
	//Metropolis-Part
	hmc_float compare_prob;
	if(deltaH < 0) {
		compare_prob = exp(deltaH);
	} else {
		compare_prob = 1.0;
	}
	logger.info() << "\tdeltaH = " << deltaH << "\tAcc-Prop = " << compare_prob;
	hmc_observables tmp;
	if(rnd <= compare_prob) {
		tmp.accept = 1;
		tmp.plaq = plaq_new;
		tmp.tplaq = tplaq_new;
		tmp.splaq = splaq_new;
		tmp.poly = poly_new;
		tmp.deltaH = deltaH;
		tmp.prob = compare_prob;
		if(meta::get_use_rectangles(get_parameters()) ) tmp.rectangles = rect_new / get_rect_norm(get_parameters());
	} else {
		tmp.accept = 0;
		tmp.plaq = plaq;
		tmp.tplaq = tplaq;
		tmp.splaq = splaq;
		tmp.poly = poly;
		tmp.deltaH = deltaH;
		tmp.prob = compare_prob;
		if(meta::get_use_rectangles(get_parameters()) ) tmp.rectangles = rect / get_rect_norm(get_parameters());
	}

	return tmp;
}

void Opencl_Module_Hmc::calc_spinorfield_init_energy(cl_mem dest)
{
	//Suppose the initial spinorfield is saved in phi_inv
	//  it is created in generate_gaussian_spinorfield_device
	if(get_parameters().get_use_eo() == true) {
		Opencl_Module_Fermions::set_float_to_global_squarenorm_eoprec_device(clmem_phi_inv_eo, dest);
	} else {
		Opencl_Module_Fermions::set_float_to_global_squarenorm_device(clmem_phi_inv, dest);
	}
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

	enqueueKernel(md_update_gaugemomenta , gs2, ls2);

	if(logger.beDebug()) {
		cl_mem force_tmp = create_rw_buffer(sizeof(hmc_float));
		hmc_float resid;
		this->set_float_to_gaugemomentum_squarenorm_device(clmem_new_p, force_tmp);
		get_buffer_from_device(force_tmp, &resid, sizeof(hmc_float));
		logger.debug() <<  "\tupdated gaugemomenta energy:\t" << resid;
		int clerr = clReleaseMemObject(force_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		if(resid != resid) {
			throw Print_Error_Message("calculation of gm gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
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
	this->get_work_sizes(set_zero_gaugemomentum, this->get_device_type(), &ls2, &gs2, &num_groups);
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

	if(logger.beDebug()) {
		cl_mem gauge_force_tmp = create_rw_buffer(sizeof(hmc_float));
		hmc_float gauge_force_energy = 0.;
		this->set_float_to_gaugemomentum_squarenorm_device(clmem_force, gauge_force_tmp);
		get_buffer_from_device(gauge_force_tmp, &gauge_force_energy, sizeof(hmc_float));

		//logger.debug() <<  "\t\t\tgauge force:\t" << gauge_force_energy;

		int clerr = clReleaseMemObject(gauge_force_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		if(gauge_force_energy != gauge_force_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

	//recalculate force with local buffer to get only this contribution to the force vector
	if(logger.beDebug()) {
		int gaugemomentum_size = get_gaugemomentum_buffer_size();
		cl_mem force2;
		force2 = create_rw_buffer(gaugemomentum_size);
		//init new buffer to zero
		this->get_work_sizes(md_update_gaugemomenta, this->get_device_type(), &ls2, &gs2, &num_groups);
		clerr = clSetKernelArg(set_zero_gaugemomentum, 0, sizeof(cl_mem), &force2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel( set_zero_gaugemomentum , gs2, ls2);

		//re-calculate force
		clerr = clSetKernelArg(gauge_force, 1, sizeof(cl_mem), &force2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		enqueueKernel( gauge_force , gs2, ls2);

		cl_mem check_force_tmp = create_rw_buffer(sizeof(hmc_float));
		hmc_float check_force_energy = 0.;
		this->set_float_to_gaugemomentum_squarenorm_device(force2, check_force_tmp);
		get_buffer_from_device(check_force_tmp, &check_force_energy, sizeof(hmc_float));
		//logger.debug() <<  "\t\t\t\tforce contribution:\t" << check_force_energy;
		int clerr = clReleaseMemObject(check_force_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		if(check_force_energy != check_force_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
		clerr = clReleaseMemObject(force2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}

}

void Opencl_Module_Hmc::gauge_force_device(cl_mem gf, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gauge_force, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gauge_force, 0, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(gauge_force, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( gauge_force , gs2, ls2);
}


void Opencl_Module_Hmc::gauge_force_tlsym_device()
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gauge_force_tlsym, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gauge_force_tlsym, 0, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(gauge_force_tlsym, 1, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( gauge_force_tlsym , gs2, ls2);

	if(logger.beDebug()) {
		cl_mem gauge_force_tlsym_tmp = create_rw_buffer(sizeof(hmc_float));
		hmc_float gauge_force_tlsym_energy = 0.;
		this->set_float_to_gaugemomentum_squarenorm_device(clmem_force, gauge_force_tlsym_tmp);
		get_buffer_from_device(gauge_force_tlsym_tmp, &gauge_force_tlsym_energy, sizeof(hmc_float));

		logger.debug() <<  "\t\t\tgauge force tlsym:\t" << gauge_force_tlsym_energy;

		int clerr = clReleaseMemObject(gauge_force_tlsym_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		if(gauge_force_tlsym_energy != gauge_force_tlsym_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

	//recalculate force with local buffer to get only this contribution to the force vector
	if(logger.beDebug()) {
		int gaugemomentum_size = get_gaugemomentum_buffer_size();
		cl_mem force2;
		force2 = create_rw_buffer(gaugemomentum_size);
		//init new buffer to zero
		this->get_work_sizes(md_update_gaugemomenta, this->get_device_type(), &ls2, &gs2, &num_groups);
		clerr = clSetKernelArg(set_zero_gaugemomentum, 0, sizeof(cl_mem), &force2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel( set_zero_gaugemomentum , gs2, ls2);

		//re-calculate force
		clerr = clSetKernelArg(gauge_force_tlsym, 1, sizeof(cl_mem), &force2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		enqueueKernel( gauge_force_tlsym , gs2, ls2);

		cl_mem check_force_tmp = create_rw_buffer(sizeof(hmc_float));
		hmc_float check_force_energy = 0.;
		this->set_float_to_gaugemomentum_squarenorm_device(force2, check_force_tmp);
		get_buffer_from_device(check_force_tmp, &check_force_energy, sizeof(hmc_float));
		logger.debug() <<  "\t\t\t\tforce contribution:\t" << check_force_energy;
		int clerr = clReleaseMemObject(check_force_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		if(check_force_energy != check_force_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
		clerr = clReleaseMemObject(force2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}

}

void Opencl_Module_Hmc::fermion_force_device(cl_mem Y, cl_mem X, hmc_float kappa)
{
	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = get_parameters().get_kappa();
	else kappa_tmp = kappa;

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(fermion_force, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//fermion_force(field, Y, X, out);
	int clerr = clSetKernelArg(fermion_force, 0, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 1, sizeof(cl_mem), &Y);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 2, sizeof(cl_mem), &X);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 3, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 4, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( fermion_force , gs2, ls2);

	if(logger.beDebug()) {
		cl_mem noneo_force_tmp = create_rw_buffer(sizeof(hmc_float));
		hmc_float noneo_force_energy = 0.;
		this->set_float_to_gaugemomentum_squarenorm_device(clmem_force, noneo_force_tmp);
		get_buffer_from_device(noneo_force_tmp, &noneo_force_energy, sizeof(hmc_float));
		//logger.debug() <<  "\t\t\tnon-eo force:\t" << noneo_force_energy;
		int clerr = clReleaseMemObject(noneo_force_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		if(noneo_force_energy != noneo_force_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

	//recalculate force with local buffer to get only this contribution to the force vector
	if(logger.beDebug()) {
		int gaugemomentum_size = get_gaugemomentum_buffer_size();
		cl_mem force2;
		force2 = create_rw_buffer(gaugemomentum_size);
		//init new buffer to zero
		this->get_work_sizes(md_update_gaugemomenta, this->get_device_type(), &ls2, &gs2, &num_groups);
		clerr = clSetKernelArg(set_zero_gaugemomentum, 0, sizeof(cl_mem), &force2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel( set_zero_gaugemomentum , gs2, ls2);

		//re-calculate force
		clerr = clSetKernelArg(fermion_force, 3, sizeof(cl_mem), &force2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		enqueueKernel( fermion_force , gs2, ls2);

		cl_mem check_force_tmp = create_rw_buffer(sizeof(hmc_float));
		hmc_float check_force_energy = 0.;
		this->set_float_to_gaugemomentum_squarenorm_device(force2, check_force_tmp);
		get_buffer_from_device(check_force_tmp, &check_force_energy, sizeof(hmc_float));
		//logger.debug() <<  "\t\t\t\tforce contribution:\t" << check_force_energy;
		int clerr = clReleaseMemObject(check_force_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		if(check_force_energy != check_force_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
		clerr = clReleaseMemObject(force2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}

}

//the argument kappa is set to ARG_DEF as default
void Opencl_Module_Hmc::fermion_force_eo_device(cl_mem Y, cl_mem X, int evenodd, hmc_float kappa)
{
	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = get_parameters().get_kappa();
	else kappa_tmp = kappa;

	//fermion_force(field, Y, X, out);
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(fermion_force_eo, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(fermion_force_eo, 0, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force_eo, 1, sizeof(cl_mem), &Y);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force_eo, 2, sizeof(cl_mem), &X);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force_eo, 3, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force_eo, 4, sizeof(int), &evenodd);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force_eo, 5, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( fermion_force_eo , gs2, ls2);

	if(logger.beDebug()) {
		cl_mem force_tmp = create_rw_buffer(sizeof(hmc_float));
		hmc_float resid;
		this->set_float_to_gaugemomentum_squarenorm_device(clmem_force, force_tmp);
		get_buffer_from_device(force_tmp, &resid, sizeof(hmc_float));
		//logger.debug() <<  "\t\t\teoprec force:\t" << resid;

		int clerr = clReleaseMemObject(force_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		if(resid != resid) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

	//recalculate force with local buffer, giving only this contribution to the force
	if(logger.beDebug()) {
		int gaugemomentum_size = get_gaugemomentum_buffer_size();
		cl_mem force2;
		force2 = create_rw_buffer(gaugemomentum_size);
		//init new buffer to zero
		this->get_work_sizes(md_update_gaugemomenta, this->get_device_type(), &ls2, &gs2, &num_groups);
		clerr = clSetKernelArg(set_zero_gaugemomentum, 0, sizeof(cl_mem), &force2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel( set_zero_gaugemomentum , gs2, ls2);

		//re-calculate force
		clerr = clSetKernelArg(fermion_force_eo, 3, sizeof(cl_mem), &force2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		enqueueKernel( fermion_force_eo , gs2, ls2);

		cl_mem check_force_tmp = create_rw_buffer(sizeof(hmc_float));
		hmc_float check_force_energy = 0.;
		this->set_float_to_gaugemomentum_squarenorm_device(force2, check_force_tmp);
		get_buffer_from_device(check_force_tmp, &check_force_energy, sizeof(hmc_float));
		//logger.debug() <<  "\t\t\t\tforce contribution:\t" << check_force_energy;
		int clerr = clReleaseMemObject(check_force_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
		if(check_force_energy != check_force_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
		clerr = clReleaseMemObject(force2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}
}

void Opencl_Module_Hmc::stout_smeared_fermion_force_device(cl_mem * gf_intermediate)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(stout_smear_fermion_force, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
}

void Opencl_Module_Hmc::set_float_to_gaugemomentum_squarenorm_device(cl_mem clmem_in, cl_mem out)
{
	//__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gaugemomentum_squarenorm, this->get_device_type(), &ls2, &gs2, &num_groups);

	int global_buf_size_float = sizeof(hmc_float) * num_groups;
	cl_mem  clmem_global_squarenorm_buf_glob = create_rw_buffer(global_buf_size_float);

	int clerr = clSetKernelArg(gaugemomentum_squarenorm, 0, sizeof(cl_mem), &clmem_in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(gaugemomentum_squarenorm, 1, sizeof(cl_mem), &clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(gaugemomentum_squarenorm, 2, sizeof(hmc_float) * ls2, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	enqueueKernel(gaugemomentum_squarenorm  , gs2, ls2);

	clerr = clSetKernelArg(global_squarenorm_reduction, 0, sizeof(cl_mem), &clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(global_squarenorm_reduction, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	enqueueKernel(global_squarenorm_reduction, gs2, ls2);

	clReleaseMemObject(clmem_global_squarenorm_buf_glob);
}

size_t Opencl_Module_Hmc::get_gaugemomentum_buffer_size()
{
	if(!gaugemomentum_buf_size) {
		if(use_soa) {
			gaugemomentum_buf_size = calculateStride(meta::get_vol4d(parameters) * NDIM, sizeof(hmc_float)) * sizeof(ae);
		} else {
			gaugemomentum_buf_size = meta::get_vol4d(parameters) * NDIM * sizeof(ae);
		}
	}
	return gaugemomentum_buf_size;
}

void Opencl_Module_Hmc::importGaugemomentumBuffer(const cl_mem dest, const ae * const data)
{
	cl_int clerr;
	if(use_soa) {
		const size_t aos_size = meta::get_vol4d(parameters) * NDIM * sizeof(ae);
		cl_mem tmp = create_ro_buffer(aos_size);
		clerr = clEnqueueWriteBuffer(get_queue(), tmp, CL_TRUE, 0, aos_size, data, 0, 0, NULL);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);

		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(gaugemomentum_convert_to_soa, this->get_device_type(), &ls2, &gs2, &num_groups);
		clerr = clSetKernelArg(gaugemomentum_convert_to_soa, 0, sizeof(cl_mem), &dest);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gaugemomentum_convert_to_soa, 1, sizeof(cl_mem), &tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel(gaugemomentum_convert_to_soa, gs2, ls2);

		clReleaseMemObject(tmp);
	} else {
		clerr = clEnqueueWriteBuffer(get_queue(), dest, CL_TRUE, 0, get_gaugemomentum_buffer_size(), data, 0, 0, NULL);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);
	}
}

void Opencl_Module_Hmc::exportGaugemomentumBuffer(ae * const dest, const cl_mem buf)
{
	cl_int clerr;
	if(use_soa) {
		const size_t aos_size = meta::get_vol4d(parameters) * NDIM * sizeof(ae);
		cl_mem tmp = create_wo_buffer(aos_size);

		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(gaugemomentum_convert_from_soa, this->get_device_type(), &ls2, &gs2, &num_groups);
		clerr = clSetKernelArg(gaugemomentum_convert_from_soa, 0, sizeof(cl_mem), &tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gaugemomentum_convert_from_soa, 1, sizeof(cl_mem), &buf);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel(gaugemomentum_convert_from_soa, gs2, ls2);

		clerr = clEnqueueReadBuffer(get_queue(), tmp, CL_TRUE, 0, aos_size, dest, 0, 0, NULL);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);

		clReleaseMemObject(tmp);
	} else {
		clerr = clEnqueueReadBuffer(get_queue(), buf, CL_TRUE, 0, get_gaugemomentum_buffer_size(), dest, 0, 0, NULL);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);
	}
}
