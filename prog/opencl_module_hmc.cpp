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

	int clerr = CL_SUCCESS;

	int spinorfield_size = sizeof(spinor) * SPINORFIELDSIZE;
	int gaugemomentum_size = sizeof(ae) * GAUGEMOMENTASIZE2;
	int float_size = sizeof(hmc_float);
	hmc_complex one = hmc_complex_one;
	hmc_complex minusone = hmc_complex_minusone;

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
		S = EOPREC_SPINORFIELDSIZE;
	else
		S = SPINORFIELDSIZE;
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
