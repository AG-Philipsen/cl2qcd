#include "opencl_hmc.h"
#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

hmc_error Opencl_hmc::fill_kernels_file ()
{
	//give a list of all kernel-files
	Opencl_fermions::fill_kernels_file();

	cl_kernels_file.push_back("types_hmc.h");
	cl_kernels_file.push_back("operations_gaugemomentum.cl");
	cl_kernels_file.push_back("operations_force.cl");
 	cl_kernels_file.push_back("molecular_dynamics.cl");

	return HMC_SUCCESS;  
}

hmc_error Opencl_hmc::fill_collect_options(stringstream* collect_options)
{

	Opencl_fermions::fill_collect_options(collect_options);
	*collect_options <<  " -DBETA=" << get_parameters()->get_beta() << " -DGAUGEMOMENTASIZE=" << GAUGEMOMENTASIZE2;
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::fill_buffers()
{
	Opencl_fermions::fill_buffers();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::fill_kernels()
{
	Opencl_fermions::fill_kernels();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::init(cl_device_type wanted_device_type, usetimer* timer, inputparameters* parameters)
{
	hmc_error err = Opencl_fermions::init(wanted_device_type, timer, parameters);
	err |= init_hmc_variables(parameters, timer);
	return err;
}

hmc_error Opencl_hmc::init_hmc_variables(inputparameters* parameters, usetimer * timer)
{
	(*timer).reset();
	
	//CP: this is copied from opencl_fermions
		// decide on work-sizes
	size_t local_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		local_work_size = NUMTHREADS; /// @todo have local work size depend on kernel properties (and device? autotune?)
	else
		local_work_size = 1; // nothing else makes sens on CPU

	size_t global_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		global_work_size = 4 * NUMTHREADS * max_compute_units; /// @todo autotune
	else
		global_work_size = max_compute_units;

	const cl_uint num_groups = (global_work_size + local_work_size - 1) / local_work_size;
	global_work_size = local_work_size * num_groups;

	(*timer).reset();

	logger.trace()<< "init HMC variables...";
	int clerr = CL_SUCCESS;

	int spinorfield_size = sizeof(spinor)*SPINORFIELDSIZE;
	int eoprec_spinorfield_size = sizeof(spinor)*EOPREC_SPINORFIELDSIZE;
	int gaugemomentum_size = sizeof(ae)*GAUGEMOMENTASIZE2;
	int complex_size = sizeof(hmc_complex);
	int float_size = sizeof(hmc_float);
	int global_buf_size = complex_size * num_groups;
	int global_buf_size_float = float_size * num_groups;
	hmc_complex one = hmc_complex_one;
	hmc_complex minusone = hmc_complex_minusone;
	hmc_float tmp;
	

	/** @todo insert variables needed */
	//init mem-objects

	logger.trace() << "Create buffer for phi_inv...";
	clmem_phi_inv = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
	if(clerr!=CL_SUCCESS) {
		cout<<"creating clmem_inout failed, aborting..."<<endl;
		exit(HMC_OCLERROR);
	}
	logger.trace() << "Create buffer for new gaugefield...";
	clmem_new_u = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(s_gaugefield), 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.trace() << "Create buffer for gaugemomentum p...";
	clmem_p = clCreateBuffer(context, CL_MEM_READ_WRITE, gaugemomentum_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.trace() << "Create buffer for gaugemomentum new_p...";
	clmem_new_p = clCreateBuffer(context, CL_MEM_READ_WRITE, gaugemomentum_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.trace() << "Create buffer for initial spinorfield energy...";
	clmem_energy_init = clCreateBuffer(context, CL_MEM_READ_WRITE, float_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.trace() << "Create buffer for deltaH...";
	clmem_deltah = clCreateBuffer(context, CL_MEM_READ_WRITE, float_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	
	//init kernels
	logger.debug() << "Create kernel generate_gaussian_spinorfield...";
	generate_gaussian_spinorfield = clCreateKernel(clprogram, "generate_gaussian_spinorfield", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.debug() << "Create kernel generate_gaussian_gaugemomenta...";
	generate_gaussian_spinorfield = clCreateKernel(clprogram, "generate_gaussian_gaugemomenta", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.debug() << "Create kernel md_update_gaugefield...";
	md_update_gaugefield = clCreateKernel(clprogram, "md_update_gaugefield", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.debug() << "Create kernel md_update_gaugemomenta...";
	md_update_gaugefield = clCreateKernel(clprogram, "md_update_gaugemomenta", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.debug() << "Create kernel gauge_force...";
	gauge_force = clCreateKernel(clprogram, "gauge_force", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.debug() << "Create kernel fermion_force...";
	fermion_force = clCreateKernel(clprogram, "fermion_force", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	s_gauge = clCreateKernel(clprogram, "s_gauge", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	/** @todo this can most likely be deleted */
	s_fermion = clCreateKernel(clprogram, "s_fermion", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}


	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::finalize_hmc(){

	logger.debug() << "release HMC-variables.." ;
	if(clReleaseMemObject(clmem_energy_init)!=CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_deltah)!=CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_p)!=CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_new_p)!=CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_new_u)!=CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_phi_inv)!=CL_SUCCESS) return HMC_RELEASEVARIABLEERR;

	logger.debug() << "release HMC-kernels.." ;
	if(clReleaseKernel(generate_gaussian_spinorfield)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(s_gauge)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(s_fermion)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	/** @todo these kernels somehow cant be released, perhaps just because they are empty.. */
	if(clReleaseKernel(generate_gaussian_gaugemomenta)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(md_update_gaugefield)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(md_update_gaugemomenta)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(gauge_force)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(fermion_force)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;

	return HMC_SUCCESS;
}

////////////////////////////////////////////////////
//Methods needed for the HMC-algorithm

hmc_error Opencl_hmc::generate_gaussian_gaugemomenta_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::generate_gaussian_spinorfield_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	
	
	(*timer).add();
	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::md_update_spinorfield_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	
	
	(*timer).add();
// #ifdef _FERMIONS_
// //phi = Q+ chi
// hmc_error md_update_spinorfield(hmc_spinor_field * in, hmc_spinor_field * out, hmc_gaugefield * field, inputparameters * parameters){
// 	Qplus(parameters, in, field, out);
// 	return HMC_SUCCESS;
// }
// #endif

	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::leapfrog_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	
	
	(*timer).add();
	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::force_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	
	
	(*timer).add();
	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::hamiltonian_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	//Without fermions here!!! H = S_gauge + S_gaugemomenta
// hmc_float hamiltonian(hmc_gaugefield * field, hmc_float beta, hmc_algebraelement2 * p){
// 	hmc_float result;
// 	(result) = 0.;
// 	(result) += s_gauge(field, beta);
// 	//s_gm = 1/2*squarenorm(Pl)
// 	hmc_float s_gm;
// 	gaugemomenta_squarenorm(p, &s_gm);
// 	result += 0.5*s_gm;
// 	
// 	return result;
// }
	(*timer).reset();
	
	
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::calc_spinorfield_init_energy_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	
	
	(*timer).add();
	return HMC_SUCCESS;
}


////////////////////////////////////////////////////
//Methods to copy new and old fields... these can be optimized!!
hmc_error Opencl_hmc::copy_gaugefield_old_new_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	
	
	(*timer).add();
	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::copy_gaugemomenta_old_new_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::copy_gaugefield_new_old_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	
	
	(*timer).add();
	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::copy_gaugemomenta_new_old_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::get_deltah_from_device(hmc_float * out, const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	hmc_float tmp;
	copy_float_from_device(clmem_deltah, &tmp, timer);
	(*out) = tmp;
	return HMC_SUCCESS;	
}
