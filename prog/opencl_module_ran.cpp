#include "opencl_module_ran.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

void Opencl_Module_Ran::init_random_arrays()
{
	// Prepare random number arrays, for each task and device separately
	if(get_device_type() == CL_DEVICE_TYPE_GPU)
		num_rndstates = 5120;
	else
		num_rndstates = 64;
	rndarray = new hmc_ocl_ran [num_rndstates];
	sizeof_rndarray = sizeof(hmc_ocl_ran) * num_rndstates;
	init_random_seeds(rndarray, "rand_seeds", num_rndstates);
	return;
}


void Opencl_Module_Ran::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module::fill_collect_options(collect_options);
	return;
}


void Opencl_Module_Ran::fill_buffers()
{

	Opencl_Module::fill_buffers();
	init_random_arrays();

	logger.trace() << "Create buffer for random numbers...";
	clmem_rndarray = create_rw_buffer(sizeof(hmc_ocl_ran) * get_num_rndstates());
	this->copy_rndarray_to_device(rndarray);
	return;
}

void Opencl_Module_Ran::fill_kernels()
{
	Opencl_Module::fill_kernels();
	return;
}

void Opencl_Module_Ran::clear_kernels()
{
	Opencl_Module::clear_kernels();
	return;
}

void Opencl_Module_Ran::clear_buffers()
{
	Opencl_Module::clear_buffers();

	delete [] rndarray;

	cl_int clerr = CL_SUCCESS;

	clerr = clReleaseMemObject(clmem_rndarray);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);

	return;
}


void Opencl_Module_Ran::copy_rndarray_to_device(hmc_ocl_ran* rndarray)
{
	cl_int clerr = clEnqueueWriteBuffer(get_queue(), clmem_rndarray, CL_TRUE, 0, sizeof(hmc_ocl_ran) * get_num_rndstates(), rndarray, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);
	return;
}

void Opencl_Module_Ran::copy_rndarray_from_device(hmc_ocl_ran* rndarray)
{
	cl_int clerr = clEnqueueReadBuffer(get_queue(), clmem_rndarray, CL_TRUE, 0, sizeof(hmc_ocl_ran) * get_num_rndstates(), rndarray, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueReadBuffer", __FILE__, __LINE__);

	return;
}

int Opencl_Module_Ran::get_num_rndstates()
{
	return num_rndstates;
}

cl_mem* Opencl_Module_Ran::get_clmem_rndarray()
{
	return &clmem_rndarray;
}

void Opencl_Module_Ran::get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups){
  Opencl_Module::get_work_sizes(kernel, dev_type, ls, gs, num_groups);
  return;
}
