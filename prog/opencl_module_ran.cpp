#include "opencl_module_ran.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

void Opencl_Module_Ran::init_random_arrays()
{
#ifdef USE_PRNG_NR3
	// Prepare random number arrays, for each task and device separately
	if(get_device_type() == CL_DEVICE_TYPE_GPU)
		num_rndstates = 5120;
	else
		num_rndstates = 64;
	rndarray = new nr3_state_dev[num_rndstates];
	sizeof_rndarray = sizeof(nr3_state_dev) * num_rndstates;
	nr3_init_seeds(rndarray, "rand_seeds", num_rndstates);
#else // USE_PRNG_NR3
#error No implemented PRNG selected
#endif // USE_PRNG_NR3
}


void Opencl_Module_Ran::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module::fill_collect_options(collect_options);
	if(get_parameters()->get_use_same_rnd_numbers() ) *collect_options <<  " -D_SAME_RND_NUMBERS_ ";
#ifdef USE_PRNG_NR3
	*collect_options << " -DUSE_PRNG_NR3";
#else // USE_PRNG_NR3
#error No implemented PRNG selected
#endif // USE_PRNG_NR3
}


void Opencl_Module_Ran::fill_buffers()
{

	Opencl_Module::fill_buffers();
	init_random_arrays();

#ifdef USE_PRNG_NR3
	logger.trace() << "Create buffer for random numbers...";
	clmem_rndarray = create_rw_buffer(sizeof(nr3_state_dev) * get_num_rndstates());
	this->copy_rndstate_to_device(rndarray);
#else // USE_PRNG_NR3
#error No implemented PRNG selected
#endif // USE_PRNG_NR3
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

#ifdef USE_PRNG_NR3
	delete [] rndarray;
#else // USE_PRNG_NR3
#error No implemented PRNG selected
#endif // USE_PRNG_NR3

	cl_int clerr = CL_SUCCESS;

	clerr = clReleaseMemObject(clmem_rndarray);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);

	return;
}

#ifdef USE_PRNG_NR3
void Opencl_Module_Ran::copy_rndstate_to_device(nr3_state_dev* rndarray)
{
	cl_int clerr = clEnqueueWriteBuffer(get_queue(), clmem_rndarray, CL_TRUE, 0, sizeof(nr3_state_dev) * get_num_rndstates(), rndarray, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);
	return;
}

void Opencl_Module_Ran::copy_rndstate_from_device(nr3_state_dev* rndarray)
{
	cl_int clerr = clEnqueueReadBuffer(get_queue(), clmem_rndarray, CL_TRUE, 0, sizeof(nr3_state_dev) * get_num_rndstates(), rndarray, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueReadBuffer", __FILE__, __LINE__);

	return;
}
#endif // USE_PRNG_NR3

int Opencl_Module_Ran::get_num_rndstates()
{
	return num_rndstates;
}

cl_mem* Opencl_Module_Ran::get_clmem_rndarray()
{
	return &clmem_rndarray;
}

void Opencl_Module_Ran::get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups)
{
	Opencl_Module::get_work_sizes(kernel, dev_type, ls, gs, num_groups);
	return;
}
