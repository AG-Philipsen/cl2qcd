#include "opencl_module_ran.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

void Opencl_Module_Ran::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module::fill_collect_options(collect_options);
	if(get_parameters()->get_use_same_rnd_numbers() ) *collect_options <<  " -D_SAME_RND_NUMBERS_ ";
#ifdef USE_PRNG_NR3
	*collect_options << " -DUSE_PRNG_NR3";
#elif defined(USE_PRNG_RANLUX)
	*collect_options << " -DUSE_PRNG_RANLUX -DRANLUXCL_MAXWORKITEMS=" << num_rndstates;
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}


void Opencl_Module_Ran::fill_buffers()
{

	Opencl_Module::fill_buffers();

#ifdef USE_PRNG_NR3
	// Prepare random number arrays, for each task and device separately
	if(get_device_type() == CL_DEVICE_TYPE_GPU)
		num_rndstates = 5120;
	else
		num_rndstates = 64;
	rndarray = new nr3_state_dev[num_rndstates];
	sizeof_rndarray = sizeof(nr3_state_dev) * num_rndstates;
	nr3_init_seeds(rndarray, "rand_seeds", num_rndstates);

	logger.trace() << "Create buffer for random numbers...";
	clmem_rndarray = create_rw_buffer(sizeof(nr3_state_dev) * get_num_rndstates());
	this->copy_rndstate_to_device(rndarray);
#elif defined(USE_PRNG_RANLUX)
	// make num of random states equal to default num of global threads
	// TODO make this somewhat more automatic (avoid code duplication)
	if(this->get_device_type() == CL_DEVICE_TYPE_GPU)
		num_rndstates = 4 * Opencl_Module::get_numthreads() * get_max_compute_units();
	else
		num_rndstates = get_max_compute_units();

	logger.trace() << "Create buffer for random numbers...";
	clmem_rndarray = create_rw_buffer(7 * num_rndstates * sizeof(cl_float4));
	// kernels are not filled yet, so delay filling until kernel creation
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}

void Opencl_Module_Ran::fill_kernels()
{
	Opencl_Module::fill_kernels();

#ifdef USE_PRNG_NR3
	prng_code = ClSourcePackage() << "random.cl";
#elif defined(USE_PRNG_RANLUX)
	prng_code = ClSourcePackage() << "ranluxcl/ranluxcl.cl" << "random.cl";
	cl_kernel init_kernel = createKernel("prng_ranlux_init") << basic_opencl_code << prng_code << "random_ranlux_init.cl";
	cl_int clerr;
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(init_kernel, this->get_device_type(), &ls, &gs, &num_groups);
	cl_uint seed = get_parameters()->get_host_seed() + 1 + device_rank; // +1 ensures that seed is not equal even if host and device seed are both 0
	if(seed > (10e9 / gs)) { // see ranluxcl source as to why
		/// @todo upgrade to newer ranluxcl to avoid this restcition
		throw Invalid_Parameters("Host seed is too large!", "<< 10e9", get_parameters()->get_host_seed());
	}
	clerr = clSetKernelArg(init_kernel, 0, sizeof(cl_uint), &seed);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(init_kernel, 1, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	enqueueKernel(init_kernel, gs, ls);
	clerr = clFinish(get_queue());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clFinish", __FILE__, __LINE__);
	clerr = clReleaseKernel(init_kernel);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
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
#elif defined(USE_PRNG_RANLUX)
	// nothing to do
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX

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
