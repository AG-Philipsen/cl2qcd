#include "opencl_module_ran.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

void Opencl_Module_Ran::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module::fill_collect_options(collect_options);
	if(get_parameters().get_use_same_rnd_numbers() ) *collect_options <<  " -D_SAME_RND_NUMBERS_ ";
#ifdef USE_PRNG_NR3
	*collect_options << " -DUSE_PRNG_NR3";
#elif defined(USE_PRNG_RANLUX)
	*collect_options << " -DUSE_PRNG_RANLUX -DRANLUXCL_MAXWORKITEMS=" << prng_buffer.get_elements();
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}


void Opencl_Module_Ran::fill_buffers()
{

	Opencl_Module::fill_buffers();

#ifdef USE_PRNG_NR3
	// Prepare random number arrays, for each task and device separately
	const size_t num_rndstates = prng_buffer.get_elements();
	rndarray = new nr3_state_dev[num_rndstates];
	nr3_init_seeds(rndarray, "rand_seeds", num_rndstates);
	prng_buffer.load(rndarray);
#elif defined(USE_PRNG_RANLUX)
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
	this->get_work_sizes(init_kernel, &ls, &gs, &num_groups);
	// FIXME reenable device rank
	cl_uint seed = get_parameters().get_host_seed() + 1;// + device_rank; // +1 ensures that seed is not equal even if host and device seed are both 0
	if(seed > (10e9 / gs)) { // see ranluxcl source as to why
		/// @todo upgrade to newer ranluxcl to avoid this restcition
		throw Invalid_Parameters("Host seed is too large!", "<< 10e9", get_parameters().get_host_seed());
	}
	clerr = clSetKernelArg(init_kernel, 0, sizeof(cl_uint), &seed);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(init_kernel, 1, sizeof(cl_mem), prng_buffer);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	get_device()->enqueue_kernel(init_kernel, gs, ls);
	clerr = clFinish(get_queue());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clFinish", __FILE__, __LINE__);
	clerr = clReleaseKernel(init_kernel);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
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
}

#ifdef USE_PRNG_NR3
void Opencl_Module_Ran::copy_rndstate_to_device(nr3_state_dev* rndarray) const
{
	prng_buffer.load(rndarray);
}

void Opencl_Module_Ran::copy_rndstate_from_device(nr3_state_dev* rndarray) const
{
	prng_buffer.dump(rndarray);
}
#endif // USE_PRNG_NR3

const hardware::buffers::PRNGBuffer& Opencl_Module_Ran::get_prng_buffer() const noexcept
{
	return prng_buffer;
}
