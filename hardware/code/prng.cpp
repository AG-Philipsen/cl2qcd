/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "prng.hpp"

using namespace std;

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params)
{
	std::ostringstream options;
	options.precision(16);
	if(params.get_use_same_rnd_numbers() ) options <<  " -D _SAME_RND_NUMBERS_ ";
#ifdef USE_PRNG_RANLUX
	options << "-D USE_PRNG_RANLUX -D RANLUXCL_MAXWORKITEMS=" << hardware::buffers::get_prng_buffer_size(device, params.get_use_same_rnd_numbers());
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
	return options.str();
}

hardware::code::PRNG::PRNG(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module(params, device)
{
#ifdef USE_PRNG_RANLUX
	logger.debug() << "Creating PRNG kernels...";
	// the ranluxcl lies in the main directory, other than the remaining kernel
	prng_code = ClSourcePackage(collect_build_options(get_device(), get_parameters())) << "../ranluxcl/ranluxcl.cl" << "random.cl";
	init_kernel = createKernel("prng_ranlux_init") << ClSourcePackage("-I " + std::string(SOURCEDIR) + " -D _INKERNEL_") << "globaldefs.h" << "types.h" << "opencl_header.cl" <<  prng_code << "random_ranlux_init.cl";
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}

hardware::code::PRNG::~PRNG()
{
#ifdef USE_PRNG_RANLUX
	logger.debug() << "Clearing PRNG kernels...";
	cl_int clerr = clReleaseKernel(init_kernel);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}

ClSourcePackage hardware::code::PRNG::get_sources() const noexcept
{
	return prng_code;
}

#ifdef USE_PRNG_RANLUX
void hardware::code::PRNG::initialize(const hardware::buffers::PRNGBuffer * buffer, cl_uint seed) const
{
	cl_int clerr;
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(init_kernel, &ls, &gs, &num_groups);
	// we need a custom global size
	gs = buffer->get_elements();
	if(seed > (10e9 / gs)) { // see ranluxcl source as to why
		/// @todo upgrade to newer ranluxcl to avoid this restcition
		throw Invalid_Parameters("Host seed is too large!", "<< 10e9", (int)get_parameters().get_host_seed());
	}
	clerr = clSetKernelArg(init_kernel, 0, sizeof(cl_uint), &seed);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(init_kernel, 1, sizeof(cl_mem), buffer->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	get_device()->enqueue_kernel(init_kernel, gs, ls);
}
#endif /* USE_PRNG_RANLUX */
