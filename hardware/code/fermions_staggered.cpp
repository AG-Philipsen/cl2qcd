/** @file
 * Implementation of the hardware::code::fermions_staggered class
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
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

#include "fermions_staggered.hpp"

#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include "spinors_staggered.hpp"
#include "spinors.hpp"

#include <cassert>
#include <cmath>
#include "flopUtilities.hpp"

using namespace std;

void hardware::code::Fermions_staggered::fill_kernels()
{
	sources = get_basic_sources() <<  "operations_geometry.cl" << "operations_complex.h"  << "operations_matrix_su3.cl" << "operations_matrix.cl" << "operations_gaugefield.cl" << "types_fermions.h" << "operations_su3vec.cl" << "operations_staggered.cl";
	if(kernelParameters->getUseEo()) {
		sources = sources << "spinorfield_staggered_eo.cl";;
	} else {
		sources = sources << "spinorfield_staggered.cl";
	}

	logger.debug() << "Creating Fermions_staggered kernels...";

	if(kernelParameters->getFermact() == common::action::rooted_stagg) {
	      if(kernelParameters->getUseEo()){
			M_staggered = 0;
			D_KS_eo = createKernel("D_KS_eo") << sources << "fermionmatrix_staggered_eo_DKS_local.cl" << "fermionmatrix_staggered_eo_DKS.cl";
	      } else {
			D_KS_eo = 0;
			M_staggered = createKernel("M_staggered") << sources << "fermionmatrix_staggered_DKS_local.cl" << "fermionmatrix_staggered_M.cl";
	      }
	} else {
		throw Print_Error_Message("Fermions_staggered module asked to be built but action set not to rooted_stagg! Aborting... ", __FILE__, __LINE__);
	}
}

void hardware::code::Fermions_staggered::clear_kernels()
{
	cl_uint clerr = CL_SUCCESS;
	
	logger.debug() << "Clearing Fermions_staggered kernels...";

	if(kernelParameters->getUseEo()){
		clerr = clReleaseKernel(D_KS_eo);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	} else {
		clerr = clReleaseKernel(M_staggered);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
}

void hardware::code::Fermions_staggered::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
    Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);
    if(kernel == D_KS_eo){
        if(*ls > 64) {
            *ls = 64;
            *num_groups = (*gs)/(*ls);
        }
        return;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//explicit fermionmatrix-kernel calling functions
void hardware::code::Fermions_staggered::M_staggered_device(const hardware::buffers::Plain<su3vec> * in, const hardware::buffers::Plain<su3vec> * out, const hardware::buffers::SU3 * gf, hmc_float mass) const
{
	//get mass
	hmc_float mass_tmp;
	if(mass == ARG_DEF) mass_tmp = kernelParameters->getMass();
	else mass_tmp = mass;

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_staggered, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_staggered, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_staggered, 1, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_staggered, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_staggered, 3, sizeof(hmc_float), &mass_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( M_staggered, gs2, ls2);
}


void hardware::code::Fermions_staggered::D_KS_eo_device(const hardware::buffers::SU3vec * in, const hardware::buffers::SU3vec * out, const hardware::buffers::SU3 * gf, int evenodd) const
{
	cl_int eo = evenodd;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(D_KS_eo, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(D_KS_eo, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(D_KS_eo, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(D_KS_eo, 2, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(D_KS_eo, 3, sizeof(cl_int), &eo);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(D_KS_eo, gs2, ls2);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////



size_t hardware::code::Fermions_staggered::get_read_write_size(const std::string& in) const
{
	//Depending on the compile-options, one has different sizes...
	size_t D = kernelParameters->getFloatSize();
	//this returns the number of entries in an su3-matrix
	size_t R = kernelParameters->getMatSize();
	//this is the number of su3vec in the system (or number of sites)
	size_t S = kernelParameters->getSpinorFieldSize();
	size_t Seo = kernelParameters->getEoprecSpinorFieldSize();
	//factor for complex numbers
	int C = 2;
	//NOTE: 1 spinor has NC = 3 complex entries
	if (in == "M_staggered") {
		//this kernel reads 9 su3vec, 8 su3matrices and writes 1 su3vec per site:
		const unsigned int dirs = 4;
		return (C * NC * (2 * dirs + 1 + 1) + C * 2 * dirs * R) * D * S; //1632 bytes * S
	}
	if (in == "D_KS_eo") {
		//this kernel reads 8 su3vec (not that in the site of the output),
		//8 su3matrices and writes 1 su3vec:
		const unsigned int dirs = 4;
		return (C * NC * (2 * dirs + 1) + C * 2 * dirs * R) * D * Seo; //1584 bytes * Seo
	}
	return 0;
}

/**
 * This function returns the number of flops that is needed to make the standard staggered
 * Dirac operator act onto a field (582 flops). Since in general the mass in the simulation is not
 * zero (and there are no check on that), the mass-term in the d_slash is always taken
 * into account (for standard Dirac operator we intend M = D_KS + m). Working with evenodd
 * preconditioning, then the mass term is not put in the kernel and this is taken into account
 * in the number of flops calculation.
 * 
 * @attention In this function (as in the whole code) the staggered phases are not included
 *            in links and, then, we have some flops in addition to take them into account
 *            as well as for the boundary conditions. Usually, in the community, in all
 *            codes with staggered fermions, staggered phases are calculated only once and
 *            are put in links (including also boundary conditions). In this way staggered phases
 *            and BC are completely forgotten. Nevertheless, here we prefer to calculate them
 *            every time that we need them, not putting them into links. There are basically
 *            two reasons to do that. \n
 *            \arg Flops needed to calculate staggered phases and BC do not lower the
 *                 performance, since they are few flops and "Flops don't count".\n
 *            \arg It is really important not to change gauge part of the code, so that
 *                 both staggered and Wilson codes share the same gauge tools.\n
 *            Hence, our code will perform more floating point operations in the kernels than
 *            codes that include phases in links. But this does not mean that we must count
 *            these additional operations in the number of total flop, because otherwise
 *            our performance is not comparable to other codes. To understand why, one can give
 *            the following argument. The best benchmark would be in term of needed time to execute
 *            the kernel: if one finds a better way to implement something, the time needed is
 *            smaller and the performance higher. But if changing the code, one adapt also the
 *            number of floating point operations, then the measured performance does not reflect
 *            the truth. It can happen that the measured performance does not improve when the
 *            time to execute the kernel has lowered. Another example is the following: suppose
 *            you do a lot of unnecessary operations in the kernel and take them into account.
 *            Since LQCD is really memory limited you will have more or less the same time with
 *            and without those operations. But if you count them, then your benchmark will look
 *            amazing, when the time needed to execute the kernel is probably the same as other
 *            codes.
 */
/*
 * Considering that:
 * 
 *  - su3matrix times su3vec is 66 flops
 *  - real times complex is 2 flops
 *  - complex times complex is 6 flops
 *  - real times su3vec is 6 flops
 *  - complex times su3vec is 18 flops
 *  - su3vec plus su3vec is 6 flops = NC*2
 * 
 * The way of counting flops for the dslash can be summarized as follows:
 * 
 *  + For each direction we have twice su3matrix times su3vec
 *  + We have to sum the 2 resulting su3vec in each direction
 *  + The mass term is a real times su3vec
 *  + Then we have to sum the 5 su3vec (one for each direction + the mass term)
 * 
 */
static int flop_dslash_staggered_per_site()
{
	return NDIM * (2 * hardware::code::getFlopSu3MatrixTimesSu3Vec() + NC * 2) + NC * 2 + 4 * NC * 2;
}

/**
 * This function is the same as flop_dslash_staggered_per_site, except the fact that here
 * the mass term is not taken into account (return 570 flops)
 */
static int flop_dks_staggered_per_site()
{
	return NDIM * (2 * hardware::code::getFlopSu3MatrixTimesSu3Vec() + NC * 2) + 3 * NC * 2;
}

uint64_t hardware::code::Fermions_staggered::get_flop_size(const std::string& in) const
{
	size_t S = kernelParameters->getSpinorFieldSize();
	size_t Seo = kernelParameters->getEoprecSpinorFieldSize();
	if (in == "M_staggered") {
		//this kernel performs one dslash on each site. To do that it also calculates
		//the staggered phases and the BC but we do not count this flop (see explanation above)
		return S * flop_dslash_staggered_per_site(); // S * 582 flop
	}
	if (in == "D_KS_eo") {
		//this kernel performs one dks on each site. To do that it also calculates
		//the staggered phases and the BC but we do not count this flop (see explanation above)
		return Seo * flop_dks_staggered_per_site(); // Seo * 570 flop
	}
	return 0;
}

void hardware::code::Fermions_staggered::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	if(M_staggered)
		Opencl_Module::print_profiling(filename, M_staggered);
	if(D_KS_eo)
		Opencl_Module::print_profiling(filename, D_KS_eo);
}

hardware::code::Fermions_staggered::Fermions_staggered(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device)
	: Opencl_Module(kernelParameters, device)
{
	fill_kernels();
}

hardware::code::Fermions_staggered::~Fermions_staggered()
{
	clear_kernels();
}

ClSourcePackage hardware::code::Fermions_staggered::get_sources() const noexcept
{
	return sources;
}
