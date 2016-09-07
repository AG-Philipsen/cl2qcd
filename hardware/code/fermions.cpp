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

#include "fermions.hpp"

#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include "spinors.hpp"

#include <cassert>
#include <cmath>
#include "flopUtilities.hpp"

using namespace std;

void hardware::code::Fermions::fill_kernels()
{
	sources = get_basic_sources() << "operations_geometry.cl" << "operations_complex.h" << "types_fermions.h" << "operations_matrix_su3.cl" << "operations_matrix.cl" << "operations_gaugefield.cl" << "operations_su3vec.cl" << "operations_spinor.cl" << "spinorfield.cl";
	if(kernelParameters->getUseEo()) {
		sources = sources << "operations_spinorfield_eo.cl";
	}

	logger.debug() << "Creating Fermions kernels...";
	
	if(kernelParameters->getFermact() == common::action::wilson) {
		M_wilson = createKernel("M_wilson") << sources << "fermionmatrix.cl" << "fermionmatrix_m.cl";
	} else if(kernelParameters->getFermact() == common::action::twistedmass) {
		M_tm_plus = createKernel("M_tm_plus") << sources << "fermionmatrix.cl" << "fermionmatrix_m_tm_plus.cl";
		M_tm_minus = createKernel("M_tm_minus") << sources << "fermionmatrix.cl" << "fermionmatrix_m_tm_minus.cl";
	} else if(kernelParameters->getFermact() == common::action::clover) {
		throw Print_Error_Message("no kernels for CLOVER-discretization implemented yet, aborting... ", __FILE__, __LINE__);
	} else {
		throw Print_Error_Message("there was a problem with which fermion-discretization to use, aborting... ", __FILE__, __LINE__);
	}
	
	gamma5 = createKernel("gamma5") << sources << "fermionmatrix.cl" << "fermionmatrix_gamma5.cl";
	//Kernels needed if eoprec is used
	if(kernelParameters->getUseEo() == true) {
		if(kernelParameters->getFermact() == common::action::twistedmass) {
			M_tm_sitediagonal = createKernel("M_tm_sitediagonal") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m.cl";
			M_tm_inverse_sitediagonal = createKernel("M_tm_inverse_sitediagonal") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m.cl";
			M_tm_sitediagonal_minus = createKernel("M_tm_sitediagonal_minus") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m.cl";
			M_tm_inverse_sitediagonal_minus = createKernel("M_tm_inverse_sitediagonal_minus") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m.cl";
		}
		dslash_eo = createKernel("dslash_eo") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_dslash.cl";
		_dslash_eo_boundary = createKernel("dslash_eo_boundary") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_dslash.cl";
		_dslash_eo_inner = createKernel("dslash_eo_inner") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_dslash.cl";
		gamma5_eo = createKernel("gamma5_eo") << sources << "fermionmatrix.cl" << "fermionmatrix_eo_gamma5.cl";
		//merged kernels
		if (kernelParameters->getUseMergeKernelsFermion() == true) {
			dslash_AND_M_tm_inverse_sitediagonal_eo = createKernel("dslash_AND_M_tm_inverse_sitediagonal_eo") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_dslash_AND_M_tm_inverse_sitediagonal.cl";
			dslash_AND_M_tm_inverse_sitediagonal_minus_eo = createKernel("dslash_AND_M_tm_inverse_sitediagonal_minus_eo") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_dslash_AND_M_tm_inverse_sitediagonal_minus.cl";
			M_tm_sitediagonal_AND_gamma5_eo = createKernel("M_tm_sitediagonal_AND_gamma5_eo") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m_merged.cl";
			M_tm_sitediagonal_minus_AND_gamma5_eo = createKernel("M_tm_sitediagonal_minus_AND_gamma5_eo") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m_merged.cl";
			saxpy_AND_gamma5_eo = createKernel("saxpy_AND_gamma5_eo") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_saxpy_AND_gamma5_eo.cl";
		}
	}
}

void hardware::code::Fermions::clear_kernels()
{
	cl_uint clerr = CL_SUCCESS;
	
	logger.debug() << "Clearing Fermions kernels...";

	if(M_wilson) {
		clerr = clReleaseKernel(M_wilson);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(M_tm_plus) {
		clerr = clReleaseKernel(M_tm_plus);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(M_tm_minus) {
		clerr = clReleaseKernel(M_tm_minus);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(gamma5) {
		clerr = clReleaseKernel(gamma5);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}

	if(kernelParameters->getUseEo()) {
		if(dslash_eo)
		{
			clerr = clReleaseKernel(dslash_eo);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		}
		if(_dslash_eo_boundary)
		{
			clerr = clReleaseKernel(_dslash_eo_boundary);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		}
		if(_dslash_eo_inner)
		{
			clerr = clReleaseKernel(_dslash_eo_inner);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		}
		if (gamma5_eo)
		{
			clerr = clReleaseKernel(gamma5_eo);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		}
		if(kernelParameters->getFermact() == common::action::twistedmass) {
			if(M_tm_sitediagonal)
			{
				clerr = clReleaseKernel(M_tm_sitediagonal);
				if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
			}
			if (M_tm_inverse_sitediagonal)
			{
				clerr = clReleaseKernel(M_tm_inverse_sitediagonal);
				if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
			}
			if(M_tm_sitediagonal_minus)
			{
				clerr = clReleaseKernel(M_tm_sitediagonal_minus);
				if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
			}
			if(M_tm_inverse_sitediagonal_minus)
			{
				clerr = clReleaseKernel(M_tm_inverse_sitediagonal_minus);
				if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
			}
		}
		if (kernelParameters->getUseMergeKernelsFermion() == true)
		{
			if(saxpy_AND_gamma5_eo)
			{
				clerr = clReleaseKernel(saxpy_AND_gamma5_eo);
				if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
			}
			if(dslash_AND_M_tm_inverse_sitediagonal_eo)
			{
				clerr = clReleaseKernel(dslash_AND_M_tm_inverse_sitediagonal_eo);
				if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
			}
			if(dslash_AND_M_tm_inverse_sitediagonal_minus_eo)
			{
				clerr = clReleaseKernel(dslash_AND_M_tm_inverse_sitediagonal_minus_eo);
				if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
			}
			if(M_tm_sitediagonal_AND_gamma5_eo)
			{
				clerr = clReleaseKernel(M_tm_sitediagonal_AND_gamma5_eo);
				if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
			}
			if(M_tm_sitediagonal_minus_AND_gamma5_eo)
			{
				clerr = clReleaseKernel(M_tm_sitediagonal_minus_AND_gamma5_eo);
				if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
			}
		}
	}
}

void hardware::code::Fermions::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);
	
	//Query specific sizes for kernels if needed
	if(kernel == saxpy_AND_gamma5_eo){
		if(*ls > 64) {
			*ls = 64;
			*num_groups = (*gs)/(*ls);
		}
		return;
	}
}


//explicit fermionmatrix-kernel calling functions
void hardware::code::Fermions::M_wilson_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa) const
{
	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = kernelParameters->getKappa();
	else kappa_tmp = kappa;

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_wilson, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_wilson, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_wilson, 1, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_wilson, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_wilson, 3, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( M_wilson, gs2, ls2);
}

void hardware::code::Fermions::M_tm_plus_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa , hmc_float mubar ) const
{
	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = kernelParameters->getKappa();
	else kappa_tmp = kappa;

	//get mu
	hmc_float mubar_tmp;
	if(mubar == ARG_DEF) mubar_tmp = kernelParameters->getMuBar();
	else mubar_tmp = mubar;

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_tm_plus, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_plus, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_plus, 1, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_plus, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_plus, 3, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_plus, 4, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( M_tm_plus, gs2, ls2);
}

void hardware::code::Fermions::M_tm_minus_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa , hmc_float mubar ) const
{
	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = kernelParameters->getKappa();
	else kappa_tmp = kappa;

	//get mu
	hmc_float mubar_tmp;
	if(mubar == ARG_DEF) mubar_tmp = kernelParameters->getMuBar();
	else mubar_tmp = mubar;

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_tm_minus, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_minus, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_minus, 1, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_minus, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_minus, 3, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_minus, 4, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( M_tm_minus, gs2, ls2);
}

void hardware::code::Fermions::gamma5_device(const hardware::buffers::Plain<spinor> * inout) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gamma5, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gamma5, 0, sizeof(cl_mem), inout->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(gamma5 , gs2, ls2);
}

//explicit eoprec fermionmatrix functions
void hardware::code::Fermions::gamma5_eo_device(const hardware::buffers::Spinor * inout) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gamma5_eo, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gamma5_eo, 0, sizeof(cl_mem), inout->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( gamma5_eo, gs2, ls2);
}

void hardware::code::Fermions::saxpy_AND_gamma5_eo_device(const hardware::buffers::Spinor * x, const hardware::buffers::Spinor * y, const hmc_complex alpha, const hardware::buffers::Spinor * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy_AND_gamma5_eo, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy_AND_gamma5_eo, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_AND_gamma5_eo, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_AND_gamma5_eo, 2, sizeof(hmc_float), &alpha.re);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_AND_gamma5_eo, 3, sizeof(hmc_float), &alpha.im);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_AND_gamma5_eo, 4, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( saxpy_AND_gamma5_eo, gs2, ls2);
}

void hardware::code::Fermions::dslash_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa) const
{
	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = kernelParameters->getKappa();
	else kappa_tmp = kappa;

	cl_int eo = evenodd;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(dslash_eo, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(dslash_eo, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_eo, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_eo, 2, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_eo, 3, sizeof(cl_int), &eo);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_eo, 4, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(dslash_eo , gs2, ls2);
}

void hardware::code::Fermions::dslash_eo_boundary(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa) const
{
	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = kernelParameters->getKappa();
	else kappa_tmp = kappa;

	cl_int eo = evenodd;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(_dslash_eo_boundary, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(_dslash_eo_boundary, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(_dslash_eo_boundary, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(_dslash_eo_boundary, 2, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(_dslash_eo_boundary, 3, sizeof(cl_int), &eo);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(_dslash_eo_boundary, 4, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(_dslash_eo_boundary , gs2, ls2);
}

void hardware::code::Fermions::dslash_eo_inner(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa) const
{
	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = kernelParameters->getKappa();
	else kappa_tmp = kappa;

	cl_int eo = evenodd;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(_dslash_eo_inner, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(_dslash_eo_inner, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(_dslash_eo_inner, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(_dslash_eo_inner, 2, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(_dslash_eo_inner, 3, sizeof(cl_int), &eo);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(_dslash_eo_inner, 4, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(_dslash_eo_inner , gs2, ls2);
}

void hardware::code::Fermions::dslash_AND_M_tm_inverse_sitediagonal_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa, hmc_float mubar) const
{
	cl_int eo = evenodd;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(dslash_AND_M_tm_inverse_sitediagonal_eo, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_eo, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_eo, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_eo, 2, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_eo, 3, sizeof(cl_int), &eo);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_eo, 4, sizeof(hmc_float), &kappa);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_eo, 5, sizeof(hmc_float), &mubar);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(dslash_AND_M_tm_inverse_sitediagonal_eo , gs2, ls2);
}

void hardware::code::Fermions::dslash_AND_M_tm_inverse_sitediagonal_minus_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa, hmc_float mubar) const
{
	cl_int eo = evenodd;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, 2, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, 3, sizeof(cl_int), &eo);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, 4, sizeof(hmc_float), &kappa);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, 5, sizeof(hmc_float), &mubar);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(dslash_AND_M_tm_inverse_sitediagonal_minus_eo , gs2, ls2);
}

void hardware::code::Fermions::M_tm_inverse_sitediagonal_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar) const
{
	//get mu
	hmc_float mubar_tmp;
	if(mubar == ARG_DEF) mubar_tmp = kernelParameters->getMuBar();
	else mubar_tmp = mubar;

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_tm_inverse_sitediagonal, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_inverse_sitediagonal, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_inverse_sitediagonal, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_inverse_sitediagonal, 2, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( M_tm_inverse_sitediagonal, gs2, ls2);
}

void hardware::code::Fermions::M_tm_sitediagonal_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar) const
{
	//get mu
	hmc_float mubar_tmp;
	if(mubar == ARG_DEF) mubar_tmp = kernelParameters->getMuBar();
	else mubar_tmp = mubar;

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_tm_sitediagonal, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_sitediagonal, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal, 2, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(M_tm_sitediagonal , gs2, ls2);
}

void hardware::code::Fermions::M_tm_sitediagonal_AND_gamma5_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar) const
{
	//get mu
	hmc_float mubar_tmp;
	if(mubar == ARG_DEF) mubar_tmp = kernelParameters->getMuBar();
	else mubar_tmp = mubar;

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_tm_sitediagonal_AND_gamma5_eo, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_sitediagonal_AND_gamma5_eo, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal_AND_gamma5_eo, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal_AND_gamma5_eo, 2, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(M_tm_sitediagonal_AND_gamma5_eo , gs2, ls2);
}

void hardware::code::Fermions::M_tm_sitediagonal_minus_AND_gamma5_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar) const
{
	//get mu
	hmc_float mubar_tmp;
	if(mubar == ARG_DEF) mubar_tmp = kernelParameters->getMuBar();
	else mubar_tmp = mubar;

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_tm_sitediagonal_minus_AND_gamma5_eo, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_sitediagonal_minus_AND_gamma5_eo, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal_minus_AND_gamma5_eo, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal_minus_AND_gamma5_eo, 2, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(M_tm_sitediagonal_minus_AND_gamma5_eo , gs2, ls2);
}

void hardware::code::Fermions::M_tm_inverse_sitediagonal_minus_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar) const
{
	//get mu
	hmc_float mubar_tmp;
	if(mubar == ARG_DEF) mubar_tmp = kernelParameters->getMuBar();
	else mubar_tmp = mubar;

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_tm_inverse_sitediagonal_minus, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_inverse_sitediagonal_minus, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_inverse_sitediagonal_minus, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_inverse_sitediagonal_minus, 2, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( M_tm_inverse_sitediagonal_minus, gs2, ls2);
}

void hardware::code::Fermions::M_tm_sitediagonal_minus_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar) const
{
	//get mu
	hmc_float mubar_tmp;
	if(mubar == ARG_DEF) mubar_tmp = kernelParameters->getMuBar();
	else mubar_tmp = mubar;

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_tm_sitediagonal_minus, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_sitediagonal_minus, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal_minus, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal_minus, 2, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(M_tm_sitediagonal_minus , gs2, ls2);
}

size_t hardware::code::Fermions::get_read_write_size(const std::string& in) const
{
	//Depending on the compile-options, one has different sizes...
	size_t D = kernelParameters->getFloatSize();
	//this returns the number of entries in an su3-matrix
	size_t R = kernelParameters->getMatSize();
	//this is the number of spinors in the system (or number of sites)
	size_t S = kernelParameters->getSpinorFieldSize();
	size_t Seo = kernelParameters->getEoprecSpinorFieldSize();
	//factor for complex numbers
	int C = 2;
	//this is the same as in the function above
	//NOTE: 1 spinor has NC*NDIM = 12 complex entries
	if (in == "M_wilson") {
		//this kernel reads 9 spinors, 8 su3matrices and writes 1 spinor:
		return (C * 12 * (9 + 1) + C * 8 * R) * D * S;
	}
	if (in == "gamma5") {
		//this kernel reads 1 spinor and writes 1 spinor:
		return 2 * C * 12 * D * S;
	}
	if (in == "M_tm_plus") {
		//this kernel reads 9 spinors, 8 su3matrices and writes 1 spinor:
		return (C * 12 * (9 + 1) + C * 8 * R) * D * S;
	}
	if (in == "M_tm_minus") {
		//this kernel reads 9 spinors, 8 su3matrices and writes 1 spinor:
		return (C * 12 * (9 + 1) + C * 8 * R) * D * S;
	}
	if (in == "gamma5_eo") {
		//this kernel reads 1 spinor and writes 1 spinor:
		return 48 * D * Seo;
	}
	if (in == "M_tm_sitediagonal") {
		//this kernel reads 1 spinor and writes 1 spinor:
		return 48 * D * Seo;
	}
	if (in == "M_tm_inverse_sitediagonal") {
		//this kernel reads 1 spinor and writes 1 spinor:
		return 48 * D * Seo;
	}
	if (in == "M_tm_sitediagonal_minus") {
		//this kernel reads 1 spinor and writes 1 spinor:
		return 48 * D * Seo;
	}
	if (in == "M_tm_inverse_sitediagonal_minus") {
		//this kernel reads 1 spinor and writes 1 spinor:
		return 48 * D * Seo;
	}
	if (in == "dslash_eo") {
		//this kernel reads 8 spinors, 8 su3matrices and writes 1 spinor:
		const unsigned int dirs = 4;
		return (C * 12 * (2 * dirs + 1) + C * 2 * dirs * R) * D * Seo;
	}
	if (in == "dslash_AND_M_tm_inverse_sitediagonal_eo") {
		//the dslash kernel reads 8 spinors, 8 su3matrices and writes 1 spinor:
		const unsigned int dirs = 4;
		//the diag kernel reads 1 spinor and writes 1 spinor:
		//the merged kernel reads 8 spinors, 8 su3matrices and writes 1 spinor, thus it is the same as the dslash
		return  (C * 12 * (2 * dirs + 1) + C * 2 * dirs * R) * D * Seo;
	}
	if (in == "dslash_AND_M_tm_inverse_sitediagonal_minus_eo") {
		//the dslash kernel reads 8 spinors, 8 su3matrices and writes 1 spinor:
		const unsigned int dirs = 4;
		//the diag kernel reads 1 spinor and writes 1 spinor:
		//the merged kernel reads 8 spinors, 8 su3matrices and writes 1 spinor, thus it is the same as the dslash
		return  (C * 12 * (2 * dirs + 1) + C * 2 * dirs * R) * D * Seo;
	}
	if (in == "M_tm_sitediagonal_AND_gamma5_eo") {
		//this kernel reads 1 spinor and writes 1 spinor:
		return 48 * D * Seo;
	}
	if (in == "M_tm_sitediagonal_minus_AND_gamma5_eo") {
		//this kernel reads 1 spinor and writes 1 spinor:
		return 48 * D * Seo;
	}
	if (in == "saxpy_AND_gamma5_eo") {
		//saxpy reads 2 spinor, 1 complex number and writes 1 spinor per site
		//the gamma5 does not affect this here.
		return C * D * Seo * (12 * (2 + 1) + 1);
	}
	return 0;
}

static int flop_dslash_per_site()
{
	/** @NOTE: this is the "original" dslash without any simplifications, counting everything "full". this is a much too hight number!!
	   *  //this kernel performs for each eo site a 2*NDIM sum over (1 + gamma_mu) * su3matrix * spinor
	   *  //return  NDIM * 2 * ( get_parameters().get_flop_su3_spinor() + get_parameters().get_flop_gamma_spinor() ) ;
	  @NOTE: However, in 0911.3191 the dslash_eo is quoted to perform 1320 flop per site
	   *  If I count our implementation of the dslash-kernel, I get 1632 flop per site:
	   *  //the kernel performs 6 su3vec_acc, 2 su3_su3vec and 2 su3vec_complex in NDIM * 2 directions per site
	  */
	return NDIM * 2 * (hardware::code::getFlopSu3MatrixTimesSu3Vec() * 2 + 6 * NC * 2 + 2 * NC * hardware::code::getFlopComplexMult() );
	/// I do not know for sure, but if one leaves out the 2 su3vec_complex operations in each directions one has almost 1320 flop per site. This would correspond to have the kappa in the diagonal matrix still.

}

uint64_t hardware::code::Fermions::get_flop_size(const std::string& in) const
{
	size_t S = kernelParameters->getSpinorFieldSize();
	size_t Seo = kernelParameters->getEoprecSpinorFieldSize();
	if (in == "M_wilson") {
		//this kernel performs one dslash on each site and adds this to a spinor
		return S * (flop_dslash_per_site() + NC * NDIM * getFlopComplexMult() + NC * NDIM * 2 );
	}
	if (in == "gamma5") {
		//this kernel performs ND*NC*2/2 real mults
		return S * NDIM * NC;
	}
	if (in == "M_tm_plus") {
		//this kernel performs ND*NC complex mults and one dslash on each site and adds the results
		return S * (flop_dslash_per_site() + NC * NDIM * getFlopComplexMult() + NC * NDIM * 2 );
	}
	if (in == "M_tm_minus") {
		//this kernel performs ND*NC complex mults and one dslash on each site and adds the results
		return S * (flop_dslash_per_site() + NC * NDIM * getFlopComplexMult() + NC * NDIM * 2 );
	}
	if (in == "gamma5_eo") {
		//this kernel performs ND*NC*2/2 real mults
		return Seo * NDIM * NC;
	}
	if (in == "M_tm_sitediagonal") {
		//this kernel performs ND*NC complex mults
		return Seo * ( NC * NDIM * getFlopComplexMult() );
	}
	if (in == "M_tm_inverse_sitediagonal") {
		//this kernel performs ND*NC complex mults and ND*NC*2 real mults
		return Seo * ( NC * NDIM * getFlopComplexMult() + NC * NDIM * 2  );
	}
	if (in == "M_tm_sitediagonal_minus") {
		//this kernel performs ND*NC complex mults
		return Seo * ( NC * NDIM * getFlopComplexMult() );
	}
	if (in == "M_tm_inverse_sitediagonal_minus") {
		//this kernel performs ND*NC complex mults and ND*NC*2 real mults
		return Seo * ( NC * NDIM * getFlopComplexMult() + NC * NDIM * 2 );
	}
	if (in == "dslash_eo") {
		return Seo * flop_dslash_per_site();
	}
	if (in == "dslash_AND_M_tm_inverse_sitediagonal_eo") {
		return Seo * flop_dslash_per_site() + Seo * ( NC * NDIM * getFlopComplexMult() + NC * NDIM * 2  );
	}
	if (in == "dslash_AND_M_tm_inverse_sitediagonal_minus_eo") {
		return Seo * flop_dslash_per_site() + Seo * ( NC * NDIM * getFlopComplexMult() + NC * NDIM * 2  );
	}
	if (in == "M_tm_sitediagonal_AND_gamma5_eo") {
		//this kernel performs ND*NC complex mults and  ND*NC*2/2 real mults
		return Seo * ( NC * NDIM * getFlopComplexMult() ) + Seo * NDIM * NC ;
	}
	if (in == "M_tm_sitediagonal_minus_AND_gamma5_eo") {
		//this kernel performs ND*NC complex mults  ND*NC*2/2 real mults
		return Seo * ( NC * NDIM * getFlopComplexMult() ) +  Seo * NDIM * NC;
	}
	if (in == "saxpy_AND_gamma5_eo") {
		//saxpy performs on each site spinor_times_complex and spinor_add
		//gamma5 performs ND*NC*2/2 real mults
		return Seo * NDIM * NC * (1+  getFlopComplexMult() + 2);
	}
	return 0;
}

void hardware::code::Fermions::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, M_wilson);
	Opencl_Module::print_profiling(filename, gamma5);
	Opencl_Module::print_profiling(filename, M_tm_plus);
	Opencl_Module::print_profiling(filename, M_tm_minus);
	Opencl_Module::print_profiling(filename, gamma5_eo);
	Opencl_Module::print_profiling(filename, M_tm_sitediagonal);
	Opencl_Module::print_profiling(filename, M_tm_inverse_sitediagonal);
	Opencl_Module::print_profiling(filename, M_tm_sitediagonal_minus);
	Opencl_Module::print_profiling(filename, M_tm_inverse_sitediagonal_minus);
	Opencl_Module::print_profiling(filename, dslash_eo);
	Opencl_Module::print_profiling(filename, _dslash_eo_boundary);
	Opencl_Module::print_profiling(filename, _dslash_eo_inner);
	Opencl_Module::print_profiling(filename, dslash_AND_M_tm_inverse_sitediagonal_eo);
	Opencl_Module::print_profiling(filename, dslash_AND_M_tm_inverse_sitediagonal_minus_eo);
	Opencl_Module::print_profiling(filename, M_tm_sitediagonal_AND_gamma5_eo);
	Opencl_Module::print_profiling(filename, M_tm_sitediagonal_minus_AND_gamma5_eo);
	Opencl_Module::print_profiling(filename, saxpy_AND_gamma5_eo);
}
hardware::code::Fermions::Fermions(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device)
	: Opencl_Module(kernelParameters, device),
	  M_wilson(0),
	  gamma5(0),
	  M_tm_plus(0),
	  M_tm_minus(0),
	  gamma5_eo(0),
	  M_tm_sitediagonal(0),
	  M_tm_inverse_sitediagonal(0),
	  M_tm_sitediagonal_minus(0),
	  M_tm_inverse_sitediagonal_minus(0),
	  dslash_eo(0),
	  _dslash_eo_boundary(0),
          _dslash_eo_inner(0),
	  dslash_AND_M_tm_inverse_sitediagonal_eo(0),
	  dslash_AND_M_tm_inverse_sitediagonal_minus_eo(0),
	  M_tm_sitediagonal_AND_gamma5_eo(0),
	  M_tm_sitediagonal_minus_AND_gamma5_eo(0),
	  saxpy_AND_gamma5_eo(0)
{
	fill_kernels();
}

hardware::code::Fermions::~Fermions()
{
	clear_kernels();
}

ClSourcePackage hardware::code::Fermions::get_sources() const noexcept
{
	return sources;
}
