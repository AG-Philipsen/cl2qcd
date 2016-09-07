/** @file
 * Implementation of the hardware::code::Spinors_staggered class
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

#include "spinors_staggered.hpp"

#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include <cassert>

#include "flopUtilities.hpp"
#include "gaugefield.hpp"
#include "prng.hpp"
#include "spinors.hpp"

using namespace std;

void hardware::code::Spinors_staggered::fill_kernels()
{
	if(kernelParameters->getFermact() != common::action::rooted_stagg){
		throw Print_Error_Message("Fermions_staggered module asked to be built but action set not to rooted_stagg! Aborting... ", __FILE__, __LINE__);
	}
  
	basic_fermion_code = get_basic_sources() <<  "operations_geometry.cl" << "operations_complex.h"  << "types_fermions.h" << "operations_su3vec.cl";
	if(kernelParameters->getUseEo()) {
		basic_fermion_code = basic_fermion_code << "spinorfield_staggered_eo.cl";;
	} else {
		basic_fermion_code = basic_fermion_code << "spinorfield_staggered.cl";
	}
	
	ClSourcePackage prng_code = get_device()->getPrngCode()->get_sources();
	
	logger.debug() << "Creating Spinors_staggered kernels...";
	
	//Reductions are really small kernels, so few needed options loaded by hands
	global_squarenorm_reduction_stagg = createKernel("global_squarenorm_reduction") << ClSourcePackage("-I " + std::string(SOURCEDIR) + " -D _INKERNEL_" + ((kernelParameters->getPrecision() == 64) ? (std::string(" -D _USEDOUBLEPREC_") + " -D _DEVICE_DOUBLE_EXTENSION_KHR_") : "")) << "types.h" << "spinorfield_staggered_squarenorm_reduction.cl";
	scalar_product_reduction_stagg = createKernel("scalar_product_reduction") << ClSourcePackage("-I " + std::string(SOURCEDIR) + " -D _INKERNEL_" + ((kernelParameters->getPrecision() == 64) ? (std::string(" -D _USEDOUBLEPREC_") + " -D _DEVICE_DOUBLE_EXTENSION_KHR_") : "")) << "types.h" << "operations_complex.h" << "spinorfield_staggered_scalar_product_reduction.cl";
	scalar_product_real_reduction_stagg = createKernel("scalar_product_real_reduction") << ClSourcePackage("-I " + std::string(SOURCEDIR) + " -D _INKERNEL_" + ((kernelParameters->getPrecision() == 64) ? (std::string(" -D _USEDOUBLEPREC_") + " -D _DEVICE_DOUBLE_EXTENSION_KHR_") : "")) << "types.h" << "operations_complex.h" << "spinorfield_staggered_scalar_product_reduction.cl";

	//In staggered formulation either eo or non-eo kernels are built!!
	if(kernelParameters->getUseEo()){
        scalar_product_stagg = 0;
        set_zero_spinorfield_stagg = 0;
        set_cold_spinorfield_stagg = 0;
        set_gaussian_spinorfield_stagg = 0;
        sax_stagg = 0;
        saxpy_stagg = 0;
        saxpbypz_stagg = 0;
		//Functionalities to convert eo from/to non eo
		convert_from_eoprec_stagg = createKernel("convert_from_eoprec_staggered") << basic_fermion_code << "spinorfield_staggered_eo_convert.cl";
		convert_to_eoprec_stagg = createKernel("convert_to_eoprec_staggered") << basic_fermion_code << "spinorfield_staggered_eo_convert.cl";
		//Functionalities to switch from AoS to SoA and viceversa
		convert_staggered_field_to_SoA_eo = createKernel("convert_staggered_field_to_SoA_eo") << basic_fermion_code << "spinorfield_staggered_eo_convert.cl";
		convert_staggered_field_from_SoA_eo = createKernel("convert_staggered_field_from_SoA_eo") << basic_fermion_code << "spinorfield_staggered_eo_convert.cl";
		//Squarenorm
		global_squarenorm_stagg_eoprec = createKernel("global_squarenorm_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_squarenorm.cl";
		//Scalar Product
		scalar_product_stagg_eoprec = createKernel("scalar_product_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_scalar_product.cl";
		scalar_product_real_part_stagg_eoprec = createKernel("scalar_product_real_part_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_scalar_product_real_part.cl";
		//Setting fields
		set_zero_spinorfield_stagg_eoprec = createKernel("set_zero_spinorfield_stagg_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_set_zero.cl";
		set_cold_spinorfield_stagg_eoprec = createKernel("set_cold_spinorfield_stagg_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_set_cold.cl";
		set_gaussian_spinorfield_stagg_eoprec = createKernel("set_gaussian_spinorfield_stagg_eoprec") << basic_fermion_code << prng_code << "spinorfield_staggered_eo_gaussian.cl";
		//Fields algebra operations
		sax_cplx_stagg_eoprec = createKernel("sax_cplx_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_sax_cplx.cl";
		sax_real_stagg_eoprec = createKernel("sax_real_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_sax_real.cl";
		sax_real_vec_stagg_eoprec = createKernel("sax_real_vec_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_sax_real_vec.cl";
		saxpy_cplx_stagg_eoprec = createKernel("saxpy_cplx_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_saxpy_cplx.cl";
		saxpy_real_stagg_eoprec = createKernel("saxpy_real_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_saxpy_real.cl";
		saxpy_real_vec_stagg_eoprec = createKernel("saxpy_real_vec_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_saxpy_real_vec.cl";
		saxpby_cplx_stagg_eoprec = createKernel("saxpby_cplx_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_saxpby_cplx.cl";
		saxpby_real_stagg_eoprec = createKernel("saxpby_real_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_saxpby_real.cl";
		saxpby_real_vec_stagg_eoprec = createKernel("saxpby_real_vec_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_saxpby_real_vec.cl";
		saxpbypz_cplx_stagg_eoprec = createKernel("saxpbypz_cplx_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_saxpbypz_cplx.cl";
		//Fields algebra operations with non buffer arguments
		sax_cplx_arg_stagg_eoprec = createKernel("sax_cplx_arg_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_sax_cplx.cl";
		sax_real_arg_stagg_eoprec = createKernel("sax_real_arg_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_sax_real.cl";
		saxpy_cplx_arg_stagg_eoprec = createKernel("saxpy_cplx_arg_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_saxpy_cplx.cl";
		saxpy_real_arg_stagg_eoprec = createKernel("saxpy_real_arg_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_saxpy_real.cl";
		saxpby_cplx_arg_stagg_eoprec = createKernel("saxpby_cplx_arg_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_saxpby_cplx.cl";
		saxpby_real_arg_stagg_eoprec = createKernel("saxpby_real_arg_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_saxpby_real.cl";
		saxpbypz_cplx_arg_stagg_eoprec = createKernel("saxpbypz_cplx_arg_staggered_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_saxpbypz_cplx.cl";
		//Mixed kernels
		sax_vectorized_and_squarenorm_reduction = createKernel("sax_vectorized_and_squarenorm_reduction") << basic_fermion_code << "spinorfield_staggered_eo_sax_AND_squarenorm.cl";
		sax_vectorized_and_squarenorm_eoprec = createKernel("sax_vectorized_and_squarenorm_eoprec") << basic_fermion_code << "spinorfield_staggered_eo_sax_AND_squarenorm.cl";
	} else {
		convert_from_eoprec_stagg = 0;
		convert_to_eoprec_stagg = 0;
		convert_staggered_field_to_SoA_eo = 0;
		convert_staggered_field_from_SoA_eo = 0;
		global_squarenorm_stagg_eoprec = 0;
		scalar_product_stagg_eoprec = 0;
		scalar_product_real_part_stagg_eoprec = 0;
		set_zero_spinorfield_stagg_eoprec = 0;
		set_cold_spinorfield_stagg_eoprec = 0;
		sax_cplx_stagg_eoprec = 0;
		sax_real_stagg_eoprec = 0;
		sax_real_vec_stagg_eoprec = 0;
		saxpy_cplx_stagg_eoprec = 0;
		saxpy_real_stagg_eoprec = 0;
		saxpy_real_vec_stagg_eoprec = 0;
		saxpby_cplx_stagg_eoprec = 0;
		saxpby_real_stagg_eoprec = 0;
		saxpby_real_vec_stagg_eoprec = 0;
		saxpbypz_cplx_stagg_eoprec = 0;
		sax_cplx_arg_stagg_eoprec = 0;
		sax_real_arg_stagg_eoprec = 0;
		saxpy_cplx_arg_stagg_eoprec = 0;
		saxpy_real_arg_stagg_eoprec = 0;
		saxpby_cplx_arg_stagg_eoprec = 0;
		saxpby_real_arg_stagg_eoprec = 0;
		saxpbypz_cplx_arg_stagg_eoprec = 0;
		set_gaussian_spinorfield_stagg_eoprec = 0;
		sax_vectorized_and_squarenorm_reduction = 0;
		sax_vectorized_and_squarenorm_eoprec = 0;
		//Scalar Product
		scalar_product_stagg = createKernel("scalar_product_staggered") << basic_fermion_code << "spinorfield_staggered_scalar_product.cl";
		//Setting fields
		set_zero_spinorfield_stagg = createKernel("set_zero_spinorfield_stagg") << basic_fermion_code << "spinorfield_staggered_set_zero.cl";
		set_cold_spinorfield_stagg = createKernel("set_cold_spinorfield_stagg") << basic_fermion_code << "spinorfield_staggered_set_cold.cl";
		set_gaussian_spinorfield_stagg = createKernel("set_gaussian_spinorfield_stagg") << basic_fermion_code << prng_code << "spinorfield_staggered_gaussian.cl";
		//Fields algebra operations
		sax_stagg = createKernel("sax_staggered") << basic_fermion_code << "spinorfield_staggered_sax.cl";
		saxpy_stagg = createKernel("saxpy_staggered") << basic_fermion_code << "spinorfield_staggered_saxpy.cl";
		saxpbypz_stagg = createKernel("saxpbypz_staggered") << basic_fermion_code << "spinorfield_staggered_saxpbypz.cl";        
	}
	//Squarenorm non_eo always built because needed in tests for conversion from eo to non eo
	global_squarenorm_stagg = createKernel("global_squarenorm_staggered") << basic_fermion_code << "spinorfield_staggered_squarenorm.cl";
}

void hardware::code::Spinors_staggered::clear_kernels()
{
	cl_int clerr = CL_SUCCESS;
	
	logger.debug() << "Clearing Spinors_staggered kernels...";
	
	//Reductions
	clerr = clReleaseKernel(global_squarenorm_reduction_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(scalar_product_reduction_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(scalar_product_real_reduction_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	
	if(kernelParameters->getUseEo()){
		//Functionalities to convert eo from/to non eo
		clerr = clReleaseKernel(convert_from_eoprec_stagg);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(convert_to_eoprec_stagg);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		//Functionalities to switch from AoS to SoA and viceversa
		clerr = clReleaseKernel(convert_staggered_field_to_SoA_eo);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(convert_staggered_field_from_SoA_eo);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		//Squarenorm
		clerr = clReleaseKernel(global_squarenorm_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		//Scalar Product
		clerr = clReleaseKernel(scalar_product_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(scalar_product_real_part_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		//Setting fields
		clerr = clReleaseKernel(set_zero_spinorfield_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(set_cold_spinorfield_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(set_gaussian_spinorfield_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		//Fields algebra operations
		clerr = clReleaseKernel(sax_cplx_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(sax_real_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(sax_real_vec_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpy_cplx_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpy_real_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpy_real_vec_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpby_cplx_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpby_real_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpby_real_vec_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpbypz_cplx_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		//Fields algebra operations with non buffers arguments
		clerr = clReleaseKernel(sax_cplx_arg_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(sax_real_arg_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpy_cplx_arg_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpy_real_arg_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpby_cplx_arg_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpby_real_arg_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpbypz_cplx_arg_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		//Mixed kernels
		clerr = clReleaseKernel(sax_vectorized_and_squarenorm_reduction);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(sax_vectorized_and_squarenorm_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	} else {
		//Scalar Product
		clerr = clReleaseKernel(scalar_product_stagg);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		//Setting fields
		clerr = clReleaseKernel(set_zero_spinorfield_stagg);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(set_cold_spinorfield_stagg);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(set_gaussian_spinorfield_stagg);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		//Fields algebra operations
		clerr = clReleaseKernel(sax_stagg);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpy_stagg);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpbypz_stagg);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
    //Squarenorm always to be released since always built because needed in tests for conversion from eo to non eo
    clerr = clReleaseKernel(global_squarenorm_stagg);
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
}


void hardware::code::Spinors_staggered::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);

	// kernels that use random numbers must not exceed the size of the random state array
	if(kernel == set_gaussian_spinorfield_stagg
	   || kernel == set_gaussian_spinorfield_stagg_eoprec) {
		if(*gs > hardware::buffers::get_prng_buffer_size(get_device(), kernelParameters->getUseSameRndNumbers())) {
			*gs = hardware::buffers::get_prng_buffer_size(get_device(), kernelParameters->getUseSameRndNumbers());
			logger.trace() << "I changed gs without changing neither ls nor num_groups (in Spinors_staggered::get_work_sizes)!!!";
		}
	}
	
	//Query specific sizes for kernels if needed
	//Whenever ls id manually modified, it is crucial to modify num_groups accordingly!
	if(kernel == global_squarenorm_stagg || kernel == scalar_product_stagg
	   || kernel == global_squarenorm_stagg_eoprec || kernel == scalar_product_stagg_eoprec
	   || kernel == sax_vectorized_and_squarenorm_eoprec || kernel == scalar_product_real_part_stagg_eoprec) {
	  if(*ls > 64) {
	    *ls = 64;
	    *num_groups = (*gs)/(*ls);
	  }
	  return;
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void hardware::code::Spinors_staggered::global_squarenorm_reduction(const hardware::buffers::Plain<hmc_float> * out, const hardware::buffers::Plain<hmc_float> * tmp_buf) const
{
	cl_int clerr = clSetKernelArg(global_squarenorm_reduction_stagg, 0, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(global_squarenorm_reduction_stagg, 1, sizeof(cl_mem), tmp_buf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	cl_uint elems = tmp_buf->get_elements();
	clerr = clSetKernelArg(global_squarenorm_reduction_stagg, 2, sizeof(cl_uint), &elems);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(global_squarenorm_reduction_stagg, 1, 1);
}

void hardware::code::Spinors_staggered::set_float_to_global_squarenorm_device(const hardware::buffers::Plain<su3vec> * a, const hardware::buffers::Plain<hmc_float> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(global_squarenorm_stagg, &ls2, &gs2, &num_groups);

	hardware::buffers::Plain<hmc_float> tmp(num_groups, get_device());

	//set arguments
	int clerr = clSetKernelArg(global_squarenorm_stagg, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	//CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(global_squarenorm_stagg, 1, sizeof(cl_mem), tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(global_squarenorm_stagg, 2, sizeof(hmc_float) * ls2, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(global_squarenorm_stagg, gs2, ls2);

	global_squarenorm_reduction(out, &tmp);
}


void hardware::code::Spinors_staggered::set_complex_to_scalar_product_device(const hardware::buffers::Plain<su3vec> * a, const hardware::buffers::Plain<su3vec> * b, const hardware::buffers::Plain<hmc_complex> * out) const
{
  //query work-sizes for kernel
  size_t ls2, gs2;
  cl_uint num_groups;
  this->get_work_sizes(scalar_product_stagg, &ls2, &gs2, &num_groups);

  hardware::buffers::Plain<hmc_complex> tmp(num_groups, get_device());

  //set arguments
  int clerr = clSetKernelArg(scalar_product_stagg, 0, sizeof(cl_mem), a->get_cl_buffer());
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  clerr = clSetKernelArg(scalar_product_stagg, 1, sizeof(cl_mem), b->get_cl_buffer());
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  clerr = clSetKernelArg(scalar_product_stagg, 2, sizeof(cl_mem), tmp);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  clerr = clSetKernelArg(scalar_product_stagg, 3, sizeof(hmc_complex) * ls2, static_cast<void*>(nullptr));
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  get_device()->enqueue_kernel(scalar_product_stagg , gs2, ls2);


  clerr = clSetKernelArg(scalar_product_reduction_stagg, 0, sizeof(cl_mem), tmp);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  clerr = clSetKernelArg(scalar_product_reduction_stagg, 1, sizeof(cl_mem), out->get_cl_buffer());
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  clerr = clSetKernelArg(scalar_product_reduction_stagg, 2, sizeof(cl_uint), &num_groups);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  get_device()->enqueue_kernel(scalar_product_reduction_stagg, 1, 1);
}


void hardware::code::Spinors_staggered::set_zero_spinorfield_device(const hardware::buffers::Plain<su3vec> * x) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_zero_spinorfield_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_zero_spinorfield_stagg, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(set_zero_spinorfield_stagg, gs2, ls2);
}


void hardware::code::Spinors_staggered::set_cold_spinorfield_device(const hardware::buffers::Plain<su3vec> * x) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_cold_spinorfield_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_cold_spinorfield_stagg, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(set_cold_spinorfield_stagg, gs2, ls2);
}


void hardware::code::Spinors_staggered::sax_device(const hardware::buffers::Plain<su3vec> * x, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<su3vec> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(sax_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(sax_stagg, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_stagg, 1, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_stagg, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(sax_stagg, gs2, ls2);
}


void hardware::code::Spinors_staggered::saxpy_device(const hardware::buffers::Plain<su3vec> * x, const hardware::buffers::Plain<su3vec> * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<su3vec> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy_stagg, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_stagg, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_stagg, 2, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_stagg, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpy_stagg , gs2, ls2);
}


void hardware::code::Spinors_staggered::saxpbypz_device(const hardware::buffers::Plain<su3vec> * x, const hardware::buffers::Plain<su3vec> * y, const hardware::buffers::Plain<su3vec> * z, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<hmc_complex> * beta, const hardware::buffers::Plain<su3vec> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpbypz_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpbypz_stagg, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_stagg, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_stagg, 2, sizeof(cl_mem), z->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_stagg, 3, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_stagg, 4, sizeof(cl_mem), beta->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_stagg, 5, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpbypz_stagg, gs2, ls2);
}


void hardware::code::Spinors_staggered::set_gaussian_spinorfield_device(const hardware::buffers::Plain<su3vec> * in, const hardware::buffers::PRNGBuffer * prng) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_gaussian_spinorfield_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_gaussian_spinorfield_stagg, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(set_gaussian_spinorfield_stagg, 1, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(set_gaussian_spinorfield_stagg  , gs2, ls2);

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> force_tmp(1, get_device());
		hmc_float resid;
		get_device()->getSpinorStaggeredCode()->set_float_to_global_squarenorm_device(in, &force_tmp);
		force_tmp.dump(&resid);
		logger.debug() <<  "\tinit gaussian spinorfield:\t" << resid;
		if(resid != resid) {
			throw Print_Error_Message("calculation of gaussian spinorfield gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}


void hardware::code::Spinors_staggered::convert_staggered_field_to_SoA_eo_device(const hardware::buffers::SU3vec * out, const hardware::buffers::Plain<su3vec> * in) const
{
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convert_staggered_field_to_SoA_eo, &ls2, &gs2, &num_groups);

	//set arguments
	int clerr = clSetKernelArg(convert_staggered_field_to_SoA_eo, 0, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_staggered_field_to_SoA_eo, 1, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	get_device()->enqueue_kernel(convert_staggered_field_to_SoA_eo, gs2, ls2);
}


void hardware::code::Spinors_staggered::convert_staggered_field_from_SoA_eo_device(const hardware::buffers::Plain<su3vec> * out, const hardware::buffers::SU3vec * in) const
{
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convert_staggered_field_from_SoA_eo, &ls2, &gs2, &num_groups);

	//set arguments
	int clerr = clSetKernelArg(convert_staggered_field_from_SoA_eo, 0, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_staggered_field_from_SoA_eo, 1, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(convert_staggered_field_from_SoA_eo, gs2, ls2);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////           EVEN-ODD PRECONDITIONING             //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hardware::code::Spinors_staggered::convert_from_eoprec_device(const hardware::buffers::SU3vec * in1, const hardware::buffers::SU3vec * in2, const hardware::buffers::Plain<su3vec> * out) const
{
	using namespace hardware::buffers;

	const size_4 mem_size = get_device()->getLocalLatticeMemoryExtents();

	// check buffer sizes
	const size_t in_size = in1->get_elements();
	if(in_size != get_eoprec_spinorfieldsize(mem_size) || in2->get_elements() != in_size) {
		throw std::invalid_argument("input buffers must be of size VOL4D / 2");
	}
	const size_t out_size = out->get_elements();
	if(out_size != get_spinorfieldsize(mem_size)) {
		throw std::invalid_argument("output buffer must be of size VOL4D");
	}

	const hardware::buffers::Buffer * tmp1, * tmp2;
	if(in1->is_soa()) {
		Plain<su3vec> * tmp = new Plain<su3vec>(in_size, get_device());
		convert_staggered_field_from_SoA_eo_device(tmp, in1);
		tmp1 = tmp;
	} else {
		tmp1 = in1;
	}
	if(in2->is_soa()) {
		Plain<su3vec> * tmp = new Plain<su3vec>(in_size, get_device());
		convert_staggered_field_from_SoA_eo_device(tmp, in2);
		tmp2 = tmp;
	} else {
		tmp2 = in2;
	}

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convert_from_eoprec_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(convert_from_eoprec_stagg, 0, sizeof(cl_mem), tmp1->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_from_eoprec_stagg, 1, sizeof(cl_mem), tmp2->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_from_eoprec_stagg, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(convert_from_eoprec_stagg, gs2, ls2);

	if(tmp1 != in1) {
		delete tmp1;
	}
	if(tmp2 != in2) {
		delete tmp2;
	}
}

void hardware::code::Spinors_staggered::convert_to_eoprec_device(const hardware::buffers::SU3vec * out1, const hardware::buffers::SU3vec * out2, const hardware::buffers::Plain<su3vec> * in) const
{
	using namespace hardware::buffers;

	const size_4 mem_size = get_device()->getLocalLatticeMemoryExtents();

	// check buffer sizes
	const size_t out_size = out1->get_elements();
	if(out_size != get_eoprec_spinorfieldsize(mem_size) || out2->get_elements() != out_size) {
		throw std::invalid_argument("output buffers must be of size VOL4D / 2");
	}
	const size_t in_size = in->get_elements();
	if(in_size != get_spinorfieldsize(mem_size)) {
		throw std::invalid_argument("input buffer must be of size VOL4D");
	}

	const hardware::buffers::Buffer * tmp1, * tmp2;
	if(out1->is_soa()) {
		tmp1 = new Plain<su3vec>(out_size, get_device());
	} else {
		tmp1 = out1;
	}
	if(out2->is_soa()) {
		tmp2 = new Plain<su3vec>(out_size, get_device());
	} else {
		tmp2 = out2;
	}

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convert_to_eoprec_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(convert_to_eoprec_stagg, 0, sizeof(cl_mem), tmp1->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_to_eoprec_stagg, 1, sizeof(cl_mem), tmp2->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_to_eoprec_stagg, 2, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(convert_to_eoprec_stagg , gs2, ls2);

	if(out1->is_soa()) {
		convert_staggered_field_to_SoA_eo_device(out1, static_cast<const Plain<su3vec> *>(tmp1));
		delete tmp1;
	}
	if(out2->is_soa()) {
		convert_staggered_field_to_SoA_eo_device(out2, static_cast<const Plain<su3vec> *>(tmp2));
		delete tmp2;
	}
}


void hardware::code::Spinors_staggered::set_float_to_global_squarenorm_eoprec_device(const hardware::buffers::SU3vec * a, const hardware::buffers::Plain<hmc_float> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(global_squarenorm_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	hardware::buffers::Plain<hmc_float> tmp(num_groups, get_device());

	int clerr = clSetKernelArg(global_squarenorm_stagg_eoprec, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(global_squarenorm_stagg_eoprec, 1, sizeof(cl_mem), tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(global_squarenorm_stagg_eoprec, 2, sizeof(hmc_float) * ls2, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( global_squarenorm_stagg_eoprec, gs2, ls2);

	global_squarenorm_reduction(out, &tmp);
}


void hardware::code::Spinors_staggered::set_complex_to_scalar_product_eoprec_device(const hardware::buffers::SU3vec * a, const hardware::buffers::SU3vec * b, const hardware::buffers::Plain<hmc_complex> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(scalar_product_stagg_eoprec, &ls2, &gs2, &num_groups);

	hardware::buffers::Plain<hmc_complex> scalar_product_buf(num_groups, get_device());
	assert(scalar_product_buf.get_elements() == num_groups);

	//set arguments
	int clerr = clSetKernelArg(scalar_product_stagg_eoprec, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_stagg_eoprec, 1, sizeof(cl_mem), b->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_stagg_eoprec, 2, sizeof(cl_mem), scalar_product_buf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_stagg_eoprec, 3, sizeof(hmc_complex) * ls2, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(scalar_product_stagg_eoprec, gs2, ls2);

	clerr = clSetKernelArg(scalar_product_reduction_stagg, 0, sizeof(cl_mem), scalar_product_buf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_reduction_stagg, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_reduction_stagg, 2, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(scalar_product_reduction_stagg, 1, 1);
}

void hardware::code::Spinors_staggered::set_float_to_scalar_product_real_part_eoprec_device(const hardware::buffers::SU3vec * a, const hardware::buffers::SU3vec * b, const hardware::buffers::Plain<hmc_float> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(scalar_product_real_part_stagg_eoprec, &ls2, &gs2, &num_groups);

	hardware::buffers::Plain<hmc_float> scalar_product_buf(num_groups, get_device());
	assert(scalar_product_buf.get_elements() == num_groups);

	//set arguments
	int clerr = clSetKernelArg(scalar_product_real_part_stagg_eoprec, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_real_part_stagg_eoprec, 1, sizeof(cl_mem), b->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_real_part_stagg_eoprec, 2, sizeof(cl_mem), scalar_product_buf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_real_part_stagg_eoprec, 3, sizeof(hmc_float) * ls2, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(scalar_product_real_part_stagg_eoprec, gs2, ls2);

	clerr = clSetKernelArg(scalar_product_real_reduction_stagg, 0, sizeof(cl_mem), scalar_product_buf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_real_reduction_stagg, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_real_reduction_stagg, 2, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	get_device()->enqueue_kernel(scalar_product_real_reduction_stagg, 1, 1);
}

void hardware::code::Spinors_staggered::set_zero_spinorfield_eoprec_device(const hardware::buffers::SU3vec * x) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_zero_spinorfield_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_zero_spinorfield_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(set_zero_spinorfield_stagg_eoprec, gs2, ls2);
}


void hardware::code::Spinors_staggered::set_cold_spinorfield_eoprec_device(const hardware::buffers::SU3vec * x) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_cold_spinorfield_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_cold_spinorfield_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(set_cold_spinorfield_stagg_eoprec , gs2, ls2);
}


void hardware::code::Spinors_staggered::sax_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(sax_cplx_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(sax_cplx_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_cplx_stagg_eoprec, 1, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_cplx_stagg_eoprec, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(sax_cplx_stagg_eoprec, gs2, ls2);
}

void hardware::code::Spinors_staggered::sax_eoprec_device(const hardware::buffers::SU3vec * x, const hmc_complex alpha, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(sax_cplx_arg_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(sax_cplx_arg_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_cplx_arg_stagg_eoprec, 1, sizeof(hmc_float), &alpha.re);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(sax_cplx_arg_stagg_eoprec, 2, sizeof(hmc_float), &alpha.im);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_cplx_arg_stagg_eoprec, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(sax_cplx_arg_stagg_eoprec, gs2, ls2);
}

void hardware::code::Spinors_staggered::sax_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::Plain<hmc_float> * alpha, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(sax_real_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(sax_real_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_real_stagg_eoprec, 1, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_real_stagg_eoprec, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(sax_real_stagg_eoprec, gs2, ls2);
}

void hardware::code::Spinors_staggered::sax_eoprec_device(const hardware::buffers::SU3vec * x, const hmc_float alpha, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(sax_real_arg_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(sax_real_arg_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_real_arg_stagg_eoprec, 1, sizeof(hmc_float), &alpha);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_real_arg_stagg_eoprec, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(sax_real_arg_stagg_eoprec, gs2, ls2);
}

void hardware::code::Spinors_staggered::sax_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::Plain<hmc_float> * alpha, const int index_alpha, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(sax_real_vec_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(sax_real_vec_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_real_vec_stagg_eoprec, 1, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(sax_real_vec_stagg_eoprec, 2, sizeof(cl_int), &index_alpha);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_real_vec_stagg_eoprec, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(sax_real_vec_stagg_eoprec, gs2, ls2);
}


void hardware::code::Spinors_staggered::saxpy_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy_cplx_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy_cplx_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_cplx_stagg_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_cplx_stagg_eoprec, 2, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_cplx_stagg_eoprec, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpy_cplx_stagg_eoprec, gs2, ls2);
}

void hardware::code::Spinors_staggered::saxpy_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hardware::buffers::Plain<hmc_float> * alpha, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy_real_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy_real_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_real_stagg_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_real_stagg_eoprec, 2, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_real_stagg_eoprec, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpy_real_stagg_eoprec, gs2, ls2);
}

void hardware::code::Spinors_staggered::saxpy_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hmc_complex alpha, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy_cplx_arg_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy_cplx_arg_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_cplx_arg_stagg_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_cplx_arg_stagg_eoprec, 2, sizeof(hmc_float), &alpha.re);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(saxpy_cplx_arg_stagg_eoprec, 3, sizeof(hmc_float), &alpha.im);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_cplx_arg_stagg_eoprec, 4, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpy_cplx_arg_stagg_eoprec, gs2, ls2);
}

void hardware::code::Spinors_staggered::saxpy_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hmc_float alpha, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy_real_arg_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy_real_arg_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_real_arg_stagg_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_real_arg_stagg_eoprec, 2, sizeof(hmc_float), &alpha);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_real_arg_stagg_eoprec, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpy_real_arg_stagg_eoprec, gs2, ls2);
}

void hardware::code::Spinors_staggered::saxpy_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hardware::buffers::Plain<hmc_float> * alpha, const int index_alpha, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy_real_vec_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy_real_vec_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_real_vec_stagg_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_real_vec_stagg_eoprec, 2, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(saxpy_real_vec_stagg_eoprec, 3, sizeof(cl_int), &index_alpha);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_real_vec_stagg_eoprec, 4, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpy_real_vec_stagg_eoprec, gs2, ls2);
}


void hardware::code::Spinors_staggered::saxpby_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<hmc_complex> * beta, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpby_cplx_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpby_cplx_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_cplx_stagg_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_cplx_stagg_eoprec, 2, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_cplx_stagg_eoprec, 3, sizeof(cl_mem), beta->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_cplx_stagg_eoprec, 4, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpby_cplx_stagg_eoprec, gs2, ls2);
}

void hardware::code::Spinors_staggered::saxpby_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hardware::buffers::Plain<hmc_float> * alpha, const hardware::buffers::Plain<hmc_float> * beta, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpby_real_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpby_real_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_real_stagg_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_real_stagg_eoprec, 2, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_real_stagg_eoprec, 3, sizeof(cl_mem), beta->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_real_stagg_eoprec, 4, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpby_real_stagg_eoprec, gs2, ls2);
}

void hardware::code::Spinors_staggered::saxpby_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hmc_complex alpha, const hmc_complex beta, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpby_cplx_arg_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpby_cplx_arg_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_cplx_arg_stagg_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_cplx_arg_stagg_eoprec, 2, sizeof(hmc_float), &alpha.re);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(saxpby_cplx_arg_stagg_eoprec, 3, sizeof(hmc_float), &alpha.im);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_cplx_arg_stagg_eoprec, 4, sizeof(hmc_float), &beta.re);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(saxpby_cplx_arg_stagg_eoprec, 5, sizeof(hmc_float), &beta.im);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_cplx_arg_stagg_eoprec, 6, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpby_cplx_arg_stagg_eoprec, gs2, ls2);
}

void hardware::code::Spinors_staggered::saxpby_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hmc_float alpha, const hmc_float beta, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpby_real_arg_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpby_real_arg_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_real_arg_stagg_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_real_arg_stagg_eoprec, 2, sizeof(hmc_float), &alpha);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_real_arg_stagg_eoprec, 3, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(saxpby_real_arg_stagg_eoprec, 4, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpby_real_arg_stagg_eoprec, gs2, ls2);
}

void hardware::code::Spinors_staggered::saxpby_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hardware::buffers::Plain<hmc_float> * alpha, const hardware::buffers::Plain<hmc_float> * beta, const int index_alpha, const int index_beta, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpby_real_vec_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpby_real_vec_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_real_vec_stagg_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_real_vec_stagg_eoprec, 2, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_real_vec_stagg_eoprec, 3, sizeof(cl_mem), beta->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_real_vec_stagg_eoprec, 4, sizeof(cl_int), &index_alpha);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(saxpby_real_vec_stagg_eoprec, 5, sizeof(cl_int), &index_beta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpby_real_vec_stagg_eoprec, 6, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpby_real_vec_stagg_eoprec, gs2, ls2);
}


void hardware::code::Spinors_staggered::saxpbypz_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hardware::buffers::SU3vec * z, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<hmc_complex> * beta, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpbypz_cplx_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpbypz_cplx_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_cplx_stagg_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_cplx_stagg_eoprec, 2, sizeof(cl_mem), z->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_cplx_stagg_eoprec, 3, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_cplx_stagg_eoprec, 4, sizeof(cl_mem), beta->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_cplx_stagg_eoprec, 5, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpbypz_cplx_stagg_eoprec, gs2, ls2);
}

void hardware::code::Spinors_staggered::saxpbypz_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hardware::buffers::SU3vec * z, const hmc_complex alpha, const hmc_complex beta, const hardware::buffers::SU3vec * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpbypz_cplx_arg_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpbypz_cplx_arg_stagg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_cplx_arg_stagg_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_cplx_arg_stagg_eoprec, 2, sizeof(cl_mem), z->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_cplx_arg_stagg_eoprec, 3, sizeof(hmc_float), &alpha.re);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(saxpbypz_cplx_arg_stagg_eoprec, 4, sizeof(hmc_float), &alpha.im);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_cplx_arg_stagg_eoprec, 5, sizeof(hmc_float), &beta.re);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_cplx_arg_stagg_eoprec, 6, sizeof(hmc_float), &beta.im);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(saxpbypz_cplx_arg_stagg_eoprec, 7, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpbypz_cplx_arg_stagg_eoprec, gs2, ls2);
}


void hardware::code::Spinors_staggered::set_gaussian_spinorfield_eoprec_device(const hardware::buffers::SU3vec
* in, const hardware::buffers::PRNGBuffer * prng) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_gaussian_spinorfield_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_gaussian_spinorfield_stagg_eoprec, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(set_gaussian_spinorfield_stagg_eoprec, 1, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(set_gaussian_spinorfield_stagg_eoprec, gs2, ls2);

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> force_tmp(1, get_device());
		hmc_float resid;
		get_device()->getSpinorStaggeredCode()->set_float_to_global_squarenorm_eoprec_device(in, &force_tmp);
		force_tmp.dump(&resid);
		logger.debug() <<  "\tinit gaussian spinorfield energy:\t" << resid;
		if(resid != resid) {
			throw Print_Error_Message("calculation of gaussian spinorfield gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

}


void hardware::code::Spinors_staggered::sax_vectorized_squarenorm_reduction(const hardware::buffers::Plain<hmc_float> * out, const hardware::buffers::Plain<hmc_float> * tmp_buf, const int numeqs) const
{
  cl_int clerr = clSetKernelArg(sax_vectorized_and_squarenorm_reduction, 0, sizeof(cl_mem), out->get_cl_buffer());
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
  
  clerr = clSetKernelArg(sax_vectorized_and_squarenorm_reduction, 1, sizeof(cl_mem), tmp_buf->get_cl_buffer());
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
  
  cl_uint elems = tmp_buf->get_elements();
  clerr = clSetKernelArg(sax_vectorized_and_squarenorm_reduction, 2, sizeof(cl_uint), &elems);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
  
  clerr = clSetKernelArg(sax_vectorized_and_squarenorm_reduction, 3, sizeof(cl_int), &numeqs);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  get_device()->enqueue_kernel(sax_vectorized_and_squarenorm_reduction, 1, 1);
}

void hardware::code::Spinors_staggered::sax_vectorized_and_squarenorm_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::Plain<hmc_float> * alpha, const int numeqs, const hardware::buffers::Plain<hmc_float> * out) const
{
	if(numeqs > std::max(kernelParameters->getMetroApproxOrd(), kernelParameters->getMdApproxOrd())){
	  throw Invalid_Parameters("In sax_vectorized_and_squarenorm_eoprec_device numeqs to big!!", "numeqs < " + to_string(std::max(kernelParameters->getMetroApproxOrd(), kernelParameters->getMdApproxOrd())), numeqs);
	}
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(sax_vectorized_and_squarenorm_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(sax_vectorized_and_squarenorm_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_vectorized_and_squarenorm_eoprec, 1, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_vectorized_and_squarenorm_eoprec, 2, sizeof(cl_int), &numeqs);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	hardware::buffers::Plain<hmc_float> tmp(num_groups*numeqs, get_device());
	clerr = clSetKernelArg(sax_vectorized_and_squarenorm_eoprec, 3, sizeof(cl_mem), tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(sax_vectorized_and_squarenorm_eoprec, 4, sizeof(hmc_float) * ls2 * numeqs, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	get_device()->enqueue_kernel(sax_vectorized_and_squarenorm_eoprec, gs2, ls2);
	//get_device()->enqueue_kernel(sax_vectorized_and_squarenorm_eoprec, 1, 1);
	
	sax_vectorized_squarenorm_reduction(out, &tmp, numeqs);

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


size_t hardware::code::Spinors_staggered::get_read_write_size(const std::string& in) const
{
	//Depending on the compile-options, one has different sizes...
	size_t D = kernelParameters->getFloatSize();
	size_t S = kernelParameters->getSpinorFieldSize();
	size_t Seo = kernelParameters->getEoprecSpinorFieldSize();
	//factor for complex numbers
	int C = 2;
	//NOTE: 1 spinor has NC*NSPIN = 3*1 = 3 complex entries
	if (in == "global_squarenorm_staggered") {
		//this kernel reads 1 su3vec and writes 1 real number
		/// @NOTE: here, the local reduction is not taken into account
		return D * S * (C * NC + 1);
	}
	if (in == "global_squarenorm_reduction") {
		//this kernel reads NUM_GROUPS real numbers and writes 1 real number
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(global_squarenorm_stagg, &ls2, &gs2, &num_groups);
		return D * (num_groups + 1);
	}
	if (in == "scalar_product_staggered") {
		//this kernel reads 2 spinors and writes 1 complex number
		/// @NOTE: here, the local reduction is not taken into account
		return C * D * S * (NC * 2 + 1);
	}
	if (in == "scalar_product_reduction") {
		//this kernel reads NUM_GROUPS complex numbers and writes 1 complex number
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(scalar_product_stagg, &ls2, &gs2, &num_groups);
		return C * D * (num_groups + 1);
	}
	if (in == "scalar_product_real_reduction") {
		//this kernel reads NUM_GROUPS complex numbers and writes 1 real number
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(scalar_product_stagg, &ls2, &gs2, &num_groups);
		return D * (C * num_groups + 1);
	}
	if (in == "set_zero_spinorfield_stagg") {
		//this kernel writes 1 su3vec
		return C * D * S * NC;
	}
	if (in == "set_cold_spinorfield_stagg") {
		//this kernel writes 1 su3vec
		return C * D * S * NC;
	}
	if (in == "set_gaussian_spinorfield_stagg") {
		//this kernel writes 1 su3vec per site
		return ( NC * C ) * D * S;
	}
	if (in == "sax_staggered") {
		//this kernel reads 1 su3vec, 1 complex number and writes 1 su3vec per site
		return C * D * S * (NC * (1 + 1) + 1);
	}
	if (in == "saxpy_staggered") {
		//this kernel reads 2 su3vec, 2 complex number and writes 1 su3vec per site
		return C * D * S * (NC * (2 + 1) + 2);
	}
	if (in == "saxpbypz_staggered") {
		//this kernel reads 3 su3vec, 2 complex number and writes 1 su3vec per site
		return C * D * S * (NC * (3 + 1) + 2);
	}
	if(in == "convert_staggered_field_to_SoA_eo") {
		//this kernel reads 1 su3vec and writes 1 su3vec per site (eo)
		return C * Seo * D * NC * (1 + 1);
	}
	if(in == "convert_staggered_field_from_SoA_eo") {
		//this kernel reads 1 su3vec and writes 1 su3vec per site (eo)
		return C * Seo * D * NC * (1 + 1);
	}
	if (in == "convert_from_eoprec_staggered") {
		//this kernel reads 2 su3vec and writes 2 su3vec per site (eo)
		//(actually it write 1 su3vec per site but on the whole lattice so 2*Seo)
		return C * NC * D * Seo * (2 + 2);
	}
	if (in == "convert_to_eoprec_staggered") {
		//this kernel reads 2 su3vec and writes 2 su3vec per site (eo)
		//(actually it reads 1 su3vec per site but on the whole lattice so 2*Seo)
		return C * NC * D * Seo * (2 + 2);
	}
	if (in == "global_squarenorm_staggered_eoprec") {
		//this kernel reads 1 su3vec and writes 1 real number (eo)
		/// @NOTE: here, the local reduction is not taken into account
		return D * Seo * (C * NC  + 1);
	}
	if (in == "scalar_product_staggered_eoprec") {
		//this kernel reads 2 su3vec and writes 1 complex number (eo)
		/// @NOTE: here, the local reduction is not taken into account
		return C * D * Seo * (2 * NC  + 1);
	}
	if (in == "scalar_product_real_part_staggered_eoprec") {
		//this kernel reads 2 su3vec and writes 1 real number (eo)
		/// @NOTE: here, the local reduction is not taken into account
		return D * Seo * (C * 2 * NC  + 1);
	}
	if (in == "set_zero_spinorfield_stagg_eoprec") {
		//this kernel writes 1 su3vec per site (eo)
		return C * NC * D * Seo;
	}
	if (in == "set_cold_spinorfield_stagg_eoprec") {
		//this kernel writes 1 su3vec per site (eo)
		return C * NC * D * Seo;
	}
	if (in == "sax_cplx_staggered_eoprec" || in == "sax_cplx_arg_staggered_eoprec") {
		//this kernel reads 1 su3vec, 1 complex number and writes 1 su3vec per site (eo)
		return C * D * Seo * (NC * (1 + 1) + 1);
	}
	
	if (in == "sax_real_staggered_eoprec" || in == "sax_real_arg_staggered_eoprec"
	                                      || in == "sax_real_vec_staggered_eoprec") {
		//this kernel reads 1 su3vec, 1 real number and writes 1 su3vec per site (eo)
		return D * Seo * (C * NC * (1 + 1) + 1);
	}
	
	if (in == "saxpy_cplx_staggered_eoprec" || in == "saxpy_cplx_arg_staggered_eoprec") {
		//this kernel reads 2 su3vec, 1 complex number and writes 1 su3vec per site (eo)
		return C * D * Seo * (NC * (2 + 1) + 1);
	}
	if (in == "saxpy_real_staggered_eoprec" || in == "saxpy_real_arg_staggered_eoprec"
	                                        || in == "saxpy_real_vec_staggered_eoprec") {
		//this kernel reads 2 su3vec, 1 real number and writes 1 su3vec per site (eo)
		return D * Seo * (C * NC * (2 + 1) + 1);
	}
	if (in == "saxpby_cplx_staggered_eoprec" || in == "saxpby_cplx_arg_staggered_eoprec") {
		//this kernel reads 2 su3vec, 2 complex number and writes 1 su3vec per site (eo)
		return C * D * Seo * (NC * (2 + 1) + 2);
	}
	if (in == "saxpby_real_staggered_eoprec" || in == "saxpby_real_arg_staggered_eoprec"
	                                         || in == "saxpby_real_vec_staggered_eoprec") {
		//this kernel reads 2 su3vec, 2 complex number and writes 1 su3vec per site (eo)
		return D * Seo * (C * NC * (2 + 1) + 2);
	}
	if (in == "saxpbypz_cplx_staggered_eoprec" || in == "saxpbypz_cplx_arg_staggered_eoprec") {
		//this kernel reads 3 su3vec, 2 complex number and writes 1 su3vec per site (eo)
		return C * D * Seo * (NC * (3 + 1) + 2);
	}
	if (in == "set_gaussian_spinorfield_stagg_eoprec") {
		//this kernel writes 1 su3vec
		return NC * C * D * Seo;
	}
	if (in == "sax_vectorized_and_squarenorm_eoprec" || in == "sax_vectorized_and_squarenorm_reduction") {
		//This if should not be entered since here we do not have the number of eqs.
		//The calculation of bytes should be done in the cgm or wherever the kernel is used
		throw Print_Error_Message("Spinors_staggered::get_read_write_size entered with in=sax_vectorized_and_squarenorm_*!", __FILE__, __LINE__);
	}
	
	logger.warn() << "No if entered in get_read_write_size(), in = " << in << ". Returning 0 bytes...";
	return 0;
}

uint64_t hardware::code::Spinors_staggered::get_flop_size(const std::string& in) const
{
	uint64_t S = kernelParameters->getSpinorFieldSize();
	uint64_t Seo = kernelParameters->getEoprecSpinorFieldSize();
	//this is the same as in the function above
	if (in == "global_squarenorm_staggered") {
		//this kernel performs spinor_squarenorm on each site and then adds S-1 complex numbers
		//Note that the sum of S numbers can be done in several way: the least efficient way is
		//to add the first 2 numbers, then the third, then the fourth and so on, performing S-1
		//additions. This is not what is done in the code but it is a good estimation because
		//we know for sure that the code will be a bit faster (it is somehow a boundary estimate)
		return S * getFlopSu3VecSquareNorm() + (S - 1) * 2;
	}
	if (in == "global_squarenorm_reduction") {
		//This if should not be entered since the sum of the site squarenorms
		//has already taken into account with the (S-1)*2 term in the previous if
		throw Print_Error_Message("Spinors_staggered::get_flop_size entered with in=global_squarenorm_reduction!", __FILE__, __LINE__);
	}
	if (in == "scalar_product_staggered") {
		//this kernel performs su3vec_scalarproduct on each site and then adds S-1 complex numbers
		//Note that the sum of S numbers can be done in several way: the least efficient way is
		//to add the first 2 numbers, then the third, then the fourth and so on, performing S-1
		//additions. This is not what is done in the code but it is a good estimation because
		//we know for sure that the code will be a bit faster (it is somehow a boundary estimate)
		return S * getFlopSu3VecTimesSu3Vec() + (S - 1) * 2;
	}
	if (in == "scalar_product_reduction" || in == "scalar_product_real_reduction") {
		//This if should not be entered since the sum of the site squarenorms
		//has already taken into account with the (S-1)*2 term in the previous if
		throw Print_Error_Message("Spinors_staggered::get_flop_size entered with in=scalar_product_reduction!", __FILE__, __LINE__);
	}
	if (in == "set_zero_spinorfield_stagg") {
		//this kernel does not do any flop
		return 0;
	}
	if (in == "set_cold_spinorfield_stagg") {
		//this kernel performs 1. / sqrt((3.f * VOL4D)) and su3vec_times_real for each site
		return S * (3 + NC * 2);
	}
	if (in == "set_gaussian_spinorfield_stagg") {
		//this kernel performs NC multiplications per site
		///@todo I did not count the gaussian normal pair production, which is very complicated...
		return NC * S;
	}
	if (in == "sax_staggered") {
		//this kernel performs on each site su3vec_times_complex
		return S * (NC * (getFlopComplexMult()));
	}
	if (in == "saxpy_staggered") {
		//this kernel performs on each site su3vec_times_complex and su3vec_acc
		return S * (NC * (getFlopComplexMult() + 2));
	}
	if (in == "saxpbypz_staggered") {
		//this kernel performs on each 2*site su3vec_times_complex and 2*su3vec_acc
		return S * (NC * 2 * ( getFlopComplexMult() + 2));
	}
	if(in == "convert_staggered_field_to_SoA_eo") {
		//this kernel does not perform any flop, he just copies memory
		return 0;
	}
	if(in == "convert_staggered_field_from_SoA_eo") {
		//this kernel does not perform any flop, he just copies memory
		return 0;
	}
	if (in == "convert_from_eoprec_staggered") {
		//this kernel does not perform any flop, he just copies memory
		return 0;
	}
	if (in == "convert_to_eoprec_staggered") {
		//this kernel does not perform any flop, he just copies memory
		return 0;
	}
	if (in == "global_squarenorm_staggered_eoprec") {
		//this kernel performs su3vec_squarenorm on each site (eo) and then adds Seo-1 complex numbers
		return Seo * getFlopSu3VecSquareNorm() + (Seo - 1) * 2;
	}
	if (in == "scalar_product_staggered_eoprec") {
		//this kernel performs su3vec*su3vec on each site (eo) and then adds Seo-1 complex numbers
		return Seo *  getFlopSu3VecTimesSu3Vec() + (Seo - 1) * 2;
	}
	if (in == "scalar_product_real_part_staggered_eoprec") {
		//this kernel performs su3vec*su3vec (only real part!) on each site (eo) and adds Seo-1 real num.
		return Seo *  getFlopSu3VecTimesSu3Vec() / 2 + (Seo - 1);
	}
	if (in == "set_zero_spinorfield_stagg_eoprec") {
		//this kernel does not do any flop
		return 0;
	}
	if (in == "set_cold_spinorfield_stagg_eoprec") {
		//this kernel performs 1. / sqrt((3.f * VOL4D)) and su3vec_times_real for each site
		return Seo * ( 3 + NC * 2);
	}
	if (in == "sax_cplx_staggered_eoprec" || in == "sax_cplx_arg_staggered_eoprec") {
		//this kernel performs on each site (eo) su3vec_times_complex
		return Seo * (NC * (getFlopComplexMult()));
	}
	if (in == "sax_real_staggered_eoprec" || in == "sax_real_arg_staggered_eoprec" 
	                                      || in == "sax_real_vec_staggered_eoprec" ) {
		//this kernel performs on each site (eo) su3vec_times_real
		return Seo * NC;
	}
	if (in == "saxpy_cplx_staggered_eoprec" || in == "saxpy_cplx_arg_staggered_eoprec") {
		//this kernel performs on each site (eo) su3vec_times_complex and su3vec_acc
		return Seo * (NC * (getFlopComplexMult() + 2));
	}
	if (in == "saxpy_real_staggered_eoprec" || in == "saxpy_real_arg_staggered_eoprec"
	                                        || in == "saxpy_real_vec_staggered_eoprec") {
		//this kernel performs on each site (eo) su3vec_times_real and su3vec_acc
		return Seo * (NC * (2 + 2));
	}
	if (in == "saxpby_cplx_staggered_eoprec" || in == "saxpby_cplx_arg_staggered_eoprec") {
		//this kernel performs on each site (eo) 2*su3vec_times_complex and 1*su3vec_acc
		return Seo * (NC * (2 * getFlopComplexMult() + 2));
	}
	if (in == "saxpby_real_staggered_eoprec" || in == "saxpby_real_arg_staggered_eoprec"
	                                         || in == "saxpby_real_vec_staggered_eoprec") {
		//this kernel performs on each site (eo) 2*su3vec_times_real and 1*su3vec_acc
		return Seo * (NC * (4 + 2));
	}
	if (in == "saxpbypz_cplx_staggered_eoprec" || in == "saxpbypz_cplx_arg_staggered_eoprec") {
		//this kernel performs on each site (eo) 2*su3vec_times_complex and 2*su3vec_acc
		return Seo * (NC * 2 * (getFlopComplexMult() + 2));
	}
	if (in == "set_gaussian_spinorfield_stagg_eoprec") {
		//this kernel performs NC multiplications per site
		///@todo I did not count the gaussian normal pair production, which is very complicated...
		return NC * Seo;
	}
	if (in == "sax_vectorized_and_squarenorm_eoprec" || in == "sax_vectorized_and_squarenorm_reduction") {
		//This if should not be entered since here we do not have the number of eqs.
		//The calculation of flops should be done in the cgm or wherever the kernel is used
		throw Print_Error_Message("Spinors_staggered::get_flop_size entered with in=sax_vectorized_and_squarenorm_*!", __FILE__, __LINE__);
	}
	
	logger.warn() << "No if entered in get_flop_size(), in = " << in << ". Returning 0 flop...";
	return 0;
}

void hardware::code::Spinors_staggered::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, global_squarenorm_stagg);
	//Opencl_Module::print_profiling(filename, global_squarenorm_reduction_stagg);
	Opencl_Module::print_profiling(filename, scalar_product_stagg);
	//Opencl_Module::print_profiling(filename, scalar_product_reduction_stagg);
	//Opencl_Module::print_profiling(filename, scalar_product_real_reduction_stagg);
	Opencl_Module::print_profiling(filename, set_zero_spinorfield_stagg);
	Opencl_Module::print_profiling(filename, set_cold_spinorfield_stagg);
	Opencl_Module::print_profiling(filename, set_gaussian_spinorfield_stagg);
	Opencl_Module::print_profiling(filename, sax_stagg);
	Opencl_Module::print_profiling(filename, saxpy_stagg);
	Opencl_Module::print_profiling(filename, saxpbypz_stagg);
	Opencl_Module::print_profiling(filename, convert_staggered_field_to_SoA_eo);
	Opencl_Module::print_profiling(filename, convert_staggered_field_from_SoA_eo);
	Opencl_Module::print_profiling(filename, convert_from_eoprec_stagg);
	Opencl_Module::print_profiling(filename, convert_from_eoprec_stagg);
	Opencl_Module::print_profiling(filename, global_squarenorm_stagg_eoprec);
	Opencl_Module::print_profiling(filename, scalar_product_stagg_eoprec);
	Opencl_Module::print_profiling(filename, scalar_product_real_part_stagg_eoprec);
	Opencl_Module::print_profiling(filename, set_zero_spinorfield_stagg_eoprec);
	Opencl_Module::print_profiling(filename, set_cold_spinorfield_stagg_eoprec);
	Opencl_Module::print_profiling(filename, set_gaussian_spinorfield_stagg_eoprec);
	Opencl_Module::print_profiling(filename, sax_cplx_stagg_eoprec);
	Opencl_Module::print_profiling(filename, sax_real_stagg_eoprec);
	Opencl_Module::print_profiling(filename, saxpy_cplx_stagg_eoprec);
	Opencl_Module::print_profiling(filename, saxpy_real_stagg_eoprec);
	Opencl_Module::print_profiling(filename, saxpby_cplx_stagg_eoprec);
	Opencl_Module::print_profiling(filename, saxpby_real_stagg_eoprec);
	Opencl_Module::print_profiling(filename, saxpbypz_cplx_stagg_eoprec);
	Opencl_Module::print_profiling(filename, sax_cplx_arg_stagg_eoprec);
	Opencl_Module::print_profiling(filename, sax_real_arg_stagg_eoprec);
	Opencl_Module::print_profiling(filename, saxpy_cplx_arg_stagg_eoprec);
	Opencl_Module::print_profiling(filename, saxpy_real_arg_stagg_eoprec);
	Opencl_Module::print_profiling(filename, saxpby_cplx_arg_stagg_eoprec);
	Opencl_Module::print_profiling(filename, saxpby_real_arg_stagg_eoprec);
	Opencl_Module::print_profiling(filename, saxpbypz_cplx_arg_stagg_eoprec);
	Opencl_Module::print_profiling(filename, sax_real_vec_stagg_eoprec);
	Opencl_Module::print_profiling(filename, saxpy_real_vec_stagg_eoprec);
	Opencl_Module::print_profiling(filename, saxpby_real_vec_stagg_eoprec);
	//Opencl_Module::print_profiling(filename, sax_vectorized_and_squarenorm_eoprec);
	//Opencl_Module::print_profiling(filename, sax_vectorized_and_squarenorm_reduction);
}

hardware::code::Spinors_staggered::Spinors_staggered(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device)
	: Opencl_Module(kernelParameters, device)
{
	fill_kernels();
}

hardware::code::Spinors_staggered::~Spinors_staggered()
{
	clear_kernels();
}

ClSourcePackage hardware::code::Spinors_staggered::get_sources() const noexcept
{
	return basic_fermion_code;
}

