/* staggeredTwoFlavourCorrelators.cpp
 *
 *@file
 * Implementation of the staggered two flavour correlator calculation.
 *
 * Copyright 2016 Alessandro Sciarra, Tim Breitenfelder
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


#include "correlator_staggered.hpp"

#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include "spinors.hpp" //for hardware::code::get_eoprec_spinorfieldsize();
#include "spinors_staggered.hpp"
#include "prng.hpp"

using namespace std;

void hardware::code::Correlator_staggered::fill_kernels()
{
	if(kernelParameters->getFermact() != common::action::rooted_stagg){
		throw Print_Error_Message("Correlator_staggered module asked to be built but action set not to rooted_stagg! Aborting... ", __FILE__, __LINE__);
	}
	
	//When some noneo method will be introduced, remove this check!
	if(!kernelParameters->getUseEo()){
		throw Print_Error_Message("Correlator_staggered module asked to be built but without even-odd preconditionig! Aborting... ", __FILE__, __LINE__);
	}
  
	basic_correlator_code = get_basic_sources() << "operations_geometry.cl" << "operations_complex.h" << "types_fermions.h" << "operations_su3vec.cl" << "spinorfield_staggered_eo.cl";
	
	ClSourcePackage prng_code = get_device()->getPrngCode()->get_sources();

	logger.debug() << "Creating Correlator_staggered kernels...";

	correlator_staggered_ps = createKernel("correlator_staggered_ps") << basic_correlator_code << prng_code << "fermionobservables/correlator_staggered_ps.cl";

	if(kernelParameters->getSourceType() == common::point)
		create_point_source_stagg_eoprec = createKernel("create_point_source_stagg_eoprec") << basic_correlator_code << prng_code << "spinorfield_staggered_eo_point_source.cl";
	else if (kernelParameters->getSourceType() == common::volume)
		create_volume_source_stagg_eoprec = createKernel("create_volume_source_stagg_eoprec") << basic_correlator_code << prng_code << "spinorfield_staggered_eo_volume_source.cl";
	else if (kernelParameters->getSourceType() == common::timeslice)
		throw Print_Error_Message("Timeslice source not implemented in Correlator_staggered module! Aborting...", __FILE__, __LINE__);
	else if (kernelParameters->getSourceType() == common::zslice)
		throw Print_Error_Message("Zslice source not implemented in Correlator_staggered module! Aborting...", __FILE__, __LINE__);

}

void hardware::code::Correlator_staggered::clear_kernels()
{
	int clerr = CL_SUCCESS;
	
	logger.debug() << "Clearing Correlator_staggered kernels...";

	if(create_volume_source_stagg_eoprec) {
		clerr = clReleaseKernel(create_volume_source_stagg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}

	if(create_point_source_stagg_eoprec) {
			clerr = clReleaseKernel(create_point_source_stagg_eoprec);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		}

	if(correlator_staggered_ps) {
			clerr = clReleaseKernel(correlator_staggered_ps);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		}
}

void hardware::code::Correlator_staggered::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);

	//LZ: should be valid for all kernels for correlators, i.e. for names that look like correlator_??_?
	string kernelname = get_kernel_name(kernel);
	if( kernelname.find("correlator") == 0 ) {
		if(get_device()->get_device_type() == CL_DEVICE_TYPE_GPU) {
			*ls = kernelParameters->getNs(); //TODO: Doesn't a huge Ns break the code since ls has a maximum depending on the device?!
			*gs = *ls;
			*num_groups = 1;
		} else {
			*ls = 1;
			*gs = *ls;
			*num_groups = 1;
		}
	}

	return;
}

void hardware::code::Correlator_staggered::create_point_source_stagg_eoprec_device(const hardware::buffers::SU3vec * inout, int i, int spacepos, int timepos) const
{
	//Make check on "i" and throw if not 0,1,2
	if(i < 0 || i > 2) throw Print_Error_Message("Color index i of staggered spinor must be between 0 and 2!", __FILE__, __LINE__);
    //query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(create_point_source_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(create_point_source_stagg_eoprec, 0, sizeof(cl_mem), inout->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_point_source_stagg_eoprec, 1, sizeof(int), &i);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_point_source_stagg_eoprec, 2, sizeof(int), &spacepos);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_point_source_stagg_eoprec, 3, sizeof(int), &timepos);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( create_point_source_stagg_eoprec, gs2, ls2);

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> sqn_tmp(1, get_device());
		hmc_float sqn;
		get_device()->getSpinorStaggeredCode()->set_float_to_global_squarenorm_eoprec_device(inout, &sqn_tmp);
		sqn_tmp.dump(&sqn);
		logger.debug() <<  "\t|source|^2:\t" << sqn;
		if(sqn != sqn) {
			throw Print_Error_Message("calculation of source gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}

void hardware::code::Correlator_staggered::create_volume_source_stagg_eoprec_device(const hardware::buffers::SU3vec * inout, const hardware::buffers::PRNGBuffer * prng) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(create_volume_source_stagg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(create_volume_source_stagg_eoprec, 0, sizeof(cl_mem), inout->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_volume_source_stagg_eoprec, 1, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(create_volume_source_stagg_eoprec, gs2, ls2);

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> sqn_tmp(1, get_device());
		hmc_float sqn;
		get_device()->getSpinorStaggeredCode()->set_float_to_global_squarenorm_eoprec_device(inout, &sqn_tmp);
		sqn_tmp.dump(&sqn);
		logger.debug() <<  "\t|source|^2:\t" << sqn;
		if(sqn != sqn) {
			throw Print_Error_Message("calculation of source gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}

void hardware::code::Correlator_staggered::pseudoScalarCorrelator(const hardware::buffers::Plain<hmc_float> * correlator, const hardware::buffers::SU3vec * invertedSourceEven,
																														  const hardware::buffers::SU3vec * invertedSourceOdd) const
{
	int clerr;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;

	this->get_work_sizes(correlator_staggered_ps, &ls2, &gs2, &num_groups);
	//set arguments
	clerr = clSetKernelArg(correlator_staggered_ps, 0, sizeof(cl_mem), correlator->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_staggered_ps, 1, sizeof(cl_mem), invertedSourceEven->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_staggered_ps, 2, sizeof(cl_mem), invertedSourceOdd->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(correlator_staggered_ps , gs2, ls2);
}

size_t hardware::code::Correlator_staggered::get_read_write_size(const std::string& in) const
{
    //Depending on the compile-options, one has different sizes...
    const size_t D = kernelParameters->getFloatSize();
    //this returns the number of entries in an su3-matrix
    const size_t S = kernelParameters->getSpinorFieldSize();
    //factor for complex numbers
    const int C = 2;
    //NOTE: one staggered spinor has NCOLOR = 3 complex entries

    if (in == "create_point_source_stagg_eoprec") {
		return module_metric_not_implemented<size_t>();
	}
	if (in == "create_volume_source_stagg_eoprec") {
		return module_metric_not_implemented<size_t>();
	}
	if (in == "create_timeslice_source_stagg_eoprec") {
		return module_metric_not_implemented<size_t>();
	}
	if (in == "create_zslice_source_stagg_eoprec") {
		return module_metric_not_implemented<size_t>();
	}

	if (in == "correlator_staggered_ps" ) {
			//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
			int size_buffer = 0;
			int num_sources = kernelParameters->getNumSources();
			if(kernelParameters->getCorrDir() == 0) size_buffer = kernelParameters->getNt();
			if(kernelParameters->getCorrDir() == 3) throw Print_Error_Message("staggered correlator calculation in z-direction not implemented yet! Aborting...");
			return num_sources * S * D * 3 * C + size_buffer * D;
	}

	logger.warn() << "No if entered in get_read_write_size(), in = " << in << ". Returning 0 bytes...";
	return 0;
}

uint64_t hardware::code::Correlator_staggered::get_flop_size(const std::string& in) const
{
	//this is the same as in the function above
	if (in == "create_point_source_stagg_eoprec") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "create_volume_source_stagg_eoprec") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "create_timeslice_source_stagg_eoprec") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "create_zslice_source_stagg_eoprec") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "correlator_staggered_ps" ) {
		return module_metric_not_implemented<uint64_t>();
	}

	logger.warn() << "No if entered in get_flop_size(), in = " << in << ". Returning 0 flop...";
	return 0;
}

void hardware::code::Correlator_staggered::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	if(create_point_source_stagg_eoprec)
	        Opencl_Module::print_profiling(filename, create_point_source_stagg_eoprec);
	if(create_volume_source_stagg_eoprec)
		Opencl_Module::print_profiling(filename, create_volume_source_stagg_eoprec);

	if(correlator_staggered_ps)
		Opencl_Module::print_profiling(filename, correlator_staggered_ps);
}

hardware::code::Correlator_staggered::Correlator_staggered(const hardware::code::OpenClKernelParametersInterface& kernelParameters , const hardware::Device * device)
	: Opencl_Module(kernelParameters, device), create_point_source_stagg_eoprec(0), create_volume_source_stagg_eoprec(0), correlator_staggered_ps(0)
{
	fill_kernels();
}

hardware::code::Correlator_staggered::~Correlator_staggered()
{
	clear_kernels();
}
