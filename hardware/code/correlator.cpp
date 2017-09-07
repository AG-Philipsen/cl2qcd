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

#include "correlator.hpp"

#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include "spinors.hpp"
#include "prng.hpp"

using namespace std;

void hardware::code::Correlator::fill_kernels()
{
	basic_correlator_code = get_basic_sources() << "operations_geometry.cl" << "operations_complex.h" << "types_fermions.h" << "operations_su3vec.cl" << "operations_spinor.cl" << "spinorfield.cl";
	
	ClSourcePackage prng_code = get_device()->getPrngCode()->get_sources();

	logger.debug() << "Creating Correlator kernels...";

	if(kernelParameters->getSourceType() == common::point)
		create_point_source = createKernel("create_point_source") << basic_correlator_code << prng_code << "spinorfield_point_source.cl";
	else if (kernelParameters->getSourceType() == common::volume)
		create_volume_source = createKernel("create_volume_source") << basic_correlator_code << prng_code << "spinorfield_volume_source.cl";
	else if (kernelParameters->getSourceType() == common::timeslice)
		create_timeslice_source = createKernel("create_timeslice_source") << basic_correlator_code << prng_code << "spinorfield_timeslice_source.cl";
	else if (kernelParameters->getSourceType() == common::zslice)
		create_zslice_source = createKernel("create_zslice_source") << basic_correlator_code << prng_code << "spinorfield_zslice_source.cl";

	if(kernelParameters->getMeasureCorrelators() ) {
		//CP: If a pointsource is chosen, the correlators have a particular simple form.
		if(kernelParameters->getSourceType() == common::point) {
			switch (kernelParameters->getCorrDir()) {
				case 0 :
					correlator_ps = createKernel("correlator_ps_t") << basic_correlator_code << "fermionobservables/correlator_ps_point.cl";
					correlator_sc = createKernel("correlator_sc_t") << basic_correlator_code << "fermionobservables/correlator_sc_point.cl";
					correlator_vx = createKernel("correlator_vx_t") << basic_correlator_code << "fermionobservables/correlator_vx_point.cl";
					correlator_vy = createKernel("correlator_vy_t") << basic_correlator_code << "fermionobservables/correlator_vy_point.cl";
					correlator_vz = createKernel("correlator_vz_t") << basic_correlator_code << "fermionobservables/correlator_vz_point.cl";
					correlator_ax = createKernel("correlator_ax_t") << basic_correlator_code << "fermionobservables/correlator_ax_point.cl";
					correlator_ay = createKernel("correlator_ay_t") << basic_correlator_code << "fermionobservables/correlator_ay_point.cl";
					correlator_az = createKernel("correlator_az_t") << basic_correlator_code << "fermionobservables/correlator_az_point.cl";
					correlator_avps = createKernel("correlator_avps_t") << basic_correlator_code << "fermionobservables/correlator_avps_point.cl";
					break;
				case 3 :
					correlator_ps = createKernel("correlator_ps_z") << basic_correlator_code << "fermionobservables/correlator_ps_point.cl";
					correlator_sc = createKernel("correlator_sc_z") << basic_correlator_code << "fermionobservables/correlator_sc_point.cl";
					correlator_vx = createKernel("correlator_vx_z") << basic_correlator_code << "fermionobservables/correlator_vx_point.cl";
					correlator_vy = createKernel("correlator_vy_z") << basic_correlator_code << "fermionobservables/correlator_vy_point.cl";
					correlator_vz = createKernel("correlator_vz_z") << basic_correlator_code << "fermionobservables/correlator_vz_point.cl";
					correlator_ax = createKernel("correlator_ax_z") << basic_correlator_code << "fermionobservables/correlator_ax_point.cl";
					correlator_ay = createKernel("correlator_ay_z") << basic_correlator_code << "fermionobservables/correlator_ay_point.cl";
					correlator_az = createKernel("correlator_az_z") << basic_correlator_code << "fermionobservables/correlator_az_point.cl";
					correlator_avps = 0;
					break;
				default:
					stringstream errmsg;
					errmsg << "Could not create correlator kernel as correlator direction " << kernelParameters->getCorrDir() << " has not been implemented.";
					throw Print_Error_Message(errmsg.str());
			}
		} else {
			std::string filename_tmp =  "fermionobservables_correlators_stochastic.cl";
			switch (kernelParameters->getCorrDir()) {
				case 0 :
					correlator_ps = createKernel("correlator_ps_t") << basic_correlator_code << filename_tmp;
					correlator_sc = createKernel("correlator_sc_t") << basic_correlator_code << filename_tmp;
					correlator_vx = createKernel("correlator_vx_t") << basic_correlator_code << filename_tmp;
					correlator_vy = createKernel("correlator_vy_t") << basic_correlator_code << filename_tmp;
					correlator_vz = createKernel("correlator_vz_t") << basic_correlator_code << filename_tmp;
					correlator_ax = createKernel("correlator_ax_t") << basic_correlator_code << filename_tmp;
					correlator_ay = createKernel("correlator_ay_t") << basic_correlator_code << filename_tmp;
					correlator_az = createKernel("correlator_az_t") << basic_correlator_code << filename_tmp;
					break;
				case 3 :
					correlator_ps = createKernel("correlator_ps_z") << basic_correlator_code << filename_tmp;
					correlator_sc = createKernel("correlator_sc_z") << basic_correlator_code << filename_tmp;
					correlator_vx = createKernel("correlator_vx_z") << basic_correlator_code << filename_tmp;
					correlator_vy = createKernel("correlator_vy_z") << basic_correlator_code << filename_tmp;
					correlator_vz = createKernel("correlator_vz_z") << basic_correlator_code << filename_tmp;
					correlator_ax = createKernel("correlator_ax_z") << basic_correlator_code << filename_tmp;
					correlator_ay = createKernel("correlator_ay_z") << basic_correlator_code << filename_tmp;
					correlator_az = createKernel("correlator_az_z") << basic_correlator_code << filename_tmp;
					break;
				default:
					stringstream errmsg;
					errmsg << "Could not create correlator kernel as correlator direction " << kernelParameters->getCorrDir() << " has not been implemented.";
					throw Print_Error_Message(errmsg.str());
			}
		}
	}
}

void hardware::code::Correlator::clear_kernels()
{
	int clerr = CL_SUCCESS;
	
	logger.debug() << "Clearing Correlator kernels...";
	
	if(correlator_ps)
		clerr = clReleaseKernel(correlator_ps);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(correlator_sc)
		clerr = clReleaseKernel(correlator_sc);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(correlator_vx)
		clerr = clReleaseKernel(correlator_vx);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(correlator_vy)
		clerr = clReleaseKernel(correlator_vy);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(correlator_vz)
		clerr = clReleaseKernel(correlator_vz);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(correlator_ax)
		clerr = clReleaseKernel(correlator_ax);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(correlator_ay)
		clerr = clReleaseKernel(correlator_ay);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(correlator_az)
		clerr = clReleaseKernel(correlator_az);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(correlator_avps)
		clerr = clReleaseKernel(correlator_avps);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(create_point_source) {
		clerr = clReleaseKernel(create_point_source);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(create_volume_source) {
		clerr = clReleaseKernel(create_volume_source);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(create_timeslice_source) {
		clerr = clReleaseKernel(create_timeslice_source);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(create_zslice_source) {
		clerr = clReleaseKernel(create_zslice_source);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
}

void hardware::code::Correlator::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);

	//LZ: should be valid for all kernels for correlators, i.e. for names that look like correlator_??_?
	string kernelname = get_kernel_name(kernel);
	if( kernelname.find("correlator") == 0 ) {
		if(get_device()->get_device_type() == CL_DEVICE_TYPE_GPU) {
			*ls = kernelParameters->getNs();
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

cl_kernel hardware::code::Correlator::get_correlator_kernel(string which) const
{
	if( which.compare("ps") == 0 ) {
		return correlator_ps;
	}
	if( which.compare("sc") == 0 ) {
		return correlator_sc;
	}
	if( which.compare("vx") == 0 ) {
		return correlator_vx;
	}
	if( which.compare("vy") == 0 ) {
		return correlator_vy;
	}
	if( which.compare("vz") == 0 ) {
		return correlator_vz;
	}
	if( which.compare("ax") == 0 ) {
		return correlator_ax;
	}
	if( which.compare("ay") == 0 ) {
		return correlator_ay;
	}
	if( which.compare("az") == 0 ) {
		return correlator_az;
	}
	if( which.compare("avps") == 0 ) {
		return correlator_avps;
	}
	throw Print_Error_Message("get_correlator_kernel failed, no appropriate kernel found");
	return 0;
}

void hardware::code::Correlator::create_point_source_device(const hardware::buffers::Plain<spinor> * inout, int i, int spacepos, int timepos) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(create_point_source, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(create_point_source, 0, sizeof(cl_mem), inout->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_point_source, 1, sizeof(int), &i);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_point_source, 2, sizeof(int), &spacepos);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_point_source, 3, sizeof(int), &timepos);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( create_point_source, gs2, ls2);

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> sqn_tmp(1, get_device());
		hmc_float sqn;
		get_device()->getSpinorCode()->set_float_to_global_squarenorm_device(inout, &sqn_tmp);
		sqn_tmp.dump(&sqn);
		logger.debug() <<  "\t|source|^2:\t" << sqn;
		if(sqn != sqn) {
			throw Print_Error_Message("calculation of source gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

}

void hardware::code::Correlator::create_volume_source_device(const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::PRNGBuffer * prng) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(create_volume_source, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(create_volume_source, 0, sizeof(cl_mem), inout->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_volume_source, 1, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( create_volume_source, gs2, ls2);

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> sqn_tmp(1, get_device());
		hmc_float sqn;
		get_device()->getSpinorCode()->set_float_to_global_squarenorm_device(inout, &sqn_tmp);
		sqn_tmp.dump(&sqn);
		logger.debug() <<  "\t|source|^2:\t" << sqn;
		if(sqn != sqn) {
			throw Print_Error_Message("calculation of source gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}

void hardware::code::Correlator::create_timeslice_source_device(const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::PRNGBuffer * prng, const int timeslice) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(create_timeslice_source, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(create_timeslice_source, 0, sizeof(cl_mem), inout->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_timeslice_source, 1, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	int tmp = timeslice;
	clerr = clSetKernelArg(create_timeslice_source, 2, sizeof(int), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( create_timeslice_source, gs2, ls2);

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> sqn_tmp(1, get_device());
		hmc_float sqn;
		get_device()->getSpinorCode()->set_float_to_global_squarenorm_device(inout, &sqn_tmp);
		sqn_tmp.dump(&sqn);
		logger.debug() <<  "\t|source|^2:\t" << sqn;
		if(sqn != sqn) {
			throw Print_Error_Message("calculation of source gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}

void hardware::code::Correlator::create_zslice_source_device(const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::PRNGBuffer * prng, const int zslice) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(create_zslice_source, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(create_zslice_source, 0, sizeof(cl_mem), inout->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_zslice_source, 1, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	int tmp = zslice;
	clerr = clSetKernelArg(create_zslice_source, 2, sizeof(int), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( create_zslice_source, gs2, ls2);

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> sqn_tmp(1, get_device());
		hmc_float sqn;
		get_device()->getSpinorCode()->set_float_to_global_squarenorm_device(inout, &sqn_tmp);
		sqn_tmp.dump(&sqn);
		logger.debug() <<  "\t|source|^2:\t" << sqn;
		if(sqn != sqn) {
			throw Print_Error_Message("calculation of source gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}

void hardware::code::Correlator::correlator(const cl_kernel correlator_kernel, const hardware::buffers::Plain<hmc_float> * correlator, const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * source) const
{
	int clerr;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(correlator_kernel, &ls2, &gs2, &num_groups);
	//set arguments
	clerr = clSetKernelArg(correlator_kernel, 0, sizeof(cl_mem), correlator->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 1, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	if(source) {
		clerr = clSetKernelArg(correlator_kernel, 2, sizeof(cl_mem), source->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	}

	get_device()->enqueue_kernel(correlator_kernel , gs2, ls2);
}

void hardware::code::Correlator::correlator(const cl_kernel correlator_kernel, const hardware::buffers::Plain<hmc_float> * correlator, const hardware::buffers::Plain<spinor> * in1, const hardware::buffers::Plain<spinor> * in2, const hardware::buffers::Plain<spinor> * in3, const hardware::buffers::Plain<spinor> * in4) const
{
	int clerr;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(correlator_kernel, &ls2, &gs2, &num_groups);
	//set arguments
	clerr = clSetKernelArg(correlator_kernel, 0, sizeof(cl_mem), correlator->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 1, sizeof(cl_mem), in1->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 2, sizeof(cl_mem), in2->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 3, sizeof(cl_mem), in3->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 4, sizeof(cl_mem), in4->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(correlator_kernel , gs2, ls2);
}

void hardware::code::Correlator::correlator(const cl_kernel correlator_kernel, const hardware::buffers::Plain<hmc_float> * correlator, const hardware::buffers::Plain<spinor> * in1, const hardware::buffers::Plain<spinor> * source1, const hardware::buffers::Plain<spinor> * in2, const hardware::buffers::Plain<spinor> * source2, const hardware::buffers::Plain<spinor> * in3, const hardware::buffers::Plain<spinor> * source3, const hardware::buffers::Plain<spinor> * in4, const hardware::buffers::Plain<spinor> * source4) const
{
	int clerr;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(correlator_kernel, &ls2, &gs2, &num_groups);
	//set arguments
	clerr = clSetKernelArg(correlator_kernel, 0, sizeof(cl_mem), correlator->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 1, sizeof(cl_mem), in1->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 2, sizeof(cl_mem), source1->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 3, sizeof(cl_mem), in2->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 4, sizeof(cl_mem), source2->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 5, sizeof(cl_mem), in3->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 6, sizeof(cl_mem), source3->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 7, sizeof(cl_mem), in4->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 8, sizeof(cl_mem), source4->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(correlator_kernel , gs2, ls2);
}
size_t hardware::code::Correlator::get_read_write_size(const std::string& in) const
{
	//Depending on the compile-options, one has different sizes...
	size_t D = kernelParameters->getFloatSize();
	//this returns the number of entries in an su3-matrix
	size_t S = kernelParameters->getSpinorFieldSize();
	//factor for complex numbers
	int C = 2;
	//this is the same as in the function above
	//NOTE: 1 spinor has NC*NDIM = 12 complex entries
	if (in == "create_point_source") {
		return module_metric_not_implemented<size_t>();
	}
	if (in == "create_volume_source") {
		return module_metric_not_implemented<size_t>();
	}
	if (in == "create_timeslice_source") {
		return module_metric_not_implemented<size_t>();
	}
	if (in == "create_zslice_source") {
		return module_metric_not_implemented<size_t>();
	}
	if (in == "correlator_ps_z" ) {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = kernelParameters->getNumSources();
		if(kernelParameters->getCorrDir() == 3) size_buffer = kernelParameters->getNs();
		if(kernelParameters->getCorrDir() == 0) size_buffer = kernelParameters->getNt();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}
	if (in == "correlator_sc_z") {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = kernelParameters->getNumSources();
		if(kernelParameters->getCorrDir() == 3) size_buffer = kernelParameters->getNs();
		if(kernelParameters->getCorrDir() == 0) size_buffer = kernelParameters->getNt();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}
	if (in == "correlator_vx_z") {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = kernelParameters->getNumSources();
		if(kernelParameters->getCorrDir() == 3) size_buffer = kernelParameters->getNs();
		if(kernelParameters->getCorrDir() == 0) size_buffer = kernelParameters->getNt();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}
	if (in == "correlator_vy_z") {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = kernelParameters->getNumSources();
		if(kernelParameters->getCorrDir() == 3) size_buffer = kernelParameters->getNs();
		if(kernelParameters->getCorrDir() == 0) size_buffer = kernelParameters->getNt();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}
	if (in == "correlator_vz_z") {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = kernelParameters->getNumSources();
		if(kernelParameters->getCorrDir() == 3) size_buffer = kernelParameters->getNs();
		if(kernelParameters->getCorrDir() == 0) size_buffer = kernelParameters->getNt();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}
	if (in == "correlator_ax_z") {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = kernelParameters->getNumSources();
		if(kernelParameters->getCorrDir() == 3) size_buffer = kernelParameters->getNs();
		if(kernelParameters->getCorrDir() == 0) size_buffer = kernelParameters->getNt();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}
	if (in == "correlator_ay_z") {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = kernelParameters->getNumSources();
		if(kernelParameters->getCorrDir() == 3) size_buffer = kernelParameters->getNs();
		if(kernelParameters->getCorrDir() == 0) size_buffer = kernelParameters->getNt();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}
	if (in == "correlator_az_z") {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = kernelParameters->getNumSources();
		if(kernelParameters->getCorrDir() == 3) size_buffer = kernelParameters->getNs();
		if(kernelParameters->getCorrDir() == 0) size_buffer = kernelParameters->getNt();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}

	return 0;
}

uint64_t hardware::code::Correlator::get_flop_size(const std::string& in) const
{
	//this is the same as in the function above
	if (in == "create_point_source") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "create_volume_source") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "create_timeslice_source") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "create_zslice_source") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "correlator_ps_z" ) {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "correlator_sc_z") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "correlator_vx_z") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "correlator_vy_z") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "correlator_vz_z") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "correlator_ax_z") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "correlator_ay_z") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "correlator_az_z") {
		return module_metric_not_implemented<uint64_t>();
	}

	return 0;
}

void hardware::code::Correlator::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	if(create_point_source)
		Opencl_Module::print_profiling(filename, create_point_source);
	if(create_volume_source)
		Opencl_Module::print_profiling(filename, create_volume_source);
	if(create_timeslice_source)
		Opencl_Module::print_profiling(filename, create_timeslice_source);
	if(create_zslice_source)
		Opencl_Module::print_profiling(filename, create_zslice_source);
	if(correlator_ps)
		Opencl_Module::print_profiling(filename, correlator_ps);
	if(correlator_sc)
		Opencl_Module::print_profiling(filename, correlator_sc);
	if(correlator_vx)
		Opencl_Module::print_profiling(filename, correlator_vx);
	if(correlator_vy)
		Opencl_Module::print_profiling(filename, correlator_vy);
	if(correlator_vz)
		Opencl_Module::print_profiling(filename, correlator_vz);
	if(correlator_ax)
		Opencl_Module::print_profiling(filename, correlator_ax);
	if(correlator_ay)
		Opencl_Module::print_profiling(filename, correlator_ay);
	if(correlator_az)
		Opencl_Module::print_profiling(filename, correlator_az);
}

hardware::code::Correlator::Correlator(const hardware::code::OpenClKernelParametersInterface& kernelParameters , const hardware::Device * device)
	: Opencl_Module(kernelParameters, device),
	  create_point_source(0), create_volume_source(0), create_timeslice_source(0), create_zslice_source(0),
	  correlator_ps(0), correlator_sc(0), correlator_vx(0), correlator_vy(0), correlator_vz(0), correlator_ax(0), correlator_ay(0), correlator_az(0),
	  correlator_avps(0), pbp_std(0), pbp_tm_one_end(0)
{
	fill_kernels();
}

hardware::code::Correlator::~Correlator()
{
	clear_kernels();
}
