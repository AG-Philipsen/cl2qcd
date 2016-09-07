/** @file
 * Implementation of the hardware::code::Gaugemomentum class
 *
 * Copyright (c) 2012 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
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

#include "gaugemomentum.hpp"

#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include <fstream>
#include <cmath>
#include "fermions.hpp"
#include "flopUtilities.hpp"
#include "spinors.hpp"
#include "prng.hpp"

using namespace std;

hardware::code::Gaugemomentum::Gaugemomentum(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device)
	: Opencl_Module(kernelParameters, device)
{
	fill_kernels();
}

hardware::code::Gaugemomentum::~Gaugemomentum()
{
	clear_kernels();
}

void hardware::code::Gaugemomentum::fill_kernels()
{
	basic_gaugemomentum_code = get_basic_sources() << "operations_geometry.cl" << "operations_complex.h" << "types_hmc.h" << "operations_gaugemomentum.cl";
	
	ClSourcePackage prng_code = get_device()->getPrngCode()->get_sources();
	
	logger.debug() << "Creating Gaugemomentum kernels...";

	gaugemomentum_squarenorm_reduction = createKernel("global_squarenorm_reduction")  << ClSourcePackage("-I " + std::string(SOURCEDIR) + " -D _INKERNEL_" + ((kernelParameters->getPrecision() == 64) ? (std::string(" -D _USEDOUBLEPREC_") + " -D _DEVICE_DOUBLE_EXTENSION_KHR_") : "")) << "types.h" << "gaugemomentum_squarenorm_reduction.cl";

	_set_zero_gaugemomentum = createKernel("set_zero_gaugemomentum") << basic_gaugemomentum_code <<  "gaugemomentum_zero.cl";
	generate_gaussian_gaugemomenta = createKernel("generate_gaussian_gaugemomenta") << basic_gaugemomentum_code << prng_code << "gaugemomentum_gaussian.cl";
	gaugemomentum_squarenorm = createKernel("gaugemomentum_squarenorm") << basic_gaugemomentum_code << "gaugemomentum_squarenorm.cl";
	gaugemomentum_saxpy = createKernel("gaugemomentum_saxpy") << basic_gaugemomentum_code  << "gaugemomentum_saxpy.cl";
	if(get_device()->get_prefers_soa()) {
		gaugemomentum_convert_to_soa = createKernel("gaugemomentum_convert_to_soa") << basic_gaugemomentum_code << "gaugemomentum_convert.cl";
		gaugemomentum_convert_from_soa = createKernel("gaugemomentum_convert_from_soa") << basic_gaugemomentum_code << "gaugemomentum_convert.cl";
	} else {
		gaugemomentum_convert_to_soa = 0;
		gaugemomentum_convert_from_soa = 0;
	}
}

void hardware::code::Gaugemomentum::clear_kernels()
{
	cl_uint clerr = CL_SUCCESS;
	
	logger.debug() << "Clearing Gaugemomentum kernels...";
	
	clerr = clReleaseKernel(gaugemomentum_squarenorm_reduction);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(generate_gaussian_gaugemomenta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(_set_zero_gaugemomentum);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(gaugemomentum_squarenorm);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(gaugemomentum_saxpy);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(get_device()->get_prefers_soa()) {
		clerr = clReleaseKernel(gaugemomentum_convert_to_soa);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(gaugemomentum_convert_from_soa);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
}

void hardware::code::Gaugemomentum::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);

	// kernels that use random numbers must not exceed the size of the random state array
	if(kernel == generate_gaussian_gaugemomenta) {
		if(*gs > hardware::buffers::get_prng_buffer_size(get_device(), kernelParameters->getUseSameRndNumbers())) {
			*gs = hardware::buffers::get_prng_buffer_size(get_device(), kernelParameters->getUseSameRndNumbers());
		}
	} else if(kernel == gaugemomentum_squarenorm) {
		if(*ls != 1) { 
		  *ls = 128;
		  *num_groups = *gs / *ls;
		}
	}

	return;
}

size_t hardware::code::Gaugemomentum::get_read_write_size(const std::string& in) const
{
	//Depending on the compile-options, one has different sizes...
	size_t D = kernelParameters->getFloatSize();
	//this is the number of links in the system (and of gaugemomenta)
	size_t G = kernelParameters->getLatticeVolume() * NDIM;
	//this is the same as in the function above
	//NOTE: 1 ae has NC*NC-1 = 8 real entries
	int A = getSu3AlgebraSize();
	if (in == "generate_gaussian_gaugemomenta") {
		//this kernel writes 1 ae
		return (A) * D * G;
	}
	if (in == "set_zero_gaugemomentum;") {
		//this kernel writes 1 ae per link
		return G * D * A;
	}
	if (in == "gaugemomentum_squarenorm") {
		//this kernel reads 1 ae and writes 1 float per link
		return G * D * ( A + 1);
	}
	if (in == "gaugemomentum_saxpy") {
		//this kernel reads 2 ae and writes 1 ae per link
		return ((2 + 1) * A) * D * G;
	}
	if (in == "gaugemomentum_convert_to_soa") {
		return module_metric_not_implemented<size_t>();
	}
	if (in == "gaugemomentum_convert_from_soa") {
		return module_metric_not_implemented<size_t>();
	}
	return 0;
}

uint64_t hardware::code::Gaugemomentum::get_flop_size(const std::string& in) const
{
	//this is the number of links in the system (and of gaugemomenta)
	uint64_t G = kernelParameters->getLatticeVolume() * NDIM;
	//NOTE: 1 ae has NC*NC-1 = 8 real entries
	uint64_t A = getSu3AlgebraSize();
	if (in == "generate_gaussian_gaugemomenta") {
		//this kernel performs 0 multiplications per site
		///@todo ? I did not count the gaussian normal pair production, which is very complicated...
		return 0;
	}
	if (in == "set_zero_gaugemomentum;") {
		//this kernel performs 0 mults
		return 0;
	}
	if (in == "gaugemomentum_squarenorm") {
		//this kernel performs 8 real mults and 8-1 real adds per ae
		return (8 + 7) * A * G;
	}
	if (in == "gaugemomentum_saxpy") {
		//this kernel performs 1 real mult and 1 real add per ae
		return (1 + 1) * A * G;
	}
	if (in == "gaugemomentum_convert_to_soa") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "gaugemomentum_convert_from_soa") {
		return module_metric_not_implemented<uint64_t>();
	}
	return 0;
}

void hardware::code::Gaugemomentum::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, generate_gaussian_gaugemomenta);
	Opencl_Module::print_profiling(filename, _set_zero_gaugemomentum);
	Opencl_Module::print_profiling(filename, gaugemomentum_squarenorm);
	Opencl_Module::print_profiling(filename, gaugemomentum_saxpy);
	Opencl_Module::print_profiling(filename, gaugemomentum_convert_to_soa);
	Opencl_Module::print_profiling(filename, gaugemomentum_convert_from_soa);
}

/////////////////////////////////////////////////
// Methods on device

void hardware::code::Gaugemomentum::generate_gaussian_gaugemomenta_device(const hardware::buffers::Gaugemomentum * in, const hardware::buffers::PRNGBuffer * prng) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(generate_gaussian_gaugemomenta, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(generate_gaussian_gaugemomenta, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(generate_gaussian_gaugemomenta, 1, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( generate_gaussian_gaugemomenta , gs2, ls2);

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> force_tmp(1, get_device());
		hmc_float resid;
		this->set_float_to_gaugemomentum_squarenorm_device(in, &force_tmp);
		force_tmp.dump(&resid);
		logger.debug() <<  "\tgaussian gaugemomenta:\t" << resid;
		if(resid != resid) {
			throw Print_Error_Message("calculation of gaussian gm gave nan! Aborting...", __FILE__, __LINE__);
		}
		if(resid == INFINITY) {
			bool writeout = false;
			if(writeout) {
				//create buffer to store ae-field
				int ae_num = in->get_elements();

				ae * ae_tmp = new ae[ae_num];

				//get buffer from device
				cout << "copy buffer to host" << endl;
				exportGaugemomentumBuffer(ae_tmp, in);

				//write out to file
				ofstream out("clmem_p_at_inf");
				if(!out) {
					cout << "Cannot open file.\n";
				}
				for(int i = 0; i < ae_num; i++) {
					out << i << "\t" << ae_tmp[i].e0 << endl;
					out << i << "\t" << ae_tmp[i].e1 << endl;
					out << i << "\t" << ae_tmp[i].e2 << endl;
					out << i << "\t" << ae_tmp[i].e3 << endl;
					out << i << "\t" << ae_tmp[i].e4 << endl;
					out << i << "\t" << ae_tmp[i].e5 << endl;
					out << i << "\t" << ae_tmp[i].e6 << endl;
					out << i << "\t" << ae_tmp[i].e7 << endl;
				}
				out.close();

				//calc sqnorm of ae_tmp
				hmc_float sqnorm = 0.;
				for(int i = 0; i < ae_num; i++) {
					sqnorm += ae_tmp[i].e0 * ae_tmp[i].e0;
					sqnorm += ae_tmp[i].e1 * ae_tmp[i].e1;
					sqnorm += ae_tmp[i].e2 * ae_tmp[i].e2;
					sqnorm += ae_tmp[i].e3 * ae_tmp[i].e3;
					sqnorm += ae_tmp[i].e4 * ae_tmp[i].e4;
					sqnorm += ae_tmp[i].e5 * ae_tmp[i].e5;
					sqnorm += ae_tmp[i].e6 * ae_tmp[i].e6;
					sqnorm += ae_tmp[i].e7 * ae_tmp[i].e7;
				}
				cout << "sqnrom: " << sqnorm << endl;
				free(ae_tmp);
			}
			throw Print_Error_Message("calculation of gaussian gm gave inf! Aborting...", __FILE__, __LINE__);
		}

	}

}

void hardware::code::Gaugemomentum::set_zero_gaugemomentum(const hardware::buffers::Gaugemomentum * buf) const
{
#ifdef CL_VERSION_1_2
	buf->clear();
#else
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(_set_zero_gaugemomentum, &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_force
	int clerr = clSetKernelArg(_set_zero_gaugemomentum, 0, sizeof(cl_mem), buf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(_set_zero_gaugemomentum , gs2, ls2);
#endif
}

void hardware::code::Gaugemomentum::global_squarenorm_reduction(const hardware::buffers::Plain<hmc_float> * out, const hardware::buffers::Plain<hmc_float> * tmp_buf) const
{
	cl_int clerr = clSetKernelArg(gaugemomentum_squarenorm_reduction, 0, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(gaugemomentum_squarenorm_reduction, 1, sizeof(cl_mem), tmp_buf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	cl_uint elems = tmp_buf->get_elements();
	clerr = clSetKernelArg(gaugemomentum_squarenorm_reduction, 2, sizeof(cl_uint), &elems);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(gaugemomentum_squarenorm_reduction, 1, 1);
}

void hardware::code::Gaugemomentum::set_float_to_gaugemomentum_squarenorm_device(const hardware::buffers::Gaugemomentum * clmem_in, const hardware::buffers::Plain<hmc_float> * out) const
{
	//auto spinor_code = get_device()->get_spinor_code();

	//__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gaugemomentum_squarenorm, &ls2, &gs2, &num_groups);

	const hardware::buffers::Plain<hmc_float> clmem_global_squarenorm_buf_glob(num_groups, get_device());

	int clerr = clSetKernelArg(gaugemomentum_squarenorm, 0, sizeof(cl_mem), clmem_in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(gaugemomentum_squarenorm, 1, sizeof(cl_mem), clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(gaugemomentum_squarenorm, 2, sizeof(hmc_float) * ls2, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	get_device()->enqueue_kernel(gaugemomentum_squarenorm, gs2, ls2);

	global_squarenorm_reduction(out, &clmem_global_squarenorm_buf_glob);
}

void hardware::code::Gaugemomentum::saxpy_device(const hardware::buffers::Gaugemomentum * x, const hardware::buffers::Gaugemomentum * y, const hardware::buffers::Plain<hmc_float> * alpha, const hardware::buffers::Gaugemomentum * out) const
{
	//__kernel void gaugemomentum_saxpy(__global ae *x, __global ae *y, hmc_float alpha, __global ae *out){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gaugemomentum_saxpy, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gaugemomentum_saxpy, 0, sizeof(hmc_float), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(gaugemomentum_saxpy, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(gaugemomentum_saxpy, 2, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(gaugemomentum_saxpy, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(gaugemomentum_saxpy , gs2, ls2);
}

void hardware::code::Gaugemomentum::importGaugemomentumBuffer(const hardware::buffers::Gaugemomentum * dest, const ae * const data) const
{
	size_t const REQUIRED_BUFFER_SIZE = get_vol4d(get_device()->getLocalLatticeMemoryExtents()) * NDIM;
	if(dest->get_elements() != REQUIRED_BUFFER_SIZE) {
		throw std::invalid_argument("Destination buffer is not of proper size");
	}

	cl_int clerr;
	if(dest->is_soa()) {
		hardware::buffers::Plain<ae> tmp(REQUIRED_BUFFER_SIZE, dest->get_device());
		tmp.load(data);

		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(gaugemomentum_convert_to_soa, &ls2, &gs2, &num_groups);
		clerr = clSetKernelArg(gaugemomentum_convert_to_soa, 0, sizeof(cl_mem), dest->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gaugemomentum_convert_to_soa, 1, sizeof(cl_mem), tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(gaugemomentum_convert_to_soa, gs2, ls2);
	} else {
		dest->load(data);
	}
}

void hardware::code::Gaugemomentum::exportGaugemomentumBuffer(ae * const dest, const hardware::buffers::Gaugemomentum * buf) const
{
	size_t const REQUIRED_BUFFER_SIZE = get_vol4d(get_device()->getLocalLatticeMemoryExtents()) * NDIM;
	if(buf->get_elements() != REQUIRED_BUFFER_SIZE) {
		throw std::invalid_argument("Source buffer is not of proper size");
	}

	cl_int clerr;
	if(buf->is_soa()) {
		hardware::buffers::Plain<ae> tmp(REQUIRED_BUFFER_SIZE, buf->get_device());

		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(gaugemomentum_convert_from_soa, &ls2, &gs2, &num_groups);
		clerr = clSetKernelArg(gaugemomentum_convert_from_soa, 0, sizeof(cl_mem), tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gaugemomentum_convert_from_soa, 1, sizeof(cl_mem), buf->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(gaugemomentum_convert_from_soa, gs2, ls2);

		tmp.dump(dest);
	} else {
		buf->dump(dest);
	}
}


ClSourcePackage hardware::code::Gaugemomentum::get_sources() const noexcept
{
	return basic_gaugemomentum_code;
}
