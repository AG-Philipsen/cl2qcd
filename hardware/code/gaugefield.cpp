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

#include "gaugefield.hpp"

#include <cmath>

#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include "../buffers/3x3.hpp"
#include "flopUtilities.hpp"
#include "../../geometry/latticeGrid.hpp"

using namespace std;

void hardware::code::Gaugefield::fill_kernels()
{
	basic_opencl_code = get_basic_sources() << "operations_geometry.cl" << "operations_complex.h"
	                    << "operations_matrix_su3.cl" << "operations_matrix.cl" << "operations_gaugefield.cl";
	
	logger.debug() << "Creating Gaugefield kernels...";
	
	plaquette = createKernel("plaquette") << basic_opencl_code << "gaugeobservables_plaquette.cl";
	plaquette_reduction = createKernel("plaquette_reduction") << basic_opencl_code << "gaugeobservables_plaquette.cl";
	if(kernelParameters->getUseRectangles() == true) {
		rectangles = createKernel("rectangles") << basic_opencl_code << "gaugeobservables_rectangles.cl";
		rectangles_reduction = createKernel("rectangles_reduction") << basic_opencl_code << "gaugeobservables_rectangles.cl";
	}
	if(get_device()->getGridSize().tExtent == 1) { //todo: work over this! should be something like "parallelized direction"
		polyakov = createKernel("polyakov") << basic_opencl_code << "gaugeobservables_polyakov.cl";
		polyakov_md_local = nullptr;
		polyakov_md_merge = nullptr;
	} else {
		polyakov = nullptr;
		polyakov_md_local = createKernel("polyakov_md_local") << basic_opencl_code << "gaugeobservables_polyakov.cl";
		polyakov_md_merge = createKernel("polyakov_md_merge") << basic_opencl_code << "gaugeobservables_polyakov.cl";
	}
	polyakov_reduction = createKernel("polyakov_reduction") << basic_opencl_code << "gaugeobservables_polyakov.cl";
	if(kernelParameters->getUseSmearing() == true) {
		stout_smear = createKernel("stout_smear") << basic_opencl_code << "operations_gaugemomentum.cl" << "stout_smear.cl";
	}
	convertGaugefieldToSOA = createKernel("convertGaugefieldToSOA") << basic_opencl_code << "gaugefield_convert.cl";
	convertGaugefieldFromSOA = createKernel("convertGaugefieldFromSOA") << basic_opencl_code << "gaugefield_convert.cl";
}

void hardware::code::Gaugefield::clear_kernels()
{
	cl_int clerr = CL_SUCCESS;
	
	logger.debug() << "Clearing Gaugefield kernels...";

	clerr = clReleaseKernel(plaquette);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(polyakov) {
		clerr = clReleaseKernel(polyakov);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(polyakov_md_local) {
		clerr = clReleaseKernel(polyakov_md_local);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(polyakov_md_merge) {
		clerr = clReleaseKernel(polyakov_md_merge);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	clerr = clReleaseKernel(plaquette_reduction);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(polyakov_reduction);
	if(kernelParameters->getUseRectangles() == true) {
		clerr = clReleaseKernel(rectangles);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(rectangles_reduction);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(kernelParameters->getUseSmearing() == true) {
		clerr = clReleaseKernel(stout_smear);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	clerr = clReleaseKernel(convertGaugefieldToSOA);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(convertGaugefieldFromSOA);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
}

void hardware::code::Gaugefield::plaquette_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> * plaq, const hardware::buffers::Plain<hmc_float> * tplaq, const hardware::buffers::Plain<hmc_float> * splaq) const
{
	using namespace hardware::buffers;

	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(plaquette, &ls, &gs, &num_groups);

	const Plain<hmc_float> clmem_plaq_buf_glob(num_groups, get_device());
	const Plain<hmc_float> clmem_tplaq_buf_glob(num_groups, get_device());
	const Plain<hmc_float> clmem_splaq_buf_glob(num_groups, get_device());

	int buf_loc_size_float = sizeof(hmc_float) * ls;

	//set arguments
	// run local plaquette calculation and first part of reduction
	int clerr = clSetKernelArg(plaquette, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 1, sizeof(cl_mem), clmem_plaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 2, sizeof(cl_mem), clmem_tplaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 3, sizeof(cl_mem), clmem_splaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 4, buf_loc_size_float, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 5, buf_loc_size_float, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 6, buf_loc_size_float, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(plaquette, gs, ls);

	// run second part of plaquette reduction

	this->get_work_sizes(plaquette_reduction, &ls, &gs, &num_groups);

	clerr = clSetKernelArg(plaquette_reduction, 0, sizeof(cl_mem), clmem_plaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 1, sizeof(cl_mem), clmem_tplaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 2, sizeof(cl_mem), clmem_splaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 3, sizeof(cl_mem), plaq->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 4, sizeof(cl_mem), tplaq->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 5, sizeof(cl_mem), splaq->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 6, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);


	///@todo improve
	ls = 1;
	gs = 1;
	get_device()->enqueue_kernel(plaquette_reduction, gs, ls);
}

void hardware::code::Gaugefield::rectangles_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> * rect) const
{
	if(rectangles == nullptr) {
		throw std::logic_error("Rectangles are not enabled.");
	}

	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(rectangles, &ls, &gs, &num_groups);

	const hardware::buffers::Plain<hmc_float> clmem_rect_buf_glob(num_groups, get_device());

	int buf_loc_size_float = sizeof(hmc_float) * ls;

	//set arguments
	// run local rectangles calculation and first part of reduction
	int clerr = clSetKernelArg(rectangles, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(rectangles, 1, sizeof(cl_mem), clmem_rect_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(rectangles, 2, buf_loc_size_float, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(rectangles, gs, ls);

	// run second part of rectangles reduction

	this->get_work_sizes(rectangles_reduction, &ls, &gs, &num_groups);

	clerr = clSetKernelArg(rectangles_reduction, 0, sizeof(cl_mem), clmem_rect_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(rectangles_reduction, 1, sizeof(cl_mem), rect->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(rectangles_reduction, 2, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);


	///@todo improve
	ls = 1;
	gs = 1;
	get_device()->enqueue_kernel(rectangles_reduction, gs, ls);

}

void hardware::code::Gaugefield::polyakov_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_complex> * pol) const
{
	if(!polyakov) {
		throw std::logic_error("This function cannot be called in multi-device environments.");
	}
	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups_pol;
	this->get_work_sizes(polyakov, &ls, &gs, &num_groups_pol);
	int buf_loc_size_complex = sizeof(hmc_complex) * ls;

	const hardware::buffers::Plain<hmc_complex> clmem_polyakov_buf_glob(num_groups_pol, get_device());

	// local polyakov compuation and first part of reduction
	int clerr = clSetKernelArg(polyakov, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov, 1, sizeof(cl_mem), clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov, 2, buf_loc_size_complex, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(polyakov, gs, ls);

	// second part of polyakov reduction
	cl_uint num_groups_reduce;
	this->get_work_sizes(polyakov_reduction, &ls, &gs, &num_groups_reduce);

	clerr = clSetKernelArg(polyakov_reduction, 0, sizeof(cl_mem), clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov_reduction, 1, sizeof(cl_mem), pol->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov_reduction, 2, sizeof(cl_uint), &num_groups_pol);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	///@todo improve
	ls = 1;
	gs = 1;
	get_device()->enqueue_kernel(polyakov_reduction, gs, ls);

}

void hardware::code::Gaugefield::polyakov_md_local_device(const hardware::buffers::Plain<Matrixsu3> * partial_results, const hardware::buffers::SU3* gf) const
{
	if(!polyakov_md_local) {
		throw std::logic_error("This function can only called in multi-device environments.");
	}

	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(polyakov_md_local, &ls, &gs, &num_groups);

	// local polyakov compuation and first part of reduction
	int clerr = clSetKernelArg(polyakov_md_local, 0, sizeof(cl_mem), partial_results->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov_md_local, 1, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(polyakov_md_local, gs, ls);
}

void hardware::code::Gaugefield::polyakov_md_merge_device(const hardware::buffers::Plain<Matrixsu3> * partial_results, const cl_uint num_slices, const hardware::buffers::Plain<hmc_complex> * pol) const
{
	if(!polyakov_md_merge) {
		throw std::logic_error("This function can only called in multi-device environments.");
	}
	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups_merge;
	this->get_work_sizes(polyakov_md_merge, &ls, &gs, &num_groups_merge);
	int buf_loc_size_complex = sizeof(hmc_complex) * ls;

	const hardware::buffers::Plain<hmc_complex> clmem_polyakov_buf_glob(num_groups_merge, get_device());

	// local polyakov compuation and first part of reduction
	int clerr = clSetKernelArg(polyakov_md_merge, 0, sizeof(cl_mem), clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov_md_merge, 1, sizeof(cl_mem), partial_results->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov_md_merge, 2, sizeof(cl_uint), &num_slices);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov_md_merge, 3, buf_loc_size_complex, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(polyakov_md_merge, gs, ls);

	// second part of polyakov reduction
	cl_uint num_groups_reduce;
	this->get_work_sizes(polyakov_reduction, &ls, &gs, &num_groups_reduce);

	clerr = clSetKernelArg(polyakov_reduction, 0, sizeof(cl_mem), clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov_reduction, 1, sizeof(cl_mem), pol->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov_reduction, 2, sizeof(cl_uint), &num_groups_merge);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	///@todo improve
	ls = 1;
	gs = 1;
	get_device()->enqueue_kernel(polyakov_reduction, gs, ls);

	get_device()->synchronize();
}

void hardware::code::Gaugefield::gaugeobservables_rectangles(const hardware::buffers::SU3 * gf, hmc_float * rect_out) const
{
	const hardware::buffers::Plain<hmc_float> rect(1, get_device());

	//measure plaquette
	rectangles_device(gf, &rect);

	//read out values
	//NOTE: these are blocking calls!
	rect.dump(rect_out);
	//NOTE: the rectangle value has not been normalized since it is mostly used for the HMC where one needs the absolute value
}

void hardware::code::Gaugefield::stout_smear_device(const hardware::buffers::SU3 * in, const hardware::buffers::SU3 * out) const
{
	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(stout_smear, &ls, &gs, &num_groups);

	int clerr = clSetKernelArg(stout_smear, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(stout_smear, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(stout_smear , gs, ls);

	return;
}

void hardware::code::Gaugefield::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);
}

size_t hardware::code::Gaugefield::get_read_write_size(const std::string& in) const
{
	//Depending on the compile-options, one has different sizes...
	size_t D = kernelParameters->getFloatSize();
	size_t R = kernelParameters->getMatSize();
	//factor for complex numbers
	int C = 2;
	const size_t VOL4D = kernelParameters->getLatticeVolume();
	if (in == "polyakov") {
		//this kernel reads NTIME*VOLSPACE=VOL4D su3matrices and writes NUM_GROUPS complex numbers
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(polyakov, &ls2, &gs2, &num_groups);
		return VOL4D * D * R + num_groups * C * D;
	}
	if (in == "polyakov_reduction") {
		//this kernel reads NUM_GROUPS complex numbers and writes 1 complex number
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(polyakov_reduction, &ls2, &gs2, &num_groups);
		return (num_groups + 1 ) * C * D;
	}
	if (in == "plaquette") {
		//this kernel reads in VOL4D * ND * (ND-1) su3matrices and writes 3*num_groups real numbers
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(plaquette, &ls2, &gs2, &num_groups);
		return C * 4 * 3 * VOL4D * D * R + 3 * D * num_groups;
	}
	if (in == "plaquette_reduction") {
		//this kernel reads 3*NUM_GROUPS real numbers and writes 3 real numbers
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(polyakov_reduction, &ls2, &gs2, &num_groups);
		return (num_groups + 1 ) * 3 * C * D;
	}
	if (in == "rectangles") {
		return module_metric_not_implemented<size_t>();
	}
	if (in == "rectangles_reduction") {
		return module_metric_not_implemented<size_t>();
	}
	if (in == "stout_smear") {
		//this kernel reads in a complete gaugefield + a staple on each site and writes out a complete gaugefield
		return VOL4D * NDIM * D * R * (6 * (NDIM - 1) + 1 + 1 );
	}
	if(in == "convertGaugefieldToSOA") {
		return 2 * kernelParameters->getLatticeVolume() * NDIM * R * C * D;
	}
	if(in == "convertGaugefieldFromSOA") {
		return 2 * kernelParameters->getLatticeVolume() * NDIM * R * C * D;
	}
	return 0;
}

uint64_t hardware::code::Gaugefield::get_flop_size(const std::string& in) const
{
	const size_t VOL4D = kernelParameters->getLatticeVolume();
	const size_t VOLSPACE = kernelParameters->getSpatialLatticeVolume();
	if (in == "polyakov") {
		//this kernel performs NTIME -1 su3matrix-multiplications, takes a complex trace and adds these real values over VOLSPACE
		return VOLSPACE * ( (kernelParameters->getNt() - 1) * getFlopSu3MatrixTimesSu3Matrix() + getFlopSu3MatrixTrace()) ;
	}
	if (in == "polyakov_reduction") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "plaquette") {
		//this kernel performs 3 su3matrix-mutliplications, a real su3 trace and sums over VOL4D and mu and nu (nu<mu)
		return VOL4D * NDIM * (NDIM - 1) * ( 3 + NC);
	}
	if (in == "plaquette_reduction") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "rectangles") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "rectangles_reduction") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "stout_smear") {
		return module_metric_not_implemented<uint64_t>();
	}
	return 0;
}

void hardware::code::Gaugefield::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, polyakov);
	Opencl_Module::print_profiling(filename, polyakov_reduction);
	Opencl_Module::print_profiling(filename, plaquette);
	Opencl_Module::print_profiling(filename, rectangles);
	Opencl_Module::print_profiling(filename, plaquette_reduction);
	Opencl_Module::print_profiling(filename, stout_smear);
	Opencl_Module::print_profiling(filename, convertGaugefieldToSOA);
	Opencl_Module::print_profiling(filename, convertGaugefieldFromSOA);
}

void hardware::code::Gaugefield::importGaugefield(const hardware::buffers::SU3 * gaugefield, const Matrixsu3 * const data) const
{
	using namespace hardware::buffers;

	logger.trace() << "Import gaugefield to get_device()";
	if(get_device()->get_prefers_soa()) {
		Plain<Matrixsu3> tmp(gaugefield->get_elements(), get_device());
		tmp.load(data);
		convertGaugefieldToSOA_device(gaugefield, &tmp);
	} else {
		gaugefield->load(data);
	}
}

void hardware::code::Gaugefield::exportGaugefield(Matrixsu3 * const dest, const hardware::buffers::SU3 * gaugefield) const
{
	using namespace hardware::buffers;

	logger.trace() << "Exporting gaugefield from get_device()";
	if(get_device()->get_prefers_soa()) {
		Plain<Matrixsu3> tmp(gaugefield->get_elements(), get_device());
		convertGaugefieldFromSOA_device(&tmp, gaugefield);
		tmp.dump(dest);
	} else {
		gaugefield->dump(dest);
	}
}

void hardware::code::Gaugefield::convertGaugefieldToSOA_device(const hardware::buffers::SU3 * out, const hardware::buffers::Plain<Matrixsu3> * in) const
{
	if(!out->is_soa()) {
		throw std::invalid_argument("Destination buffer must be a SOA buffer");
	}

	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convertGaugefieldToSOA, &ls2, &gs2, &num_groups);

	//set arguments
	int clerr = clSetKernelArg(convertGaugefieldToSOA, 0, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convertGaugefieldToSOA, 1, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(convertGaugefieldToSOA, gs2, ls2);
}

void hardware::code::Gaugefield::convertGaugefieldFromSOA_device(const hardware::buffers::Plain<Matrixsu3> * out, const hardware::buffers::SU3 * in) const
{
	if(!in->is_soa()) {
		throw std::invalid_argument("Source buffer must be a SOA buffer");
	}

	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convertGaugefieldFromSOA, &ls2, &gs2, &num_groups);

	//set arguments
	int clerr = clSetKernelArg(convertGaugefieldFromSOA, 0, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convertGaugefieldFromSOA, 1, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(convertGaugefieldFromSOA, gs2, ls2);
}

hardware::code::Gaugefield::Gaugefield( const hardware::code::OpenClKernelParametersInterface& kernelParameters , const hardware::Device * device)
	: Opencl_Module(kernelParameters, device),
	  stout_smear(0), rectangles(0), rectangles_reduction(0)
{
	fill_kernels();
}


hardware::code::Gaugefield::~Gaugefield()
{
	clear_kernels();
}

ClSourcePackage hardware::code::Gaugefield::get_sources() const noexcept
{
	return basic_opencl_code;
}
