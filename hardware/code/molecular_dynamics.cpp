/** @file
 * Implementation of the hardware::code::Molecular_Dynamics class
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

#include "molecular_dynamics.hpp"

#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include "fermions.hpp"
#include "prng.hpp"
#include "gaugemomentum.hpp"
#include "spinors.hpp"
#include <cmath>
#include "flopUtilities.hpp"

using namespace std;

static bool use_multipass_gauge_force_tlsym(const hardware::Device * device)
{
	return (device->get_name() == "Tahiti");
}

void hardware::code::Molecular_Dynamics::fill_kernels()
{
	basic_molecular_dynamics_code = get_basic_sources() << "operations_geometry.cl" << "operations_complex.h" << "types_fermions.h" << "types_hmc.h" << "operations_matrix_su3.cl" << "operations_matrix.cl" << "operations_gaugefield.cl" << "operations_su3vec.cl" << "operations_spinor.cl" << "spinorfield.cl" << "operations_gaugemomentum.cl";
	
	ClSourcePackage prng_code = get_device()->getPrngCode()->get_sources();
	
	logger.debug() << "Creating Molecular_Dynamics kernels...";

	//init kernels for HMC
	if(kernelParameters->getUseEo() == true) {
		fermion_force_eo_0 = createKernel("fermion_force_eo_0") << basic_molecular_dynamics_code  << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "force_fermion_eo.cl";
		fermion_force_eo_1 = createKernel("fermion_force_eo_1") << basic_molecular_dynamics_code  << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "force_fermion_eo.cl";
		fermion_force_eo_2 = createKernel("fermion_force_eo_2") << basic_molecular_dynamics_code  << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "force_fermion_eo.cl";
		fermion_force_eo_3 = createKernel("fermion_force_eo_3") << basic_molecular_dynamics_code  << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "force_fermion_eo.cl";
		fermion_stagg_partial_force_eo = createKernel("fermion_staggered_partial_force_eo") << basic_molecular_dynamics_code << "operations_staggered.cl" << "spinorfield_staggered_eo.cl" << "force_staggered_fermion_eo.cl";
	}
	fermion_force = createKernel("fermion_force") << basic_molecular_dynamics_code << "fermionmatrix.cl" << "force_fermion.cl";
	md_update_gaugefield = createKernel("md_update_gaugefield") << basic_molecular_dynamics_code << "md_update_gaugefield.cl";
	gauge_force = createKernel("gauge_force") << basic_molecular_dynamics_code  << "force_gauge.cl";
	if(kernelParameters->getUseRectangles() == true) {
		//at the time of writing this kernel, the OpenCL compiler crashed the kernel using optimizations
		if(gauge_force_tlsym_tmp) {
			gauge_force_tlsym = 0;
			gauge_force_tlsym_1 = createKernel("gauge_force_tlsym_multipass1_tpe") << basic_molecular_dynamics_code << "force_gauge_tlsym.cl";
			gauge_force_tlsym_2 = createKernel("gauge_force_tlsym_multipass2_tpe") << basic_molecular_dynamics_code << "force_gauge_tlsym.cl";
			gauge_force_tlsym_3 = createKernel("gauge_force_tlsym_multipass3_tpe") << basic_molecular_dynamics_code << "force_gauge_tlsym.cl";
			gauge_force_tlsym_4 = createKernel("gauge_force_tlsym_multipass4_tpe") << basic_molecular_dynamics_code << "force_gauge_tlsym.cl";
			gauge_force_tlsym_5 = createKernel("gauge_force_tlsym_multipass5_tpe") << basic_molecular_dynamics_code << "force_gauge_tlsym.cl";
			gauge_force_tlsym_6 = createKernel("gauge_force_tlsym_multipass6_tpe") << basic_molecular_dynamics_code << "force_gauge_tlsym.cl";
		} else {
			gauge_force_tlsym = createKernel("gauge_force_tlsym") << basic_molecular_dynamics_code << "force_gauge_tlsym.cl";
			gauge_force_tlsym_1 = 0;
			gauge_force_tlsym_2 = 0;
			gauge_force_tlsym_3 = 0;
			gauge_force_tlsym_4 = 0;
			gauge_force_tlsym_5 = 0;
			gauge_force_tlsym_6 = 0;
		}
	}
	if(kernelParameters->getUseSmearing() == true) {
		stout_smear_fermion_force = createKernel("stout_smear_fermion_force") << basic_molecular_dynamics_code << "force_fermion_stout_smear.cl";
	}
}

void hardware::code::Molecular_Dynamics::clear_kernels()
{
	cl_uint clerr = CL_SUCCESS;

	logger.debug() << "Clearing Molecular_Dynamics kernels.." ;
	
	if(kernelParameters->getUseEo() == true) {
		clerr = clReleaseKernel(fermion_force_eo_3);
		clerr = clReleaseKernel(fermion_force_eo_2);
		clerr = clReleaseKernel(fermion_force_eo_1);
		clerr = clReleaseKernel(fermion_force_eo_0);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(fermion_stagg_partial_force_eo);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	} else {
		clerr = clReleaseKernel(fermion_force);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	clerr = clReleaseKernel(md_update_gaugefield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(gauge_force);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(gauge_force_tlsym) {
		clerr = clReleaseKernel(gauge_force_tlsym);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(gauge_force_tlsym_1) {
		clerr = clReleaseKernel(gauge_force_tlsym_1);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(gauge_force_tlsym_2) {
		clerr = clReleaseKernel(gauge_force_tlsym_2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(gauge_force_tlsym_3) {
		clerr = clReleaseKernel(gauge_force_tlsym_3);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(gauge_force_tlsym_4) {
		clerr = clReleaseKernel(gauge_force_tlsym_4);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(gauge_force_tlsym_5) {
		clerr = clReleaseKernel(gauge_force_tlsym_5);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(gauge_force_tlsym_6) {
		clerr = clReleaseKernel(gauge_force_tlsym_6);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(kernelParameters->getUseSmearing() == true) {
		clerr = clReleaseKernel(stout_smear_fermion_force);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
}

void hardware::code::Molecular_Dynamics::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);
	return;
}


size_t hardware::code::Molecular_Dynamics::get_read_write_size(const std::string& in) const
{
//Depending on the compile-options, one has different sizes...
	size_t D = kernelParameters->getFloatSize();
	//this returns the number of entries in an su3-matrix
	size_t R = kernelParameters->getMatSize();
	//this is the number of spinors in the system (or number of sites)
	size_t S = kernelParameters->getSpinorFieldSize();
	size_t Seo = kernelParameters->getEoprecSpinorFieldSize();
	//this is the number of links in the system (and of gaugemomenta)
	size_t G = kernelParameters->getLatticeVolume() * NDIM;
	//factor for complex numbers
	int C = 2;
	//this is the same as in the function above
	//NOTE: 1 spinor has NC*NDIM = 12 complex entries
	//NOTE: 1 ae has NC*NC-1 = 8 real entries
	int A = getSu3AlgebraSize();
	if (in == "md_update_gaugefield") {
		//this kernel reads 1 ae and 1 su3 matrix and writes 1 su3 matrix for every link
		return (A + (1 + 1) * R * C) * D * G;
	}
	if (in == "gauge_force") {
		//this kernel reads ingredients for 1 staple plus 1 su3matrix and writes 1 ae for every link
		return G * D * (R * C * ( 6 * (NDIM - 1) + 1 ) + A );
	}
	if (in == "gauge_force_tlsym") {
		//this kernel reads ingredients for 1 rect-staple plus 1 su3matrix and writes 1 ae for every link
		//the rect staple is the same as the normal staple, but with 3 add. contributions and 2 add. matrices in each contribution, so instead of 2 * 3 = 6, one reads 6 * 5 matrices in each direction
		return G * D * (R * C * ( 6 * 5 * (NDIM - 1 ) + 1 ) + A );
	}
	if (in == "fermion_force") {
		//this kernel reads 16 spinors, 8 su3matrices and writes 8 ae per site
		//NOTE: the kernel now runs over all ae instead of all sites, but this must be equivalent!
		return (C * 12 * (16) + C * 8 * R + 8 * A) * D * S;
	}
	if (in == "fermion_force_eo") {
		//this kernel reads 16 spinors, 8 su3matrices and writes 8 ae per site_eo
		return (C * 12 * (16) + C * 8 * R + 8 * A) * D * Seo;
	}
	if (in == "fermion_staggered_partial_force_eo") {
		//this kernel reads 8 su3vec, 4 su3matrices and writes 4 ae per site_eo
		return (C * 3 * (8) + C * 4 * R + 4 * A) * D * Seo;
	}
	if (in == "stout_smear_fermion_force") {
		return module_metric_not_implemented<uint64_t>();
	}
	return 0;
}

uint64_t hardware::code::Molecular_Dynamics::get_flop_size(const std::string& in) const
{
	//this is the number of spinors in the system (or number of sites)
	size_t Seo = kernelParameters->getEoprecSpinorFieldSize();
	//this is the number of links in the system (and of gaugemomenta)
	uint64_t G = kernelParameters->getLatticeVolume() * NDIM;
	//NOTE: 1 ae has NC*NC-1 = 8 real entries
	uint64_t A = getSu3AlgebraSize();
	//this returns the number of entries in an su3-matrix
	uint64_t R = kernelParameters->getMatSize();
	//this is the same as in the function above
	if (in == "md_update_gaugefield") {
		//this kernel performs one exp(i ae) ( = 327 flops + 1 su3 mult ) and 1 su3 mult per link
		return (getFlopSu3MatrixTimesSu3Matrix() * ( 1 + 1)  + 327 ) * G;
	}
	if (in == "gauge_force") {
		//this kernel calculates 1 staple (= 4*ND-1 su3_su3 + 2*ND-1 su3_add), 1 su3*su3, 1 tr_lambda_u (19 flops) plus 8 add and 8 mult per ae
		return ( 4 * (NDIM - 1) * getFlopSu3MatrixTimesSu3Matrix() + 2 * (NDIM - 1) * 18 + 1 * getFlopSu3MatrixTimesSu3Matrix() + 19  + A * ( 1 + 1 )
		       ) * G;
	}
	if (in == "gauge_force_tlsym") {
		//this kernel calculates 1 rect-staple (= 24*ND-1 su3_su3 + 6*ND-1 su3_add), 1 su3*su3, 1 tr_lambda_u (19 flops) plus 8 add and 8 mult per ae
		//24 = 6 contr. per dir, 4 mat_mat per contr.
		return ( 24 * (NDIM - 1) * getFlopSu3MatrixTimesSu3Matrix() + 6 * (NDIM - 1) * 18 + 1 * getFlopSu3MatrixTimesSu3Matrix() + 19  + A * ( 1 + 1 )
		       ) * G;
	}
	if (in == "fermion_force") {
		//this kernel performs NDIM * ( 4 * su3vec_acc (6 flops) + tr(v*u) (126 flops) + tr_lambda_u(19 flops) + update_ae(8*2 flops) + su3*su3 + su3*complex (flop_complex_mult * R ) ) per site
		//NOTE: the kernel now runs over all ae instead of all sites, but this must be equivalent!
		return Seo * NDIM * ( 4 * 6 + 126 + 19 + 8 * 2 + getFlopSu3MatrixTimesSu3Matrix() + getFlopComplexMult() * R );
	}
	if (in == "fermion_force_eo") {
		//this kernel performs NDIM * ( 4 * su3vec_acc (6 flops) + tr(v*u) (126 flops) + tr_lambda_u(19 flops) + update_ae(8*2 flops) + su3*su3 + su3*complex (flop_complex_mult * R ) ) per site_eo
		return Seo * NDIM * ( 4 * 6 + 126 + 19 + 8 * 2 + getFlopSu3MatrixTimesSu3Matrix() + getFlopComplexMult() * R );
	}
	if (in == "fermion_staggered_partial_force_eo") {
		//this kernel performs (su3vec_times_real               (6 flops) +
		//                      u_times_v_dagger                (get_flop_su3vec_direct_su3vec()) +
		//                      multiply_matrix_3x3             (get_flop_su3_su3()) +
		//                      traceless_antihermitian_part    (18 flops) +
		//                      multiply_matrix3x3_by_complex   (R * get_flop_complex_mult()) +
		//                      build_ae_from_su3               (9 flops) +
		//                      update_ae                       (8*2 flops))  per direction per site_eo
		return Seo * NDIM * (6 + getFlopSu3VecDirectSu3Vec() + getFlopSu3MatrixTimesSu3Matrix() +
		                     18 + R * 18 + getFlopComplexMult() + 9 + 16);
	}
	if (in == "stout_smear_fermion_force") {
		return module_metric_not_implemented<uint64_t>();
	}
	return 0;
}

void hardware::code::Molecular_Dynamics::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, md_update_gaugefield);
	Opencl_Module::print_profiling(filename, gauge_force);
	Opencl_Module::print_profiling(filename, gauge_force_tlsym);
	Opencl_Module::print_profiling(filename, gauge_force_tlsym_1);
	Opencl_Module::print_profiling(filename, gauge_force_tlsym_2);
	Opencl_Module::print_profiling(filename, gauge_force_tlsym_3);
	Opencl_Module::print_profiling(filename, gauge_force_tlsym_4);
	Opencl_Module::print_profiling(filename, gauge_force_tlsym_5);
	Opencl_Module::print_profiling(filename, gauge_force_tlsym_6);
	Opencl_Module::print_profiling(filename, fermion_force);
	Opencl_Module::print_profiling(filename, fermion_force_eo_0);
	Opencl_Module::print_profiling(filename, fermion_force_eo_1);
	Opencl_Module::print_profiling(filename, fermion_force_eo_2);
	Opencl_Module::print_profiling(filename, fermion_force_eo_3);
	Opencl_Module::print_profiling(filename, fermion_stagg_partial_force_eo);
	Opencl_Module::print_profiling(filename, stout_smear_fermion_force);
}

void hardware::code::Molecular_Dynamics::md_update_gaugefield_device(const hardware::buffers::Gaugemomentum * gm_in, const hardware::buffers::SU3 * gf_out, hmc_float eps) const
{
	// __kernel void md_update_gaugefield(hmc_float eps, __global ae * p_in, __global ocl_s_gaugefield * u_inout){
	hmc_float tmp = eps;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(md_update_gaugefield, &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_force
	int clerr = clSetKernelArg(md_update_gaugefield, 0, sizeof(hmc_float), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(md_update_gaugefield, 1, sizeof(cl_mem), gm_in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(md_update_gaugefield, 2, sizeof(cl_mem), gf_out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( md_update_gaugefield , gs2, ls2);
}

void hardware::code::Molecular_Dynamics::gauge_force_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Gaugemomentum * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gauge_force, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gauge_force, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(gauge_force, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( gauge_force , gs2, ls2);

	if(logger.beDebug()) {
		auto gm_code = get_device()->getGaugemomentumCode();
		hardware::buffers::Plain<hmc_float> gauge_force_tmp(1, get_device());
		hmc_float gauge_force_energy = 0.;
		gm_code->set_float_to_gaugemomentum_squarenorm_device(out, &gauge_force_tmp);
		gauge_force_tmp.dump(&gauge_force_energy);

		logger.debug() <<  "\t\t\tFORCE [GAUGE]:\t" << gauge_force_energy;

		if(gauge_force_energy != gauge_force_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

}

void hardware::code::Molecular_Dynamics::gauge_force_tlsym_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Gaugemomentum * out) const
{
	if(gauge_force_tlsym_tmp) {
		// run multipass
		size_t foo, ls;
		cl_uint bla;
		size_t global_size = gauge_force_tlsym_tmp->get_elements();

		this->get_work_sizes(gauge_force_tlsym_1, &ls, &foo, &bla);
		int clerr = clSetKernelArg(gauge_force_tlsym_1, 0, sizeof(cl_mem), gf->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gauge_force_tlsym_1, 1, sizeof(cl_mem), gauge_force_tlsym_tmp->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(gauge_force_tlsym_1, global_size, ls);

		this->get_work_sizes(gauge_force_tlsym_2, &ls, &foo, &bla);
		clerr = clSetKernelArg(gauge_force_tlsym_2, 0, sizeof(cl_mem), gf->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gauge_force_tlsym_2, 1, sizeof(cl_mem), gauge_force_tlsym_tmp->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(gauge_force_tlsym_2, global_size, ls);

		this->get_work_sizes(gauge_force_tlsym_3, &ls, &foo, &bla);
		clerr = clSetKernelArg(gauge_force_tlsym_3, 0, sizeof(cl_mem), gf->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gauge_force_tlsym_3, 1, sizeof(cl_mem), gauge_force_tlsym_tmp->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(gauge_force_tlsym_3, global_size, ls);

		this->get_work_sizes(gauge_force_tlsym_4, &ls, &foo, &bla);
		clerr = clSetKernelArg(gauge_force_tlsym_4, 0, sizeof(cl_mem), gf->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gauge_force_tlsym_4, 1, sizeof(cl_mem), gauge_force_tlsym_tmp->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(gauge_force_tlsym_4, global_size, ls);

		this->get_work_sizes(gauge_force_tlsym_5, &ls, &foo, &bla);
		clerr = clSetKernelArg(gauge_force_tlsym_5, 0, sizeof(cl_mem), gf->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gauge_force_tlsym_5, 1, sizeof(cl_mem), gauge_force_tlsym_tmp->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(gauge_force_tlsym_5, global_size, ls);

		this->get_work_sizes(gauge_force_tlsym_6, &ls, &foo, &bla);
		clerr = clSetKernelArg(gauge_force_tlsym_6, 0, sizeof(cl_mem), gf->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gauge_force_tlsym_6, 1, sizeof(cl_mem), out->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gauge_force_tlsym_6, 2, sizeof(cl_mem), gauge_force_tlsym_tmp->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(gauge_force_tlsym_6, global_size, ls);

	} else {
		//query work-sizes for kernel
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(gauge_force_tlsym, &ls2, &gs2, &num_groups);
		//set arguments
		int clerr = clSetKernelArg(gauge_force_tlsym, 0, sizeof(cl_mem), gf->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		clerr = clSetKernelArg(gauge_force_tlsym, 1, sizeof(cl_mem), out->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		get_device()->enqueue_kernel( gauge_force_tlsym , gs2, ls2);
	}

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> gauge_force_tlsym_tmp(1, get_device());
		auto gm_code = get_device()->getGaugemomentumCode();
		hmc_float gauge_force_tlsym_energy = 0.;
		gm_code->set_float_to_gaugemomentum_squarenorm_device(out, &gauge_force_tlsym_tmp);
		gauge_force_tlsym_tmp.dump(&gauge_force_tlsym_energy);

		logger.debug() <<  "\t\t\tFORCE [TLSYM]:\t" << gauge_force_tlsym_energy;

		if(gauge_force_tlsym_energy != gauge_force_tlsym_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}

void hardware::code::Molecular_Dynamics::fermion_force_device(const hardware::buffers::Plain<spinor> * Y, const hardware::buffers::Plain<spinor> * X, const hardware::buffers::SU3 * gf, const hardware::buffers::Gaugemomentum * out, hmc_float kappa) const
{
	using namespace hardware::buffers;

	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = kernelParameters->getKappa();
	else kappa_tmp = kappa;

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(fermion_force, &ls2, &gs2, &num_groups);
	//set arguments
	//fermion_force(field, Y, X, out);
	int clerr = clSetKernelArg(fermion_force, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 1, sizeof(cl_mem), Y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 2, sizeof(cl_mem), X->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_force, 4, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( fermion_force , gs2, ls2);

	if(logger.beDebug()) {
		Plain<hmc_float> noneo_force_tmp(1, get_device());
		auto gm_code = get_device()->getGaugemomentumCode();
		hmc_float noneo_force_energy = 0.;
		gm_code->set_float_to_gaugemomentum_squarenorm_device(out, &noneo_force_tmp);
		noneo_force_tmp.dump(&noneo_force_energy);
		logger.debug() <<  "\t\t\tFORCE [DET] (noneo):\t" << noneo_force_energy;

		if(noneo_force_energy != noneo_force_energy) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}

//the argument kappa is set to ARG_DEF as default
void hardware::code::Molecular_Dynamics::fermion_force_eo_device(const hardware::buffers::Spinor * Y, const hardware::buffers::Spinor * X, const hardware::buffers::SU3 * gf, const hardware::buffers::Gaugemomentum * out, int evenodd, hmc_float kappa) const
{
	using namespace hardware::buffers;

	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = kernelParameters->getKappa();
	else kappa_tmp = kappa;

	if(logger.beDebug()) {
		Plain<hmc_float> tmp(1, get_device());
		auto spinor_code = get_device()->getSpinorCode();
		hmc_float resid;
		spinor_code->set_float_to_global_squarenorm_eoprec_device(Y, &tmp);
		tmp.dump(&resid);
		logger.debug() <<  "\t\t\tY:\t" << resid;
		spinor_code->set_float_to_global_squarenorm_eoprec_device(X, &tmp);
		tmp.dump(&resid);
		logger.debug() <<  "\t\t\tX:\t" << resid;

		auto gf_code = get_device()->getGaugefieldCode();
		Plain<hmc_complex> pol_dev(1, get_device());
		gf_code->polyakov_device(gf, &pol_dev);
		hmc_complex pol_host;
		pol_dev.dump(&pol_host);
		logger.debug() <<  "\t\t\tGf:\t" << pol_host.re << ' ' << pol_host.im;

		logger.debug() << "Even-odd: " << evenodd;
		logger.debug() << "Kappa: " << kappa_tmp;

		auto gm_code = get_device()->getGaugemomentumCode();
		gm_code->set_float_to_gaugemomentum_squarenorm_device(out, &tmp);
		tmp.dump(&resid);
		logger.debug() <<  "\t\t\tGm:\t" << resid;
	}

	auto run_kernel = [&](cl_kernel kernel) {
		//fermion_force(field, Y, X, out);
		//query work-sizes for kernel
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(kernel, &ls2, &gs2, &num_groups);
		//set arguments
		int clerr = clSetKernelArg(kernel, 0, sizeof(cl_mem), gf->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		clerr = clSetKernelArg(kernel, 1, sizeof(cl_mem), Y->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		clerr = clSetKernelArg(kernel, 2, sizeof(cl_mem), X->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		clerr = clSetKernelArg(kernel, 3, sizeof(cl_mem), out->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		clerr = clSetKernelArg(kernel, 4, sizeof(int), &evenodd);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		clerr = clSetKernelArg(kernel, 5, sizeof(hmc_float), &kappa_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		get_device()->enqueue_kernel( kernel , gs2, ls2);
	};

	run_kernel(fermion_force_eo_0);
	run_kernel(fermion_force_eo_1);
	run_kernel(fermion_force_eo_2);
	run_kernel(fermion_force_eo_3);

	if(logger.beDebug()) {
		Plain<hmc_float> force_tmp(1, get_device());
		auto gm_code = get_device()->getGaugemomentumCode();
		hmc_float resid;
		gm_code->set_float_to_gaugemomentum_squarenorm_device(out, &force_tmp);
		force_tmp.dump(&resid);
		logger.debug() <<  "\t\t\tFORCE [DET]:\t" << resid;

		if(resid != resid) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}


void hardware::code::Molecular_Dynamics::fermion_staggered_partial_force_device(const hardware::buffers::SU3 * gf, const hardware::buffers::SU3vec * A, const hardware::buffers::SU3vec * B, const hardware::buffers::Gaugemomentum * out, int evenodd) const
{
	using namespace hardware::buffers;

	//kernel prototype fermion_staggered_partial_force_eo(A, B, out, evenodd);
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(fermion_stagg_partial_force_eo, &ls2, &gs2, &num_groups);
	//set arguments
	
	int clerr = clSetKernelArg(fermion_stagg_partial_force_eo, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(fermion_stagg_partial_force_eo, 1, sizeof(cl_mem), A->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_stagg_partial_force_eo, 2, sizeof(cl_mem), B->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_stagg_partial_force_eo, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(fermion_stagg_partial_force_eo, 4, sizeof(int), &evenodd);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	get_device()->enqueue_kernel(fermion_stagg_partial_force_eo, gs2, ls2);

	if(logger.beDebug()) {
		Plain<hmc_float> force_tmp(1, get_device());
		auto gm_code = get_device()->getGaugemomentumCode();
		hmc_float resid;
		gm_code->set_float_to_gaugemomentum_squarenorm_device(out, &force_tmp);
		force_tmp.dump(&resid);
		logger.debug() <<  "\t\t\tFORCE_STAGG_PARTIAL:\t" << resid;

		if(resid != resid) {
			throw Print_Error_Message("calculation of force gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}


void hardware::code::Molecular_Dynamics::stout_smeared_fermion_force_device(std::vector<const hardware::buffers::SU3 *>&) const
{
	throw std::runtime_error("Not implemented!");
}

hardware::code::Molecular_Dynamics::Molecular_Dynamics(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device)
	: Opencl_Module(kernelParameters, device), md_update_gaugefield (0), gauge_force (0),
	  gauge_force_tlsym (0), fermion_force (0), fermion_force_eo_0(0), fermion_force_eo_1(0), fermion_force_eo_2(0), fermion_force_eo_3(0), stout_smear_fermion_force(0),
	  fermion_stagg_partial_force_eo(0),
	  gauge_force_tlsym_1 (0), gauge_force_tlsym_2 (0), gauge_force_tlsym_3 (0),
	  gauge_force_tlsym_4 (0), gauge_force_tlsym_5 (0), gauge_force_tlsym_6 (0),
	  gauge_force_tlsym_tmp(use_multipass_gauge_force_tlsym(device) ? new hardware::buffers::Matrix3x3(NDIM * get_vol4d(device->getLocalLatticeMemoryExtents()), device) : 0)
{
	fill_kernels();
}

hardware::code::Molecular_Dynamics::~Molecular_Dynamics()
{
	clear_kernels();

	if(gauge_force_tlsym_tmp) {
		delete gauge_force_tlsym_tmp;
	}
}
