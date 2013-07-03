/** @file
 * Implementation of the hardware::code::fermions_staggered class
 */

#include "fermions_staggered.hpp"

#include "../../logger.hpp"
#include "../../meta/util.hpp"
#include "../device.hpp"
#include "spinors_staggered.hpp"
#include "spinors.hpp"

#include <cassert>
#include <cmath>

using namespace std;

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params);

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params)
{
	std::ostringstream options;

	options.precision(16);

	//CP: These are the BCs in spatial and temporal direction
	hmc_float tmp_spatial = (params.get_theta_fermion_spatial() * PI) / ( (hmc_float) params.get_nspace());
	hmc_float tmp_temporal = (params.get_theta_fermion_temporal() * PI) / ( (hmc_float) params.get_ntime());
	//BC: on the corners in each direction: exp(i theta*PI) -> on each site
	//    exp(i theta*PI /LATEXTENSION) = cos(tmp2) + isin(tmp2)

	options << " -D SPATIAL_RE=" << cos(tmp_spatial);
	options << " -D MSPATIAL_RE=" << -cos(tmp_spatial);
	options << " -D SPATIAL_IM=" << sin(tmp_spatial);
	options << " -D MSPATIAL_IM=" << -sin(tmp_spatial);

	options << " -D TEMPORAL_RE=" << cos(tmp_temporal);
	options << " -D MTEMPORAL_RE=" << -cos(tmp_temporal);
	options << " -D TEMPORAL_IM=" << sin(tmp_temporal);
	options << " -D MTEMPORAL_IM=" << -sin(tmp_temporal);
  	
	return options.str();
}

void hardware::code::Fermions_staggered::fill_kernels()
{
	sources = get_device()->get_spinor_staggered_code()->get_sources() << "operations_staggered.cl" << ClSourcePackage(collect_build_options(get_device(), get_parameters()));

	logger.debug() << "Create staggered fermion kernels...";

	M_staggered = createKernel("M_staggered") << sources << "fermionmatrix_staggered_DKS_local.cl" << "fermionmatrix_staggered_M.cl";
	/////////////////////////////////////////////////
	/////////// EVEN-ODD PRECONDITIONING ////////////
	/////////////////////////////////////////////////
	if(get_parameters().get_use_eo()){
	  D_KS_eo = createKernel("D_KS_eo") << sources << "fermionmatrix_staggered_eo_DKS_local.cl" << "fermionmatrix_staggered_eo_DKS.cl";
	} else {
	  D_KS_eo = 0;
	}

	return;
}

void hardware::code::Fermions_staggered::clear_kernels()
{
	logger.trace() << "clearing staggered fermion kernels...";
	cl_uint clerr = CL_SUCCESS;

	if(M_staggered) {
		clerr = clReleaseKernel(M_staggered);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		/////////////////////////////////////////////////
		/////////// EVEN-ODD PRECONDITIONING ////////////
		/////////////////////////////////////////////////
		if(get_parameters().get_use_eo()){
			clerr = clReleaseKernel(D_KS_eo);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		}
	}
}

void hardware::code::Fermions_staggered::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//explicit fermionmatrix-kernel calling functions
void hardware::code::Fermions_staggered::M_staggered_device(const hardware::buffers::Plain<su3vec> * in, const hardware::buffers::Plain<su3vec> * out, const hardware::buffers::SU3 * gf, hmc_float mass) const
{
	//get mass
	hmc_float mass_tmp;
	if(mass == ARG_DEF) mass_tmp = get_parameters().get_kappa();
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
	size_t D = meta::get_float_size(get_parameters());
	//this returns the number of entries in an su3-matrix
	size_t R = meta::get_mat_size(get_parameters());
	//this is the number of spinors in the system (or number of sites)
	size_t S = get_spinorfieldsize(get_parameters());
	size_t Seo = get_eoprec_spinorfieldsize(get_parameters());
	//factor for complex numbers
	int C = 2;
	//NOTE: 1 spinor has NC = 3 complex entries
	if (in == "M_staggered") {
		//this kernel reads 9 su3vec, 8 su3matrices and writes 1 su3vec per site:
		return (C * NC * (9 + 1) + C * 8 * R) * D * S;
	}
	if (in == "D_KS_eo") {
		//Not yet counted
		return 100000000000000000000000;
	}
	return 0;
}

/**
 * 
 * This function returns the number of flops that is needed to make the standard staggered
 * Dirac operator act onto a field (606 flops). Since in general the mass in the simulation is not
 * zero (and there are no check on that), the mass-term in the d_slash is always taken
 * into account (For standard Dirac operator we intend M = D_KS + m).
 * 
 * @NOTE In this function (as in the whole code) the staggered phases are not included
 *       in links and, then, we have some flops in addition to take into account.
 * @NOTE Here, we do not count the flops that are performed for example in the mulptiplications
 *       done to impose the boundary conditions. These flops will be taken into account in the
 *       function get_flop_size.
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
 *  + We have to sum the 2 rsulting su3vec in each direction
 *  + To multiply by (0.5*eta) the result of the sum in each direction, a real times su3vec is required
 *  + The mass term is a real times su3vec
 *  + Then we have to sum the 5 su3vec (one for each direction + the mass term)
 * 
 */
static int flop_dslash_staggered_per_site(const meta::Inputparameters & parameters)
{
	return NDIM * (2 * meta::get_flop_su3_su3vec() + NC * 2 + NC * 2) + NC * 2 + 4 * NC * 2;
}

uint64_t hardware::code::Fermions_staggered::get_flop_size(const std::string& in) const
{
	size_t S = get_spinorfieldsize(get_parameters());
	size_t Seo = get_eoprec_spinorfieldsize(get_parameters());
	if (in == "M_staggered") {
		//this kernel performs one dslash on each site. To do that it also impose
		//the boundary conditions multiplying twice in each direction an su3vec by a complex
		return S * (flop_dslash_staggered_per_site(get_parameters()) +
		            2* NC * NDIM * meta::get_flop_complex_mult());
	}
	if (in == "D_KS_eo") {
		//Not yet counted
		return 100000000000000000000000;
	}
	return 0;
}

void hardware::code::Fermions_staggered::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, M_staggered);
	Opencl_Module::print_profiling(filename, D_KS_eo);
}

hardware::code::Fermions_staggered::Fermions_staggered(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module(params, device)
{
	fill_kernels();
};

hardware::code::Fermions_staggered::~Fermions_staggered()
{
	clear_kernels();
}

ClSourcePackage hardware::code::Fermions_staggered::get_sources() const noexcept
{
	return sources;
}
