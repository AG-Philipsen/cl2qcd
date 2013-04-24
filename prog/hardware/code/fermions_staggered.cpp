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

static std::string collect_build_options(hardware::Device *, const meta::Inputparameters& params)
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
	sources = get_device()->get_spinor_staggered_code()->get_sources() << ClSourcePackage(collect_build_options(get_device(), get_parameters()));

	M_staggered = 0;

	logger.debug() << "Create staggered fermion kernels...";

	M_staggered = createKernel("M_staggered") << sources << "fermionmatrix_staggered.cl" << "fermionmatrix_m_staggered.cl";

	return;
}

void hardware::code::Fermions_staggered::clear_kernels()
{
	logger.trace() << "clearing staggered fermion kernels...";
	cl_uint clerr = CL_SUCCESS;

	if(M_staggered) {
		clerr = clReleaseKernel(M_staggered);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
}

void hardware::code::Fermions_staggered::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);
}

//explicit fermionmatrix-kernel calling functions
void hardware::code::Fermions_staggered::M_staggered_device(const hardware::buffers::Plain<su3vec> * in, const hardware::buffers::Plain<su3vec> * out, const hardware::buffers::SU3 * gf, hmc_float mass) const
{
	//get mass
	hmc_float kappa_tmp;
	if(mass == ARG_DEF) kappa_tmp = get_parameters().get_kappa();
	else kappa_tmp = mass;

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

	clerr = clSetKernelArg(M_staggered, 3, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( M_staggered, gs2, ls2);
}


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
	//this is the same as in the function above
	//NOTE: 1 spinor has NC*NDIM = 12 complex entries
	if (in == "M_staggered") {
		//this kernel reads 9 spinors, 8 su3matrices and writes 1 spinor:
		return (C * 12 * (9 + 1) + C * 8 * R) * D * S;
	}
	return 0;
}

static int flop_dslash_per_site(const meta::Inputparameters & parameters)
{
	/** @NOTE: this is the "original" dslash without any simplifications, counting everything "full". this is a much too hight number!!
	   *  //this kernel performs for each eo site a 2*NDIM sum over (1 + gamma_mu) * su3matrix * spinor
	   *  //return  NDIM * 2 * ( get_parameters().get_flop_su3_spinor() + get_parameters().get_flop_gamma_spinor() ) ;
	  @NOTE: However, in 0911.3191 the dslash_eo is quoted to perform 1320 flop per site
	   *  If I count our implementation of the dslash-kernel, I get 1632 flop per site:
	   *  //the kernel performs 6 su3vec_acc, 2 su3_su3vec and 2 su3vec_complex in NDIM * 2 directions per site
	  */
	return NDIM * 2 * (meta::get_flop_su3_su3vec() * 2 + 6 * NC * 2 + 2 * NC * meta::get_flop_complex_mult() );
	/// I do not know for sure, but if one leaves out the 2 su3vec_complex operations in each directions one has almost 1320 flop per site. This would correspond to have the kappa in the diagonal matrix still.

}

uint64_t hardware::code::Fermions_staggered::get_flop_size(const std::string& in) const
{
	size_t S = get_spinorfieldsize(get_parameters());
	size_t Seo = get_eoprec_spinorfieldsize(get_parameters());
	if (in == "M_staggered") {
		//this kernel performs one dslash on each site and adds this to a spinor
		return S * (flop_dslash_per_site(get_parameters()) + NC * NDIM * meta::get_flop_complex_mult() + NC * NDIM * 2 );
	}
	return 0;
}

void hardware::code::Fermions_staggered::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, M_staggered);
}

hardware::code::Fermions_staggered::Fermions_staggered(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module(params, device),
	  M_staggered(0)
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
