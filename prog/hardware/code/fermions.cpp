#include "fermions.hpp"

#include "../../logger.hpp"
#include "../../meta/util.hpp"
#include "../device.hpp"
#include "spinors.hpp"

#include <cassert>
#include <cmath>

using namespace std;

void hardware::code::Fermions::fill_kernels()
{
	sources = get_basic_sources() << "operations_geometry.cl" << "operations_complex.cl" << "types_fermions.h" << "operations_matrix_su3.cl" << "operations_matrix.cl" << "operations_gaugefield.cl" << "operations_su3vec.cl" << "operations_spinor.cl" << "spinorfield.cl";
	if(get_parameters().get_use_eo()) {
		sources = sources << "operations_spinorfield_eo.cl";
	}

	logger.debug() << "Create Fermions kernels...";
	
	if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
		M_wilson = createKernel("M_wilson") << sources << "fermionmatrix.cl" << "fermionmatrix_m.cl";
	} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
		M_tm_plus = createKernel("M_tm_plus") << sources << "fermionmatrix.cl" << "fermionmatrix_m_tm_plus.cl";
		M_tm_minus = createKernel("M_tm_minus") << sources << "fermionmatrix.cl" << "fermionmatrix_m_tm_minus.cl";
	} else if(get_parameters().get_fermact() == meta::Inputparameters::clover) {
		throw Print_Error_Message("no kernels for CLOVER-discretization implemented yet, aborting... ", __FILE__, __LINE__);
	} else {
		throw Print_Error_Message("there was a problem with which fermion-discretization to use, aborting... ", __FILE__, __LINE__);
	}
	
	gamma5 = createKernel("gamma5") << sources << "fermionmatrix.cl" << "fermionmatrix_gamma5.cl";
	//Kernels needed if eoprec is used
	if(get_parameters().get_use_eo() == true) {
		if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
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
		if (get_parameters().get_use_merge_kernels_fermion() == true) {
			dslash_AND_gamma5_eo = createKernel("dslash_AND_gamma5_eo") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_dslash_AND_gamma5.cl";
			dslash_AND_M_tm_inverse_sitediagonal_eo = createKernel("dslash_AND_M_tm_inverse_sitediagonal_eo") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_dslash_AND_M_tm_inverse_sitediagonal.cl";
			dslash_AND_M_tm_inverse_sitediagonal_minus_eo = createKernel("dslash_AND_M_tm_inverse_sitediagonal_minus_eo") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_dslash_AND_M_tm_inverse_sitediagonal_minus.cl";
			M_tm_sitediagonal_AND_gamma5_eo = createKernel("M_tm_sitediagonal_AND_gamma5_eo") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m_merged.cl";
			M_tm_sitediagonal_minus_AND_gamma5_eo = createKernel("M_tm_sitediagonal_minus_AND_gamma5_eo") << sources << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m_merged.cl";
		}
	}
}

void hardware::code::Fermions::clear_kernels()
{
	logger.trace() << "clearing fermion kernels...";
	cl_uint clerr = CL_SUCCESS;

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

	if(get_parameters().get_use_eo()) {
		clerr = clReleaseKernel(dslash_eo);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(_dslash_eo_boundary);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(_dslash_eo_inner);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(gamma5_eo);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			clerr = clReleaseKernel(M_tm_sitediagonal);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
			clerr = clReleaseKernel(M_tm_inverse_sitediagonal);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
			clerr = clReleaseKernel(M_tm_sitediagonal_minus);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
			clerr = clReleaseKernel(M_tm_inverse_sitediagonal_minus);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		}
	}
}

void hardware::code::Fermions::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);
}


//explicit fermionmatrix-kernel calling functions
void hardware::code::Fermions::M_wilson_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa) const
{
	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = get_parameters().get_kappa();
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
	if(kappa == ARG_DEF) kappa_tmp = get_parameters().get_kappa();
	else kappa_tmp = kappa;

	//get mu
	hmc_float mubar_tmp;
	if(mubar == ARG_DEF) mubar_tmp = meta::get_mubar(get_parameters());
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
	if(kappa == ARG_DEF) kappa_tmp = get_parameters().get_kappa();
	else kappa_tmp = kappa;

	//get mu
	hmc_float mubar_tmp;
	if(mubar == ARG_DEF) mubar_tmp = meta::get_mubar(get_parameters());
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

//merged fermionmatrix-functions with eoprec
//void hardware::code::Fermions::Aee_AND_gamma5_eo(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa , hmc_float mubar )
//{
//	int even = EVEN;
//	int odd = ODD;
//
//	auto spinor_code = get_device()->get_spinor_code();
//
//	/**
//	 * This is the even-odd preconditioned fermion matrix with the
//	 * non-trivial inversion on the even sites (see DeGran/DeTar p. 174).
//	 * If one has fermionmatrix
//	 *  M = R + D,
//	 * then Aee is:
//	 * Aee = R_e - D_eo R_o_inv D_oe
//	 */
//	if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
//		//in this case, the diagonal matrix is just 1 and falls away.
//		//this case has not been adjusted for the merged kernels yet...
//		logger.warn() << "No merged kernels implemented for pure Wilson case!";
//		dslash_eo_device(in, &clmem_tmp_eo_1, gf, odd, kappa);
//		dslash_eo_device(&clmem_tmp_eo_1, out, gf, even, kappa);
//		spinor_code->saxpy_eoprec_device(out, in, &clmem_one, out);
//		gamma5_eo_device(out);
//	} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
//		/*
//		dslash_eo_device(in, &clmem_tmp_eo_1, gf, odd, kappa);
//		M_tm_inverse_sitediagonal_device(&clmem_tmp_eo_1, &clmem_tmp_eo_2, mubar);
//		*/
//		dslash_AND_M_tm_inverse_sitediagonal_eo_device(in, &clmem_tmp_eo_2, gf, odd, kappa, mubar);
//		/*
//		dslash_eo_device(&clmem_tmp_eo_2, out, gf, even, kappa);
//		gamma5_eo_device(out);
//		*/
//		dslash_AND_gamma5_eo_device(&clmem_tmp_eo_2, out, gf, even, kappa);
//		/*
//		M_tm_sitediagonal_device(in, &clmem_tmp_eo_1, mubar);
//		gamma5_eo_device(&clmem_tmp_eo_1);
//		*/
//		M_tm_sitediagonal_AND_gamma5_eo_device(in, &clmem_tmp_eo_1, mubar);
//		spinor_code->saxpy_eoprec_device(out, &clmem_tmp_eo_1, &clmem_one, out);
//	}
//}

///**
// *  This is the equivalent of the above function, but for the lower
// *  flavour, which essentially means mu -> -mu in the tm-case and
// *  no changes in the meta::Inputparameters::wilson case.
// */
//void hardware::code::Fermions::Aee_minus_AND_gamma5_eo(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa , hmc_float mubar )
//{
//	int even = EVEN;
//	int odd = ODD;
//
//	auto spinor_code = get_device()->get_spinor_code();
//
//	/**
//	 * This is the even-odd preconditioned fermion matrix with the
//	 * non-trivial inversion on the even sites (see DeGran/DeTar p. 174).
//	 * If one has fermionmatrix
//	 *  M = R + D,
//	 * then Aee is:
//	 * Aee = R_e - D_eo R_o_inv D_oe
//	 */
//	if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
//		//in this case, the diagonal matrix is just 1 and falls away.
//		//this case has not been adjusted for the merged kernels yet...
//		logger.warn() << "No merged kernels implemented for pure Wilson case!";
//		dslash_eo_device(in, &clmem_tmp_eo_1, gf, odd, kappa);
//		dslash_eo_device(&clmem_tmp_eo_1, out, gf, even, kappa);
//		spinor_code->saxpy_eoprec_device(out, in, &clmem_one, out);
//		gamma5_eo_device(out);
//	} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
//		/*
//		dslash_eo_device(in, &clmem_tmp_eo_1, gf, odd, kappa);
//		M_tm_inverse_sitediagonal_minus_device(&clmem_tmp_eo_1, &clmem_tmp_eo_2, mubar);
//		*/
//		dslash_AND_M_tm_inverse_sitediagonal_minus_eo_device(in, &clmem_tmp_eo_2, gf, odd, kappa, mubar);
//		/*
//		dslash_eo_device(&clmem_tmp_eo_2, out, gf, even, kappa);
//		gamma5_eo_device(out);
//		*/
//		dslash_AND_gamma5_eo_device(&clmem_tmp_eo_2, out, gf, even, kappa);
//		/*
//		M_tm_sitediagonal_minus_device(in, &clmem_tmp_eo_1, mubar);
//		gamma5_eo_device(&clmem_tmp_eo_1);
//		*/
//		M_tm_sitediagonal_minus_AND_gamma5_eo_device(in, &clmem_tmp_eo_1, mubar);
//		spinor_code->saxpy_eoprec_device(out, &clmem_tmp_eo_1, &clmem_one, out);
//	}
//}


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

void hardware::code::Fermions::dslash_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa) const
{
	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = get_parameters().get_kappa();
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
	if(kappa == ARG_DEF) kappa_tmp = get_parameters().get_kappa();
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
	if(kappa == ARG_DEF) kappa_tmp = get_parameters().get_kappa();
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

void hardware::code::Fermions::dslash_AND_gamma5_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa) const
{
	//get kappa
	hmc_float kappa_tmp;
	if(kappa == ARG_DEF) kappa_tmp = get_parameters().get_kappa();
	else kappa_tmp = kappa;

	cl_int eo = evenodd;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(dslash_AND_gamma5_eo, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(dslash_AND_gamma5_eo, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_gamma5_eo, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_gamma5_eo, 2, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_gamma5_eo, 3, sizeof(cl_int), &eo);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_gamma5_eo, 4, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(dslash_AND_gamma5_eo , gs2, ls2);
}

void hardware::code::Fermions::dslash_AND_M_tm_inverse_sitediagonal_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa, hmc_float mubar) const
{
	//get kappa
	hmc_float kappa_tmp, mubar_tmp;
	if(kappa == ARG_DEF) kappa_tmp = get_parameters().get_kappa();
	else kappa_tmp = kappa;
	if(mubar == ARG_DEF) mubar_tmp = meta::get_mubar(get_parameters());
	else mubar_tmp = kappa;

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

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_eo, 4, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_eo, 5, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(dslash_AND_M_tm_inverse_sitediagonal_eo , gs2, ls2);
}

void hardware::code::Fermions::dslash_AND_M_tm_inverse_sitediagonal_minus_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa, hmc_float mubar) const
{
	//get kappa
	hmc_float kappa_tmp, mubar_tmp;
	if(kappa == ARG_DEF) kappa_tmp = get_parameters().get_kappa();
	else kappa_tmp = kappa;
	if(mubar == ARG_DEF) mubar_tmp = meta::get_mubar(get_parameters());
	else mubar_tmp = kappa;

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

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, 4, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, 5, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(dslash_AND_M_tm_inverse_sitediagonal_minus_eo , gs2, ls2);
}

void hardware::code::Fermions::M_tm_inverse_sitediagonal_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar) const
{
	//get mu
	hmc_float mubar_tmp;
	if(mubar == ARG_DEF) mubar_tmp = meta::get_mubar(get_parameters());
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
	if(mubar == ARG_DEF) mubar_tmp = meta::get_mubar(get_parameters());
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
	if(mubar == ARG_DEF) mubar_tmp = meta::get_mubar(get_parameters());
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
	if(mubar == ARG_DEF) mubar_tmp = meta::get_mubar(get_parameters());
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
	if(mubar == ARG_DEF) mubar_tmp = meta::get_mubar(get_parameters());
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
	if(mubar == ARG_DEF) mubar_tmp = meta::get_mubar(get_parameters());
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
	if (in == "dslash_AND_gamma5_eo") {
		//the dslash kernel reads 8 spinors, 8 su3matrices and writes 1 spinor:
		const unsigned int dirs = 4;
		//the gamma5 kernel reads 1 spinor and writes 1 spinor:
		//the merged kernel reads 8 spinors, 8 su3matrices and writes 1 spinor, thus it is the same as the dslash
		return  (C * 12 * (2 * dirs + 1) + C * 2 * dirs * R) * D * Seo;
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

uint64_t hardware::code::Fermions::get_flop_size(const std::string& in) const
{
	size_t S = get_spinorfieldsize(get_parameters());
	size_t Seo = get_eoprec_spinorfieldsize(get_parameters());
	if (in == "M_wilson") {
		//this kernel performs one dslash on each site and adds this to a spinor
		return S * (flop_dslash_per_site(get_parameters()) + NC * NDIM * meta::get_flop_complex_mult() + NC * NDIM * 2 );
	}
	if (in == "gamma5") {
		//this kernel performs ND*NC*2/2 real mults
		return S * NDIM * NC;
	}
	if (in == "M_tm_plus") {
		//this kernel performs ND*NC complex mults and one dslash on each site and adds the results
		return S * (flop_dslash_per_site(get_parameters()) + NC * NDIM * meta::get_flop_complex_mult() + NC * NDIM * 2 );
	}
	if (in == "M_tm_minus") {
		//this kernel performs ND*NC complex mults and one dslash on each site and adds the results
		return S * (flop_dslash_per_site(get_parameters()) + NC * NDIM * meta::get_flop_complex_mult() + NC * NDIM * 2 );
	}
	if (in == "gamma5_eo") {
		//this kernel performs ND*NC*2/2 real mults
		return Seo * NDIM * NC;
	}
	if (in == "M_tm_sitediagonal") {
		//this kernel performs ND*NC complex mults
		return Seo * ( NC * NDIM * meta::get_flop_complex_mult() );
	}
	if (in == "M_tm_inverse_sitediagonal") {
		//this kernel performs ND*NC complex mults and ND*NC*2 real mults
		return Seo * ( NC * NDIM * meta::get_flop_complex_mult() + NC * NDIM * 2  );
	}
	if (in == "M_tm_sitediagonal_minus") {
		//this kernel performs ND*NC complex mults
		return Seo * ( NC * NDIM * meta::get_flop_complex_mult() );
	}
	if (in == "M_tm_inverse_sitediagonal_minus") {
		//this kernel performs ND*NC complex mults and ND*NC*2 real mults
		return Seo * ( NC * NDIM * meta::get_flop_complex_mult() + NC * NDIM * 2 );
	}
	if (in == "dslash_eo") {
		return Seo * flop_dslash_per_site(get_parameters());
	}
	if (in == "dslash_AND_gamma5_eo") {
		return Seo * flop_dslash_per_site(get_parameters()) +  Seo * NDIM * NC;
	}
	if (in == "dslash_AND_M_tm_inverse_sitediagonal_eo") {
		return Seo * flop_dslash_per_site(get_parameters()) + Seo * ( NC * NDIM * meta::get_flop_complex_mult() + NC * NDIM * 2  );
	}
	if (in == "dslash_AND_M_tm_inverse_sitediagonal_minus_eo") {
		return Seo * flop_dslash_per_site(get_parameters()) + Seo * ( NC * NDIM * meta::get_flop_complex_mult() + NC * NDIM * 2  );
	}
	if (in == "M_tm_sitediagonal_AND_gamma5_eo") {
		//this kernel performs ND*NC complex mults and  ND*NC*2/2 real mults
		return Seo * ( NC * NDIM * meta::get_flop_complex_mult() ) + Seo * NDIM * NC ;
	}
	if (in == "M_tm_sitediagonal_minus_AND_gamma5_eo") {
		//this kernel performs ND*NC complex mults  ND*NC*2/2 real mults
		return Seo * ( NC * NDIM * meta::get_flop_complex_mult() ) +  Seo * NDIM * NC;
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
	Opencl_Module::print_profiling(filename, dslash_AND_gamma5_eo);
	Opencl_Module::print_profiling(filename, dslash_AND_M_tm_inverse_sitediagonal_eo);
	Opencl_Module::print_profiling(filename, dslash_AND_M_tm_inverse_sitediagonal_minus_eo);
	Opencl_Module::print_profiling(filename, M_tm_sitediagonal_AND_gamma5_eo);
	Opencl_Module::print_profiling(filename, M_tm_sitediagonal_minus_AND_gamma5_eo);
}
hardware::code::Fermions::Fermions(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module(params, device),
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
	  dslash_AND_gamma5_eo(0),
	  dslash_AND_M_tm_inverse_sitediagonal_eo(0),
	  dslash_AND_M_tm_inverse_sitediagonal_minus_eo(0),
	  M_tm_sitediagonal_AND_gamma5_eo(0),
	  M_tm_sitediagonal_minus_AND_gamma5_eo(0)
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
