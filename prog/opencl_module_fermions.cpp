#include "opencl_module_fermions.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

/**
 * What follows are functions that call opencl_fermions-class-functions.
 * This is needed to be able to pass different fermionmatrices as
 *  arguments to class-functions.
 */
void M::operator()(cl_mem in, cl_mem out, cl_mem gf) const
{
	that->M(in, out, gf);
}
cl_ulong M::get_Flops() const
{
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			return that->get_flop_size("M_wilson");
		case TWISTEDMASS:
			return that->get_flop_size("M_tm_plus");
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
}
cl_ulong M::get_Bytes() const
{
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			return that->get_read_write_size("M_wilson");
		case TWISTEDMASS:
			return that->get_read_write_size("M_tm_plus");
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
}

void Qplus::operator()(cl_mem in, cl_mem out, cl_mem gf) const
{
	that->Qplus(in, out, gf);
}
cl_ulong Qplus::get_Flops() const
{
	cl_ulong res;
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			res = that->get_flop_size("M_wilson");
			break;
		case TWISTEDMASS:
			res = that->get_flop_size("M_tm_plus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
	res += that->get_flop_size("gamma5");
	return res;
}
cl_ulong Qplus::get_Bytes() const
{
	cl_ulong res;
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			res = that->get_read_write_size("M_wilson");
			break;
		case TWISTEDMASS:
			res = that->get_read_write_size("M_tm_plus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
	res += that->get_read_write_size("gamma5");
	return res;
}

void Qminus::operator()(cl_mem in, cl_mem out, cl_mem gf) const
{
	that->Qminus(in, out, gf);
}
cl_ulong Qminus::get_Flops() const
{
	cl_ulong res;
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			res = that->get_flop_size("M_wilson");
			break;
		case TWISTEDMASS:
			res = that->get_flop_size("M_tm_minus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
	res += that->get_flop_size("gamma5");
	return res;
}
cl_ulong Qminus::get_Bytes() const
{
	cl_ulong res;
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			res = that->get_read_write_size("M_wilson");
			break;
		case TWISTEDMASS:
			res = that->get_read_write_size("M_tm_minus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
	res += that->get_read_write_size("gamma5");
	return res;
}

void QplusQminus::operator()(cl_mem in, cl_mem out, cl_mem gf) const
{
	that->QplusQminus(in, out, gf);
}
cl_ulong QplusQminus::get_Flops() const
{
	cl_ulong res;
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			res = 2 * that->get_flop_size("M_wilson");
			break;
		case TWISTEDMASS:
			res = that->get_flop_size("M_tm_plus");
			res += that->get_flop_size("M_tm_minus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
	res += 2 * that->get_flop_size("gamma5");
	return res;
}
cl_ulong QplusQminus::get_Bytes() const
{
	cl_ulong res;
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			res = 2 * that->get_read_write_size("M_wilson");
			break;
		case TWISTEDMASS:
			res = that->get_read_write_size("M_tm_plus");
			res += that->get_read_write_size("M_tm_minus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
	res += 2 * that->get_read_write_size("gamma5");
	return res;
}

void Aee::operator()(cl_mem in, cl_mem out, cl_mem gf) const
{
	that->Aee(in, out, gf);
}
cl_ulong Aee::get_Flops() const
{
	cl_ulong res;
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			res = 2 * that->get_flop_size("dslash_eoprec");
			res += that->get_flop_size("saxpy_eoprec");
			break;
		case TWISTEDMASS:
			res = 2 * that->get_flop_size("dslash_eoprec");
			res += that->get_flop_size("M_tm_inverse_sitediagonal");
			res += that->get_flop_size("M_tm_sitediagonal");
			res += that->get_flop_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
	return res;
}
cl_ulong Aee::get_Bytes() const
{
	cl_ulong res;
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			res = 2 * that->get_read_write_size("dslash_eoprec");
			res += that->get_read_write_size("saxpy_eoprec");
			break;
		case TWISTEDMASS:
			res = 2 * that->get_read_write_size("dslash_eoprec");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal");
			res += that->get_read_write_size("M_tm_sitediagonal");
			res += that->get_read_write_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
	return res;
}

void Qplus_eoprec::operator()(cl_mem in, cl_mem out, cl_mem gf) const
{
	that->Qplus_eoprec(in, out, gf);
}
cl_ulong Qplus_eoprec::get_Flops() const
{
	cl_ulong res;
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			res = 2 * that->get_flop_size("dslash_eoprec");
			res += that->get_flop_size("saxpy_eoprec");
			break;
		case TWISTEDMASS:
			res = 2 * that->get_flop_size("dslash_eoprec");
			res += that->get_flop_size("M_tm_inverse_sitediagonal");
			res += that->get_flop_size("M_tm_sitediagonal");
			res += that->get_flop_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
	res += that->get_flop_size("gamma5_eoprec");
	return res;
}
cl_ulong Qplus_eoprec::get_Bytes() const
{
	cl_ulong res;
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			res = 2 * that->get_read_write_size("dslash_eoprec");
			res += that->get_read_write_size("saxpy_eoprec");
			break;
		case TWISTEDMASS:
			res = 2 * that->get_read_write_size("dslash_eoprec");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal");
			res += that->get_read_write_size("M_tm_sitediagonal");
			res += that->get_read_write_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
	res += that->get_flop_size("gamma5_eoprec");
	return res;
}

void Qminus_eoprec::operator()(cl_mem in, cl_mem out, cl_mem gf) const
{
	that->Qminus_eoprec(in, out, gf);
}
cl_ulong Qminus_eoprec::get_Flops() const
{
	cl_ulong res;
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			res = 2 * that->get_flop_size("dslash_eoprec");
			res += that->get_flop_size("saxpy_eoprec");
			break;
		case TWISTEDMASS:
			res = 2 * that->get_flop_size("dslash_eoprec");
			res += that->get_flop_size("M_tm_inverse_sitediagonal_minus");
			res += that->get_flop_size("M_tm_sitediagonal_minus");
			res += that->get_flop_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
	res += that->get_flop_size("gamma5_eoprec");
	return res;
}
cl_ulong Qminus_eoprec::get_Bytes() const
{
	cl_ulong res;
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			res = 2 * that->get_read_write_size("dslash_eoprec");
			res += that->get_read_write_size("saxpy_eoprec");
			break;
		case TWISTEDMASS:
			res = 2 * that->get_read_write_size("dslash_eoprec");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal_minus");
			res += that->get_read_write_size("M_tm_sitediagonal_minus");
			res += that->get_read_write_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
	res += that->get_flop_size("gamma5_eoprec");
	return res;
}

void QplusQminus_eoprec::operator()(cl_mem in, cl_mem out, cl_mem gf) const
{
	that->QplusQminus_eoprec(in, out, gf);
}
cl_ulong QplusQminus_eoprec::get_Flops() const
{
	cl_ulong res;
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			res = 2 * that->get_flop_size("dslash_eoprec");
			res += that->get_flop_size("saxpy_eoprec");
			res *= 2;
			break;
		case TWISTEDMASS:
			res = 4 * that->get_flop_size("dslash_eoprec");
			res += that->get_flop_size("M_tm_inverse_sitediagonal");
			res += that->get_flop_size("M_tm_sitediagonal");
			res += that->get_flop_size("M_tm_inverse_sitediagonal_minus");
			res += that->get_flop_size("M_tm_sitediagonal_minus");
			res += 2 * that->get_flop_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
	res += 2 * that->get_flop_size("gamma5_eoprec");
	return res;
}
cl_ulong QplusQminus_eoprec::get_Bytes() const
{
	cl_ulong res;
	switch(that->get_parameters()->get_fermact()) {
		case WILSON:
			res = 2 * that->get_read_write_size("dslash_eoprec");
			res += that->get_read_write_size("saxpy_eoprec");
			res *= 2;
			break;
		case TWISTEDMASS:
			res = 4 * that->get_read_write_size("dslash_eoprec");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal");
			res += that->get_read_write_size("M_tm_sitediagonal");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal_minus");
			res += that->get_read_write_size("M_tm_sitediagonal_minus");
			res += 2 * that->get_read_write_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "WILSON or TWISTEDMASS", that->get_parameters()->get_fermact());
	}
	res += 2 * that->get_flop_size("gamma5_eoprec");
	return res;
}


void Opencl_Module_Fermions::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module_Spinors::fill_collect_options(collect_options);

	switch (get_parameters()->get_fermact()) {
		case TWISTEDMASS :
			*collect_options << " -D_TWISTEDMASS_";
			break;
		case CLOVER :
			*collect_options << " -D_CLOVER_";
			break;
	}

	//CP: give kappa and its negative value
	hmc_float kappa_tmp = get_parameters()->get_kappa();
	*collect_options << " -DKAPPA=" << kappa_tmp;
	*collect_options << " -DMKAPPA=" << -kappa_tmp;
	//CP: These are the kappas including BC in spatial and temporal direction
	hmc_float tmp_spatial = (get_parameters()->get_theta_fermion_spatial() * PI) / ( (hmc_float) get_parameters()->get_ns());
	hmc_float tmp_temporal = (get_parameters()->get_theta_fermion_temporal() * PI) / ( (hmc_float) get_parameters()->get_nt());
	//BC: on the corners in each direction: exp(i theta) -> on each site exp(i theta*PI /LATEXTENSION) = cos(tmp2) + isin(tmp2)
	*collect_options << " -DKAPPA_SPATIAL_RE=" << kappa_tmp*cos(tmp_spatial);
	*collect_options << " -DMKAPPA_SPATIAL_RE=" << -kappa_tmp*cos(tmp_spatial);
	*collect_options << " -DKAPPA_SPATIAL_IM=" << kappa_tmp*sin(tmp_spatial);
	*collect_options << " -DMKAPPA_SPATIAL_IM=" << -kappa_tmp*sin(tmp_spatial);

	*collect_options << " -DKAPPA_TEMPORAL_RE=" << kappa_tmp*cos(tmp_temporal);
	*collect_options << " -DMKAPPA_TEMPORAL_RE=" << -kappa_tmp*cos(tmp_temporal);
	*collect_options << " -DKAPPA_TEMPORAL_IM=" << kappa_tmp*sin(tmp_temporal);
	*collect_options << " -DMKAPPA_TEMPORAL_IM=" << -kappa_tmp*sin(tmp_temporal);

	switch (get_parameters()->get_fermact()) {
		case TWISTEDMASS :
			*collect_options << " -DMU=" << get_parameters()->get_mu();
			get_parameters()->calc_mubar();
			*collect_options << " -DMUBAR=" << get_parameters()->get_mubar();
			get_parameters()->set_mubar_negative();
			*collect_options << " -DMMUBAR=" << get_parameters()->get_mubar();
			get_parameters()->set_mubar_negative();
			break;
		case CLOVER :
			*collect_options << " -DCSW=" << get_parameters()->get_csw();
			break;
	}

	return;
}


void Opencl_Module_Fermions::fill_buffers()
{

	Opencl_Module_Spinors::fill_buffers();

	int clerr = CL_SUCCESS;

	int spinorfield_size = sizeof(spinor) * get_parameters()->get_spinorfieldsize();
	int eoprec_spinorfield_size = sizeof(spinor) * get_parameters()->get_eoprec_spinorfieldsize();
	int complex_size = sizeof(hmc_complex);
	int float_size = sizeof(hmc_float);

	logger.debug() << "init general spinorfield-buffers";
	clmem_inout = create_rw_buffer(spinorfield_size);
	clmem_tmp = create_rw_buffer(spinorfield_size);
	clmem_source = create_rw_buffer(spinorfield_size);

	logger.debug() << "init solver spinorfield-buffers";
	///@todo some buffers can be saved here if only cg is used
	if(get_parameters()->get_use_eo() == false) {
		//these are only used in a non-eoprec solver
		clmem_rn = create_rw_buffer(spinorfield_size);
		clmem_rhat = create_rw_buffer(spinorfield_size);
		clmem_v = create_rw_buffer(spinorfield_size);
		clmem_p = create_rw_buffer(spinorfield_size);
		clmem_s = create_rw_buffer(spinorfield_size);
		clmem_t = create_rw_buffer(spinorfield_size);
		clmem_aux = create_rw_buffer(spinorfield_size);
	} else {
		//LZ only use the following if we want to apply even odd preconditioning
		logger.debug() << "init solver eoprec-spinorfield-buffers";
		clmem_rn_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_rhat_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_v_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_p_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_s_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_t_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_aux_eoprec = create_rw_buffer(eoprec_spinorfield_size);

	} //end if: eoprec

	if(get_parameters()->get_use_eo() == true) {
		logger.debug() << "init general eoprec-spinorfield-buffers";
		clmem_inout_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_source_even = create_rw_buffer(eoprec_spinorfield_size);
		clmem_source_odd = create_rw_buffer(eoprec_spinorfield_size);
		clmem_tmp_eoprec_1 = create_rw_buffer(eoprec_spinorfield_size);
		//this field is used only with twistedmass
		if(get_parameters()->get_fermact() == TWISTEDMASS)
			clmem_tmp_eoprec_2 = create_rw_buffer(eoprec_spinorfield_size);
	}


	logger.debug() << "create buffers for complex and real numbers";
	clmem_rho = create_rw_buffer(complex_size);
	clmem_rho_next = create_rw_buffer(complex_size);
	clmem_alpha = create_rw_buffer(complex_size);
	clmem_omega = create_rw_buffer(complex_size);
	clmem_beta = create_rw_buffer(complex_size);
	clmem_tmp1 = create_rw_buffer(complex_size);
	clmem_tmp2 = create_rw_buffer(complex_size);
	clmem_one = create_rw_buffer(complex_size);
	clmem_minusone = create_rw_buffer(complex_size);
	clmem_resid = create_rw_buffer(float_size);
	clmem_trueresid = create_rw_buffer(float_size);

	logger.debug() << "write contents to some buffers";
	hmc_complex one = hmc_complex_one;
	hmc_complex minusone = hmc_complex_minusone;
	clerr = clEnqueueWriteBuffer(get_queue(), clmem_one, CL_TRUE, 0, complex_size, &one, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);

	clerr = clEnqueueWriteBuffer(get_queue(), clmem_minusone, CL_TRUE, 0, complex_size, &minusone, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);

	hmc_complex zero = hmc_complex_zero;
	clerr = clEnqueueWriteBuffer(get_queue(), clmem_resid, CL_TRUE, 0, float_size, &zero, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);
	clerr = clEnqueueWriteBuffer(get_queue(), clmem_trueresid, CL_TRUE, 0, float_size, &zero, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);



	return;
}

void Opencl_Module_Fermions::fill_kernels()
{
	Opencl_Module_Spinors::fill_kernels();

	M_wilson = 0;
	M_tm_plus = 0;
	M_tm_minus = 0;

	logger.debug() << "Create fermion kernels...";
	if(get_parameters()->get_fermact() == WILSON) {
		M_wilson = createKernel("M_wilson") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_m.cl";
	} else if(get_parameters()->get_fermact() == TWISTEDMASS) {
		M_tm_plus = createKernel("M_tm_plus") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_m_tm_plus.cl";
		M_tm_minus = createKernel("M_tm_minus") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_m_tm_minus.cl";
	} else if(get_parameters()->get_fermact() == CLOVER) {
		throw Print_Error_Message("no kernels for CLOVER-discretization implemented yet, aborting... ", __FILE__, __LINE__);
	} else {
		throw Print_Error_Message("there was a problem with which fermion-discretization to use, aborting... ", __FILE__, __LINE__);
	}

	gamma5 = createKernel("gamma5") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_gamma5.cl";


	//Kernels needed if eoprec is used
	if(get_parameters()->get_use_eo() == true) {
		if(get_parameters()->get_fermact() == TWISTEDMASS) {
			M_tm_sitediagonal = createKernel("M_tm_sitediagonal") << basic_fermion_code << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m.cl";
			M_tm_inverse_sitediagonal = createKernel("M_tm_inverse_sitediagonal") << basic_fermion_code << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m.cl";
			M_tm_sitediagonal_minus = createKernel("M_tm_sitediagonal_minus") << basic_fermion_code << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m.cl";
			M_tm_inverse_sitediagonal_minus = createKernel("M_tm_inverse_sitediagonal_minus") << basic_fermion_code << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m.cl";
		}
		dslash_eoprec = createKernel("dslash_eoprec") << basic_fermion_code << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_dslash.cl";
		gamma5_eoprec = createKernel("gamma5_eoprec") << basic_fermion_code << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "fermionmatrix_eo_gamma5.cl";
	}
	return;
}

void Opencl_Module_Fermions::clear_kernels()
{
	Opencl_Module_Spinors::clear_kernels();

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

	if(get_parameters()->get_use_eo()) {
		clerr = clReleaseKernel(dslash_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(M_tm_sitediagonal);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(M_tm_inverse_sitediagonal);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}


	return;
}

void Opencl_Module_Fermions::clear_buffers()
{
	Opencl_Module_Spinors::clear_buffers();

	cl_uint clerr = CL_SUCCESS;


	clerr = clReleaseMemObject(clmem_inout);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_source);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);

	if(get_parameters()->get_use_eo()) {
		clerr = clReleaseMemObject(clmem_inout_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(clmem_source_even);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(clmem_source_odd);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(clmem_rn_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(clmem_rhat_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(clmem_v_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(clmem_p_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(clmem_s_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(clmem_t_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(clmem_aux_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clerr = clReleaseMemObject(clmem_tmp_eoprec_1);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		if(get_parameters()->get_fermact() == TWISTEDMASS) {
			clerr = clReleaseMemObject(clmem_tmp_eoprec_2);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		}
	}

	clerr = clReleaseMemObject(clmem_rho);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_rho_next);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_alpha);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_omega);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_beta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_tmp1);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_tmp2);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_one);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_minusone);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(clmem_resid);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);

	if(clmem_trueresid) {
		clerr = clReleaseMemObject(clmem_trueresid);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	}

	return;
}


void Opencl_Module_Fermions::get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups)
{
	Opencl_Module_Spinors::get_work_sizes(kernel, dev_type, ls, gs, num_groups);

	return;
}


//compound fermionmatrix-functions without eoprec
void Opencl_Module_Fermions::M(cl_mem in, cl_mem out, cl_mem gf)
{

	if(get_parameters()->get_fermact() == WILSON) {
		//in the pure Wilson case there is just one fermionmatrix
		M_wilson_device(in, out, gf);
	} else if(get_parameters()->get_fermact() == TWISTEDMASS) {
		M_tm_plus_device(in, out, gf);
	}
}

void Opencl_Module_Fermions::Qplus(cl_mem in, cl_mem out, cl_mem gf)
{
	if(get_parameters()->get_fermact() == WILSON) {
		//in the pure Wilson case there is just one fermionmatrix
		M_wilson_device(in, out, gf);
	} else if(get_parameters()->get_fermact() == TWISTEDMASS) {
		M_tm_plus_device(in, out, gf);
	}
	gamma5_device(out);
}

void Opencl_Module_Fermions::Qminus(cl_mem in, cl_mem out, cl_mem gf)
{
	if(get_parameters()->get_fermact() == WILSON) {
		//in the pure Wilson case there is just one fermionmatrix
		M_wilson_device(in, out, gf);
	} else if(get_parameters()->get_fermact() == TWISTEDMASS) {
		M_tm_minus_device(in, out, gf);
	}
	gamma5_device(out);
}

void Opencl_Module_Fermions::QplusQminus(cl_mem in, cl_mem out, cl_mem gf)
{
	/** @todo one could save one field here if an additional copying would be included in the end...
	 * or the field should be created in here, local */
	Qminus(in, clmem_tmp, gf);
	Qplus(clmem_tmp, out, gf);
}

//explicit fermionmatrix-kernel calling functions
void Opencl_Module_Fermions::M_wilson_device(cl_mem in, cl_mem out, cl_mem gf)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_wilson, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_wilson, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_wilson, 1, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_wilson, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( M_wilson, gs2, ls2);
}

void Opencl_Module_Fermions::M_tm_plus_device(cl_mem in, cl_mem out, cl_mem gf)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_tm_plus, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_plus, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_plus, 1, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_plus, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( M_tm_plus, gs2, ls2);
}

void Opencl_Module_Fermions::M_tm_minus_device(cl_mem in, cl_mem out, cl_mem gf)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_tm_minus, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_minus, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_minus, 1, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_minus, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( M_tm_minus, gs2, ls2);
}

void Opencl_Module_Fermions::gamma5_device(cl_mem inout)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gamma5, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gamma5, 0, sizeof(cl_mem), &inout);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(gamma5 , gs2, ls2);
}

//compound fermionmatrix-functions with eoprec
void Opencl_Module_Fermions::Aee(cl_mem in, cl_mem out, cl_mem gf)
{
	int even = EVEN;
	int odd = ODD;

	/**
	 * This is the even-odd preconditioned fermion matrix with the
	 * non-trivial inversion on the even sites (see DeGran/DeTar p. 174).
	 * If one has fermionmatrix
	 *  M = R + D,
	 * then Aee is:
	 * Aee = R_e - D_eo R_o_inv D_oe
	 */
	if(get_parameters()->get_fermact() == WILSON) {
		//in this case, the diagonal matrix is just 1 and falls away.
		dslash_eoprec_device(in, clmem_tmp_eoprec_1, gf, odd);
		dslash_eoprec_device(clmem_tmp_eoprec_1, out, gf, even);
		saxpy_eoprec_device(out, in, clmem_one, out);
	} else if(get_parameters()->get_fermact() == TWISTEDMASS) {
	  dslash_eoprec_device(in, clmem_tmp_eoprec_1, gf, odd);
	  	M_tm_inverse_sitediagonal_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2);
		dslash_eoprec_device(clmem_tmp_eoprec_2, out, gf, even);
		M_tm_sitediagonal_device(in, clmem_tmp_eoprec_1);
				saxpy_eoprec_device(out, clmem_tmp_eoprec_1, clmem_one, out);

		//without dslash
		//M_tm_inverse_sitediagonal_device(in, clmem_tmp_eoprec_1);
	        //M_tm_sitediagonal_device(in, out);
		//		saxpy_eoprec_device(out, clmem_tmp_eoprec_1, clmem_one, out);
		//saxpy_eoprec_device(clmem_tmp_eoprec_1, out , clmem_one, out); 

		//apply only one dslash
	  //	  dslash_eoprec_device(in, out, gf, odd);


		//M_tm_inverse_sitediagonal_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2);
	        //M_tm_sitediagonal_device(in, out);
		//saxpy_eoprec_device(clmem_tmp_eoprec_2, out , clmem_one, out); 

		//one more
	  //  dslash_eoprec_device(in, clmem_tmp_eoprec_1, gf, odd);
	  //	M_tm_inverse_sitediagonal_device(clmem_tmp_eoprec_1, out);

		//two more
		//dslash_eoprec_device(in, clmem_tmp_eoprec_1, gf, odd);
		//M_tm_inverse_sitediagonal_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2);
		//dslash_eoprec_device(clmem_tmp_eoprec_2, out, gf, even);

		//three more
		//dslash_eoprec_device(in, clmem_tmp_eoprec_1, gf, odd);
		//M_tm_inverse_sitediagonal_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2);
		//dslash_eoprec_device(clmem_tmp_eoprec_2, out, gf, even);

		//four more
		//dslash_eoprec_device(in, clmem_tmp_eoprec_1, gf, odd);
		//dslash_eoprec_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2, gf, even);
		//M_tm_inverse_sitediagonal_device(clmem_tmp_eoprec_2, out);

		//five more
		//dslash_eoprec_device(in, clmem_tmp_eoprec_1, gf, odd);
		//dslash_eoprec_device(clmem_tmp_eoprec_1, out, gf, even);

	}
}

/**
 *  This is the equivalent of the above function, but for the lower
 *  flavour, which essentially means mu -> -mu in the tm-case and
 *  no changes in the wilson case.
 */
void Opencl_Module_Fermions::Aee_minus(cl_mem in, cl_mem out, cl_mem gf)
{
	int even = EVEN;
	int odd = ODD;

	/**
	 * This is the even-odd preconditioned fermion matrix with the
	 * non-trivial inversion on the even sites (see DeGran/DeTar p. 174).
	 * If one has fermionmatrix
	 *  M = R + D,
	 * then Aee is:
	 * Aee = R_e - D_eo R_o_inv D_oe
	 * and Aee_minus is:
	 * Aee_minus = R_e(-mu) - D_eo R_o(-mu)_inv D_oe
	 */
	if(get_parameters()->get_fermact() == WILSON) {
		//in this case, the diagonal matrix is just 1 and falls away.
		dslash_eoprec_device(in, clmem_tmp_eoprec_1, gf, odd);
		dslash_eoprec_device(clmem_tmp_eoprec_1, out, gf, even);
		saxpy_eoprec_device(out, in, clmem_one, out);
	} else if(get_parameters()->get_fermact() == TWISTEDMASS) {
		dslash_eoprec_device(in, clmem_tmp_eoprec_1, gf, odd);
		M_tm_inverse_sitediagonal_minus_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2);
		dslash_eoprec_device(clmem_tmp_eoprec_2, out, gf, even);
		M_tm_sitediagonal_minus_device(in, clmem_tmp_eoprec_1);
		saxpy_eoprec_device(out, clmem_tmp_eoprec_1, clmem_one, out);
	}
}

void Opencl_Module_Fermions::Qplus_eoprec(cl_mem in, cl_mem out, cl_mem gf)
{
	Aee(in, out, gf);
	gamma5_eoprec_device(out);
	return;
}

void Opencl_Module_Fermions::Qminus_eoprec(cl_mem in, cl_mem out, cl_mem gf)
{
	Aee_minus(in, out, gf);
	gamma5_eoprec_device(out);
	return;
}

void Opencl_Module_Fermions::QplusQminus_eoprec(cl_mem in, cl_mem out, cl_mem gf)
{
	Qminus_eoprec(in, clmem_tmp_eoprec_1, gf);
	Qplus_eoprec(clmem_tmp_eoprec_1, out, gf);
	return;
}

//explicit eoprec fermionmatrix functions
void Opencl_Module_Fermions::gamma5_eoprec_device(cl_mem inout)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gamma5_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gamma5_eoprec, 0, sizeof(cl_mem), &inout);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( gamma5_eoprec, gs2, ls2);
}

void Opencl_Module_Fermions::dslash_eoprec_device(cl_mem in, cl_mem out, cl_mem gf, int evenodd)
{
	int eo = evenodd;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(dslash_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments

	int clerr = clSetKernelArg(dslash_eoprec, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_eoprec, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_eoprec, 2, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_eoprec, 3, sizeof(int), &eo);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(dslash_eoprec , gs2, ls2);
}

void Opencl_Module_Fermions::M_tm_inverse_sitediagonal_device(cl_mem in, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_tm_inverse_sitediagonal, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_inverse_sitediagonal, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_inverse_sitediagonal, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( M_tm_inverse_sitediagonal, gs2, ls2);
}

void Opencl_Module_Fermions::M_tm_sitediagonal_device(cl_mem in, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_tm_sitediagonal, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_sitediagonal, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(M_tm_sitediagonal , gs2, ls2);
}

void Opencl_Module_Fermions::M_tm_inverse_sitediagonal_minus_device(cl_mem in, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_tm_inverse_sitediagonal_minus, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_inverse_sitediagonal_minus, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_inverse_sitediagonal_minus, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( M_tm_inverse_sitediagonal_minus, gs2, ls2);
}

void Opencl_Module_Fermions::M_tm_sitediagonal_minus_device(cl_mem in, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(M_tm_sitediagonal_minus, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_sitediagonal_minus, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal_minus, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(M_tm_sitediagonal_minus , gs2, ls2);
}

int Opencl_Module_Fermions::bicgstab(const Matrix_Function & f, cl_mem inout, cl_mem source, cl_mem gf, hmc_float prec)
{
	int debug = 0;
	int old = 1;
	if(old == 0) {
		if(debug) cout << "debug-output at bicgstab is activated" << endl;
		int cgmax = get_parameters()->get_cgmax();
		//CP: these have to be on the host
		hmc_float resid;
		hmc_float trueresid;

		for(int iter = 0; iter < cgmax; iter++) {
			if(iter % get_parameters()->get_iter_refresh() == 0) {
				set_zero_spinorfield_device(clmem_v);
				set_zero_spinorfield_device(clmem_p);
				f(inout, clmem_rn, gf);

				saxpy_device(clmem_rn, source, clmem_one, clmem_rn);
				copy_buffer_on_device(clmem_rn, clmem_rhat, sizeof(spinor) * get_parameters()->get_spinorfieldsize());

				copy_buffer_on_device(clmem_one, clmem_alpha, sizeof(hmc_complex));
				copy_buffer_on_device(clmem_one, clmem_omega, sizeof(hmc_complex));
				copy_buffer_on_device(clmem_one, clmem_rho, sizeof(hmc_complex));

				//CP: calc initial residuum for output, this is not needed for the algorithm!!
//      set_float_to_global_squarenorm_device(clmem_rn, clmem_resid, local_work_size, global_work_size);
//			get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));
//      cout << "initial residuum at iter " << iter << " is: " << scientific << resid << endl;
				//printf("initial residuum at iter %i is %.40e\n", iter, resid);
			}

////////////////////////////////////
//collect all variables if debug is enabled
			if(debug) {
				hmc_complex omega;
				hmc_complex rho;
				hmc_complex rho_next;
				hmc_complex tmp1;
				hmc_complex tmp2;
				hmc_complex beta;
				hmc_complex alpha;

				get_buffer_from_device(clmem_omega, &omega, sizeof(hmc_complex));
				get_buffer_from_device(clmem_rho, &rho, sizeof(hmc_complex));
				get_buffer_from_device(clmem_rho_next, &rho_next, sizeof(hmc_complex));
				get_buffer_from_device(clmem_tmp1, &tmp1, sizeof(hmc_complex));
				get_buffer_from_device(clmem_tmp2, &tmp2, sizeof(hmc_complex));
				get_buffer_from_device(clmem_beta, &beta, sizeof(hmc_complex));
				get_buffer_from_device(clmem_alpha, &alpha, sizeof(hmc_complex));

				cout << "debug output at start: " << endl;
				cout << " rho: " << rho.re << "  " <<  rho.im << endl;
				cout << " alpha: " << alpha.re << "  " <<  alpha.im << endl;
				cout << " omega: " << omega.re << "  " <<  omega.im << endl;
			}
////////////////////////////////////

			set_complex_to_scalar_product_device(clmem_rhat, clmem_rn, clmem_rho_next);

			if(debug) {
				hmc_complex rho_next;
				get_buffer_from_device(clmem_rho_next, &rho_next, sizeof(hmc_complex));
				cout << "(rhat, rn): " << rho_next.re << "  " <<  rho_next.im << endl;
			}
			set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1);
			copy_buffer_on_device(clmem_rho_next, clmem_rho, sizeof(hmc_complex));
			set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2);
			set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta);

			set_complex_to_product_device(clmem_beta, clmem_omega, clmem_tmp1);
			set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2);
			saxsbypz_device(clmem_p, clmem_v, clmem_rn, clmem_beta, clmem_tmp2, clmem_p);

			if(debug) {
				hmc_complex rho_next;
				set_complex_to_scalar_product_device(clmem_p, clmem_p, clmem_tmp2);
				get_buffer_from_device(clmem_tmp2, &rho_next, sizeof(hmc_complex));
				cout << "(p,p): " << rho_next.re << "  " <<  rho_next.im << endl;
				set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2);
			}

			f(clmem_p, clmem_v, gf);

////////////////////////////////////
//collect all variables if debug is enabled
			if(debug) {
				hmc_complex omega;
				hmc_complex rho;
				hmc_complex rho_next;
				hmc_complex tmp1;
				hmc_complex tmp2;
				hmc_complex beta;
				hmc_complex alpha;

				get_buffer_from_device(clmem_omega, &omega, sizeof(hmc_complex));
				get_buffer_from_device(clmem_rho, &rho, sizeof(hmc_complex));
				get_buffer_from_device(clmem_rho_next, &rho_next, sizeof(hmc_complex));
				get_buffer_from_device(clmem_tmp1, &tmp1, sizeof(hmc_complex));
				get_buffer_from_device(clmem_tmp2, &tmp2, sizeof(hmc_complex));
				get_buffer_from_device(clmem_beta, &beta, sizeof(hmc_complex));
				get_buffer_from_device(clmem_alpha, &alpha, sizeof(hmc_complex));

				cout << "debug output after first half: " << endl;
				cout << " rho: " << rho.re << "  " <<  rho.im << endl;
				cout << " rho_next: " << rho_next.re << "  " <<  rho_next.im << endl;
				cout << " beta: " << beta.re << "  " <<  beta.im << endl;
				cout << " alpha: " << alpha.re << "  " <<  alpha.im << endl;
				cout << " omega: " << omega.re << "  " <<  omega.im << endl;
				cout << " tmp1: " << tmp1.re << "  " <<  tmp1.im << endl;
				cout << " tmp2: " << tmp2.re << "  " <<  tmp2.im << endl;
			}
////////////////////////////////////


			set_complex_to_scalar_product_device(clmem_rhat, clmem_v, clmem_tmp1);

			if(debug) {
				hmc_complex rho_next;
				get_buffer_from_device(clmem_tmp1, &rho_next, sizeof(hmc_complex));
				cout << "(rhat, v): " << rho_next.re << "  " <<  rho_next.im << endl;
			}

			set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha);
			saxpy_device(clmem_v, clmem_rn, clmem_alpha, clmem_s);

			f(clmem_s, clmem_t, gf);

			set_complex_to_scalar_product_device(clmem_t, clmem_s, clmem_tmp1);
			//!!CP: this can also be global_squarenorm
			set_complex_to_scalar_product_device(clmem_t, clmem_t, clmem_tmp2);
			set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega);

			saxpy_device(clmem_t, clmem_s, clmem_omega, clmem_rn);

			saxsbypz_device(clmem_p, clmem_s, inout, clmem_alpha, clmem_omega, inout);

			if(debug) {
				hmc_complex rho_next;
				set_complex_to_scalar_product_device(inout, inout, clmem_omega);
				get_buffer_from_device(clmem_omega, &rho_next, sizeof(hmc_complex));
				cout << "(inout,inout): " << rho_next.re << "  " <<  rho_next.im << endl;
				set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega);
			}

////////////////////////////////////
//collect all variables if debug is enabled
			if(debug) {
				hmc_complex omega;
				hmc_complex rho;
				hmc_complex rho_next;
				hmc_complex tmp1;
				hmc_complex tmp2;
				hmc_complex beta;
				hmc_complex alpha;

				get_buffer_from_device(clmem_omega, &omega, sizeof(hmc_complex));
				get_buffer_from_device(clmem_rho, &rho, sizeof(hmc_complex));
				get_buffer_from_device(clmem_rho_next, &rho_next, sizeof(hmc_complex));
				get_buffer_from_device(clmem_tmp1, &tmp1, sizeof(hmc_complex));
				get_buffer_from_device(clmem_tmp2, &tmp2, sizeof(hmc_complex));
				get_buffer_from_device(clmem_beta, &beta, sizeof(hmc_complex));
				get_buffer_from_device(clmem_alpha, &alpha, sizeof(hmc_complex));

				cout << "debug output after second half: " << endl;
				cout << " rho: " << rho.re << "  " <<  rho.im << endl;
				cout << " rho_next: " << rho_next.re << "  " <<  rho_next.im << endl;
				cout << " beta: " << beta.re << "  " <<  beta.im << endl;
				cout << " alpha: " << alpha.re << "  " <<  alpha.im << endl;
				cout << " omega: " << omega.re << "  " <<  omega.im << endl;

				cout << " tmp1: " << tmp1.re << "  " <<  tmp1.im << endl;
				cout << " tmp2: " << tmp2.re << "  " <<  tmp2.im << endl;
			}
////////////////////////////////////

			set_float_to_global_squarenorm_device(clmem_rn, clmem_resid);
			get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));

			logger.debug() << "resid: " << resid;

			if(resid < prec) {
				f(inout, clmem_aux, gf);
				saxpy_device(clmem_aux, source, clmem_one, clmem_aux);
				set_float_to_global_squarenorm_device(clmem_aux, clmem_trueresid);
				get_buffer_from_device(clmem_trueresid, &trueresid, sizeof(hmc_float));
//      cout << "\tsolver converged! residuum:\t" << resid << " is smaller than " << prec << endl;
//      cout << "\ttrueresiduum:\t" << trueresid << " has to be smaller than " << prec << endl;
				if(trueresid < prec)
					return iter;
				else {
//        cout << "trueresiduum not small enough" <<endl;
//        hmc_complex s_norm;
//        //borrow clmem_alpha for this
//        set_complex_to_scalar_product_device(clmem_s, clmem_s, clmem_alpha, localsize, globalsize);
//        get_buffer_from_device(clmem_alpha, &s_norm, sizeof(hmc_complex));
//        cout << "|s|^2: " << s_norm.re << "  " <<  s_norm.im << endl;
//        //reset value of alpha
//        set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha);
//        //check if |s|^2 is too small
//        if(s_norm.re < prec){
//          cout << "|s|^2 is too small to continue..." << endl;
// //           return;
//        }
				}
			} else {
//      printf("residuum at iter%i is:\t%.10e\n", iter, resid);//cout << "residuum:\t" << resid << endl;
			}
		}
		return -1;
	}
	//"save" version, with comments. this is called if "bicgstab_save" is choosen.
	else if (get_parameters()->get_use_bicgstab_save() == true) {
		for(int iter = 0; iter < get_parameters()->get_cgmax(); iter++) {
			if(iter % get_parameters()->get_iter_refresh() == 0) {
				set_zero_spinorfield_device(clmem_v);
				set_zero_spinorfield_device(clmem_p);
				//initial r_n
				f(inout, clmem_rn, gf);
				saxpy_device(clmem_rn, source, clmem_one, clmem_rn);
				//rhat = r_n
				copy_buffer_on_device(clmem_rn, clmem_rhat, get_parameters()->get_sf_buf_size());
				//set some constants to 1
				copy_buffer_on_device(clmem_one, clmem_alpha, sizeof(hmc_complex));
				copy_buffer_on_device(clmem_one, clmem_omega, sizeof(hmc_complex));
				copy_buffer_on_device(clmem_one, clmem_rho, sizeof(hmc_complex));
			}
			//rho_next = (rhat, rn)
			set_complex_to_scalar_product_device(clmem_rhat, clmem_rn, clmem_rho_next);
			//check if algorithm is stuck
			hmc_complex check;
			get_buffer_from_device(clmem_rho_next, &check, sizeof(hmc_complex));
			//if rho is too small the algorithm will get stuck and will never converge!!
			if(abs(check.re) < 1e-25 && abs(check.im) < 1e-25 ) {
				return -iter;
			}
			//tmp1 = rho_next/rho = (rhat, rn)/..
			set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1);
			//rho_next = rho
			copy_buffer_on_device(clmem_rho_next, clmem_rho, sizeof(hmc_complex));
			//tmp2 = alpha/omega = ...
			set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2);
			//beta = tmp1*tmp2
			set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta);

			//tmp1 = beta*omega
			set_complex_to_product_device(clmem_beta, clmem_omega, clmem_tmp1);
			//tmp2 = -tmp1
			set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2);
			//p = beta*p + tmp2*v + r_n = beta*p - beta*omega*v + r_n
			saxsbypz_device(clmem_p, clmem_v, clmem_rn, clmem_beta, clmem_tmp2, clmem_p);

			//v = A*p
			f(clmem_p, clmem_v, gf);
			//tmp1 = (rhat, v)
			set_complex_to_scalar_product_device(clmem_rhat, clmem_v, clmem_tmp1);
			//alpha = rho/tmp1 = (..)/(rhat, v)
			set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha);
			//s = - alpha * v - r_n
			saxpy_device(clmem_v, clmem_rn, clmem_alpha, clmem_s);
			//t = A s
			f(clmem_s, clmem_t, gf);
			//tmp1 = (t, s)
			set_complex_to_scalar_product_device(clmem_t, clmem_s, clmem_tmp1);
			//!!CP: this can also be global_squarenorm, but one needs a complex number here
			//tmp2 = (t,t)
			set_complex_to_scalar_product_device(clmem_t, clmem_t, clmem_tmp2);
			//omega = tmp1/tmp2 = (t,s)/(t,t)
			set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega);
			//r_n = - omega*t - s
			saxpy_device(clmem_t, clmem_s, clmem_omega, clmem_rn);
			//inout = alpha*p + omega * s + inout
			saxsbypz_device(clmem_p, clmem_s, inout, clmem_alpha, clmem_omega, inout);
			//resid = (rn,rn)
			set_float_to_global_squarenorm_device(clmem_rn, clmem_resid);
			hmc_float resid;
			get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));

//			cout << "resid at iter " << iter << " is: " << resid << endl;

			if(resid < prec) {
				//aux = A inout
				f(inout, clmem_aux, gf);
				//aux = -aux + source
				saxpy_device(clmem_aux, source, clmem_one, clmem_aux);
				//trueresid = (aux, aux)
				set_float_to_global_squarenorm_device(clmem_aux, clmem_trueresid);
				hmc_float trueresid;
				get_buffer_from_device(clmem_trueresid, &trueresid, sizeof(hmc_float));
//				cout << "\tsolver converged! residuum:\t" << resid << " is smaller than " << prec << endl;
//        cout << "\ttrueresiduum:\t" << trueresid << " has to be smaller than " << prec << endl;
				if(trueresid < prec)
					return iter;
				else {
//					cout << "trueresiduum not small enough" <<endl;
				}
			}
		}
		return -1;
	}
	//version with different structure than "save" one, similar to tmlqcd. This should be the default bicgstab.
	//  In particular this version does not perform the check if the "real" residuum is sufficiently small!
	else if (get_parameters()->get_use_bicgstab_save() != true) {
		for(int iter = 0; iter < get_parameters()->get_cgmax(); iter++) {
			if(iter % get_parameters()->get_iter_refresh() == 0) {
				//initial r_n, saved in p
				f(inout, clmem_rn, gf);
				saxpy_device(clmem_rn, source, clmem_one, clmem_p);
				//rhat = p
				copy_buffer_on_device(clmem_p, clmem_rhat, get_parameters()->get_sf_buf_size());
				//r_n = p
				copy_buffer_on_device(clmem_p, clmem_rn, get_parameters()->get_sf_buf_size());
				//rho = (rhat, rn)
				set_complex_to_scalar_product_device(clmem_rhat, clmem_rn, clmem_rho);
			}
			//resid = (rn,rn)
			set_float_to_global_squarenorm_device(clmem_rn, clmem_resid);
			hmc_float resid;
			get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));
			if(resid < prec) {
				return iter;
			}
			//v = A*p
			f(clmem_p, clmem_v, gf);
			//tmp1 = (rhat, v)
			set_complex_to_scalar_product_device(clmem_rhat, clmem_v, clmem_tmp1);
			//alpha = rho/tmp1 = (rhat, rn)/(rhat, v)
			set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha);
			//s = - alpha * v - r_n
			saxpy_device(clmem_v, clmem_rn, clmem_alpha, clmem_s);
			//t = A s
			f(clmem_s, clmem_t, gf);
			//tmp1 = (t, s)
			set_complex_to_scalar_product_device(clmem_t, clmem_s, clmem_tmp1);
			//!!CP: this can also be global_squarenorm, but one needs a complex number here
			//tmp2 = (t,t)
			set_complex_to_scalar_product_device(clmem_t, clmem_t, clmem_tmp2);
			//omega = tmp1/tmp2 = (t,s)/(t,t)
			set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega);
			//inout = alpha*p + omega * s + inout
			saxsbypz_device(clmem_p, clmem_s, inout, clmem_alpha, clmem_omega, inout);
			//r_n = - omega*t - s
			saxpy_device(clmem_t, clmem_s, clmem_omega, clmem_rn);
			//rho_next = (rhat, rn)
			set_complex_to_scalar_product_device(clmem_rhat, clmem_rn, clmem_rho_next);
			//check if algorithm is stuck
			hmc_complex check;
			get_buffer_from_device(clmem_rho_next, &check, sizeof(hmc_complex));
			//if rho is too small the algorithm will get stuck and will never converge!!
			if(abs(check.re) < 1e-25 && abs(check.im) < 1e-25 ) {
				return -iter;
			}
			//tmp1 = rho_next/rho = (rhat, rn)/..
			set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1);
			//tmp2 = alpha/omega = ...
			set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2);
			//beta = tmp1*tmp2 = alpha*rho_next / (omega*rho)
			set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta);
			//tmp1 = beta*omega = alpha* rho_next / rho
			set_complex_to_product_device(clmem_beta, clmem_omega, clmem_tmp1);
			//tmp2 = -tmp1
			set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2);
			//p = beta*p + tmp2*v + r_n = beta*p - beta*omega*v + r_n
			saxsbypz_device(clmem_p, clmem_v, clmem_rn, clmem_beta, clmem_tmp2, clmem_p);
			//rho_next = rho
			copy_buffer_on_device(clmem_rho_next, clmem_rho, sizeof(hmc_complex));
		}
		return -1;
	}
	//this return cannot be reached and is inserted to remove a warning
	return 0;
}

int Opencl_Module_Fermions::bicgstab_eoprec(const Matrix_Function & f, cl_mem inout, cl_mem source, cl_mem gf, hmc_float prec)
{
	//"save" version, with comments. this is called if "bicgstab_save" is choosen.
	if (get_parameters()->get_use_bicgstab_save()) {
		klepsydra::Monotonic timer;
		if(logger.beInfo()) {
			cl_event start_event;
			clEnqueueMarker(get_queue(), &start_event);
			clSetEventCallback(start_event, CL_COMPLETE, resetTimerOnComplete, &timer);
			clReleaseEvent(start_event);
		}
		//CP: these have to be on the host
		hmc_float resid;
		hmc_float trueresid;
		unsigned retests = 0;
		int cgmax = get_parameters()->get_cgmax();
		for(int iter = 0; iter < cgmax; iter++) {
			if(iter % get_parameters()->get_iter_refresh() == 0) {
				set_zero_spinorfield_eoprec_device(clmem_v_eoprec);
				set_zero_spinorfield_eoprec_device(clmem_p_eoprec);

				f(inout, clmem_rn_eoprec, gf);

				saxpy_eoprec_device(clmem_rn_eoprec, source, clmem_one, clmem_rn_eoprec);

				//set_float_to_global_squarenorm_eoprec_device(clmem_rn_eoprec, clmem_resid);
				//get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));
				//cout << "init residuum:\t" << resid << endl;

				copy_buffer_on_device(clmem_rn_eoprec, clmem_rhat_eoprec, get_parameters()->get_eo_sf_buf_size());

				copy_buffer_on_device(clmem_one, clmem_alpha, sizeof(hmc_complex));
				copy_buffer_on_device(clmem_one, clmem_omega, sizeof(hmc_complex));
				copy_buffer_on_device(clmem_one, clmem_rho, sizeof(hmc_complex));
			}
			set_complex_to_scalar_product_eoprec_device(clmem_rhat_eoprec, clmem_rn_eoprec, clmem_rho_next);
			//check if algorithm is stuck
			hmc_complex check;
			get_buffer_from_device(clmem_rho_next, &check, sizeof(hmc_complex));
			if(abs(check.re) < 1e-25 && abs(check.im) < 1e-25 ) {
				return -iter;
			}
			set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1);
			copy_buffer_on_device(clmem_rho_next, clmem_rho, sizeof(hmc_complex));
			set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2);
			set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta);

			set_complex_to_product_device(clmem_beta, clmem_omega, clmem_tmp1);
			set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2);
			saxsbypz_eoprec_device(clmem_p_eoprec, clmem_v_eoprec, clmem_rn_eoprec, clmem_beta, clmem_tmp2, clmem_p_eoprec);

			f(clmem_p_eoprec, clmem_v_eoprec, gf);

			set_complex_to_scalar_product_eoprec_device(clmem_rhat_eoprec, clmem_v_eoprec, clmem_tmp1);
			set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha);

			saxpy_eoprec_device(clmem_v_eoprec, clmem_rn_eoprec, clmem_alpha, clmem_s_eoprec);

			f(clmem_s_eoprec, clmem_t_eoprec, gf);

			set_complex_to_scalar_product_eoprec_device(clmem_t_eoprec, clmem_s_eoprec, clmem_tmp1);
			//!!CP: can this also be global_squarenorm??
			set_complex_to_scalar_product_eoprec_device(clmem_t_eoprec, clmem_t_eoprec, clmem_tmp2);
			set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega);

			saxpy_eoprec_device(clmem_t_eoprec, clmem_s_eoprec, clmem_omega, clmem_rn_eoprec);

			saxsbypz_eoprec_device(clmem_p_eoprec, clmem_s_eoprec, inout, clmem_alpha, clmem_omega, inout);

			set_float_to_global_squarenorm_eoprec_device(clmem_rn_eoprec, clmem_resid);
			get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));

			logger.debug() << "resid: " << resid;

			if(resid < prec) {
				++retests;

				f(inout, clmem_aux_eoprec, gf);
				saxpy_eoprec_device(clmem_aux_eoprec, source, clmem_one, clmem_aux_eoprec);

				set_float_to_global_squarenorm_eoprec_device(clmem_aux_eoprec, clmem_trueresid);
				get_buffer_from_device(clmem_trueresid, &trueresid, sizeof(hmc_float));
				//cout << "residuum:\t" << resid << "\ttrueresiduum:\t" << trueresid << endl;
				if(trueresid < prec) {
					// report on performance
					if(logger.beInfo()) {
						// we are always synchroneous here, as we had to recieve the residium from the device
						uint64_t duration = timer.getTime();

						// calculate flops
						unsigned refreshs = iter / get_parameters()->get_iter_refresh() + 1;
						cl_ulong mf_flops = f.get_Flops();

						cl_ulong total_flops = 4 * get_flop_size("scalar_product_eoprec") + 4 * get_flop_size("ratio") + 3 * get_flop_size("product") + 2 * get_flop_size("saxsbypz_eoprec") + 2 * mf_flops + 2 * get_flop_size("saxpy_eoprec") + get_flop_size("global_squarenorm_eoprec");
						total_flops *= iter;

						total_flops += refreshs * (mf_flops + get_flop_size("saxpy_eoprec"));

						total_flops += retests * (mf_flops + get_flop_size("saxpy_eoprec") + get_flop_size("global_squarenorm_eoprec"));

						// report performanc
						logger.info() << "BiCGstab_save completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f) << " Gflops.";
					}

					// we are done here
					return iter;
				} else {
//					cout << "trueresiduum not small enough" <<endl;
				}
			}
		}
		return -1;
	} else {
		//version with different structure than "save" one, similar to tmlqcd. This should be the default bicgstab.
		//  In particular this version does not perform the check if the "real" residuum is sufficiently small!
		klepsydra::Monotonic timer;
		if(logger.beInfo()) {
			cl_event start_event;
			clEnqueueMarker(get_queue(), &start_event);
			clSetEventCallback(start_event, CL_COMPLETE, resetTimerOnComplete, &timer);
			clReleaseEvent(start_event);
		}
		for(int iter = 0; iter < get_parameters()->get_cgmax(); iter++) {
			if(iter % get_parameters()->get_iter_refresh() == 0) {
				//initial r_n, saved in p
				f(inout, clmem_rn_eoprec, gf);
				saxpy_eoprec_device(clmem_rn_eoprec, source, clmem_one, clmem_p_eoprec);
				//rhat = p
				copy_buffer_on_device(clmem_p_eoprec, clmem_rhat_eoprec, get_parameters()->get_eo_sf_buf_size());
				//r_n = p
				copy_buffer_on_device(clmem_p_eoprec, clmem_rn_eoprec, get_parameters()->get_eo_sf_buf_size());
				//rho = (rhat, rn)
				set_complex_to_scalar_product_eoprec_device(clmem_rhat_eoprec, clmem_rn_eoprec, clmem_rho);
			}
			//resid = (rn,rn)
			set_float_to_global_squarenorm_eoprec_device(clmem_rn_eoprec, clmem_resid);
			hmc_float resid;
			get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));

			logger.debug() << "resid: " << resid;

			if(resid < prec) {
				// report on performance
				if(logger.beInfo()) {
					// we are always synchroneous here, as we had to recieve the residium from the device
					uint64_t duration = timer.getTime();

					// calculate flops
					unsigned refreshs = iter / get_parameters()->get_iter_refresh() + 1;
					cl_ulong mf_flops = f.get_Flops();

					cl_ulong total_flops = get_flop_size("global_squarenorm_eoprec") + 2 * mf_flops + 4 * get_flop_size("scalar_product_eoprec") + 4 * get_flop_size("ratio") + 2 * get_flop_size("saxpy_eoprec") + 2 * get_flop_size("saxsbypz_eoprec") + 3 * get_flop_size("product");
					total_flops *= iter;

					total_flops += refreshs * (mf_flops + get_flop_size("saxpy_eoprec") + get_flop_size("scalar_product_eoprec"));

					// report performanc
					logger.info() << "BiCGstab completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f) << " Gflops.";
				}

				// we are done here
				return iter;
			}
			//v = A*p
			f(clmem_p_eoprec, clmem_v_eoprec, gf);
			//tmp1 = (rhat, v)
			set_complex_to_scalar_product_eoprec_device(clmem_rhat_eoprec, clmem_v_eoprec, clmem_tmp1);
			//alpha = rho/tmp1 = (rhat, rn)/(rhat, v)
			set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha);
			//s = - alpha * v - r_n
			saxpy_eoprec_device(clmem_v_eoprec, clmem_rn_eoprec, clmem_alpha, clmem_s_eoprec);
			//t = A s
			f(clmem_s_eoprec, clmem_t_eoprec, gf);
			//tmp1 = (t, s)
			set_complex_to_scalar_product_eoprec_device(clmem_t_eoprec, clmem_s_eoprec, clmem_tmp1);
			//!!CP: this can also be global_squarenorm, but one needs a complex number here
			//tmp2 = (t,t)
			set_complex_to_scalar_product_eoprec_device(clmem_t_eoprec, clmem_t_eoprec, clmem_tmp2);
			//omega = tmp1/tmp2 = (t,s)/(t,t)
			set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega);
			//inout = alpha*p + omega * s + inout
			saxsbypz_eoprec_device(clmem_p_eoprec, clmem_s_eoprec, inout, clmem_alpha, clmem_omega, inout);
			//r_n = - omega*t - s
			saxpy_eoprec_device(clmem_t_eoprec, clmem_s_eoprec, clmem_omega, clmem_rn_eoprec);
			//rho_next = (rhat, rn)
			set_complex_to_scalar_product_eoprec_device(clmem_rhat_eoprec, clmem_rn_eoprec, clmem_rho_next);
			//check if algorithm is stuck
			hmc_complex check;
			get_buffer_from_device(clmem_rho_next, &check, sizeof(hmc_complex));
			if(abs(check.re) < 1e-25 && abs(check.im) < 1e-25 ) {
				return -iter;
			}

			//tmp1 = rho_next/rho = (rhat, rn)/..
			set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1);
			//tmp2 = alpha/omega = ...
			set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2);
			//beta = tmp1*tmp2 = alpha*rho_next / (omega*rho)
			set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta);
			//tmp1 = beta*omega = alpha* rho_next / rho
			set_complex_to_product_device(clmem_beta, clmem_omega, clmem_tmp1);
			//tmp2 = -tmp1
			set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2);
			//p = beta*p + tmp2*v + r_n = beta*p - beta*omega*v + r_n
			saxsbypz_eoprec_device(clmem_p_eoprec, clmem_v_eoprec, clmem_rn_eoprec, clmem_beta, clmem_tmp2, clmem_p_eoprec);
			//rho_next = rho
			copy_buffer_on_device(clmem_rho_next, clmem_rho, sizeof(hmc_complex));
		}
		return -1;
	}
	//this return cannot be reached and is inserted to remove a warning
	return 0;
}

int Opencl_Module_Fermions::cg(const Matrix_Function & f, cl_mem inout, cl_mem source, cl_mem gf, hmc_float prec)
{
	//CP: here I do not use clmem_rnhat anymore and saved one scalar_product (omega)
	//NOTE: here, most of the complex numbers may also be just hmc_floats. However, for this one would need some add. functions...
	for(int iter = 0; iter < get_parameters()->get_cgmax(); iter ++) {
		if(iter % get_parameters()->get_iter_refresh() == 0) {
			//rn = A*inout
			f(inout, clmem_rn, gf);
			//rn = source - A*inout
			saxpy_device(clmem_rn, source, clmem_one, clmem_rn);
			//p = rn
			copy_buffer_on_device(clmem_rn, clmem_p, get_parameters()->get_sf_buf_size());
			//omega = (rn,rn)
			set_complex_to_scalar_product_device(clmem_rn, clmem_rn, clmem_omega);
		} else {
			//update omega
			copy_buffer_on_device(clmem_rho_next, clmem_omega, sizeof(hmc_complex));
		}
		//v = A pn
		f(clmem_p, clmem_v, gf);
		//alpha = (rn, rn)/(pn, Apn) --> alpha = omega/rho
		set_complex_to_scalar_product_device(clmem_p, clmem_v, clmem_rho);
		set_complex_to_ratio_device(clmem_omega, clmem_rho, clmem_alpha);
		set_complex_to_product_device(clmem_alpha, clmem_minusone, clmem_tmp1);

		//xn+1 = xn + alpha*p = xn - tmp1*p = xn - (-tmp1)*p
		saxpy_device(clmem_p, inout, clmem_tmp1, inout);
		//rn+1 = rn - alpha*v -> rhat
		saxpy_device(clmem_v, clmem_rn, clmem_alpha, clmem_rn);

		//calc residuum
		//NOTE: for beta one needs a complex number at the moment, therefore, this is done with "rho_next" instead of "resid"
		set_complex_to_scalar_product_device(clmem_rn, clmem_rn, clmem_rho_next);
		hmc_float resid;
		get_buffer_from_device(clmem_rho_next, &resid, sizeof(hmc_float));
		//this is the orig. call
		//set_float_to_global_squarenorm_device(clmem_rn, clmem_resid);
		//get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));

		logger.debug() << "resid: " << resid;

		if(resid < prec)
			return iter;

		//beta = (rn+1, rn+1)/(rn, rn) --> alpha = rho_next/omega
		set_complex_to_scalar_product_device(clmem_rn, clmem_rn, clmem_rho_next);
		set_complex_to_ratio_device(clmem_rho_next, clmem_omega, clmem_beta);

		//pn+1 = rn+1 + beta*pn
		set_complex_to_product_device(clmem_beta, clmem_minusone, clmem_tmp2);
		saxpy_device(clmem_p, clmem_rn, clmem_tmp2, clmem_p);
	}
	return -1;
}

int Opencl_Module_Fermions::cg_eoprec(const Matrix_Function & f, cl_mem inout, cl_mem source, cl_mem gf, hmc_float prec)
{
	//this corresponds to the above function
	//NOTE: here, most of the complex numbers may also be just hmc_floats. However, for this one would need some add. functions...
	klepsydra::Monotonic timer;
	if(logger.beInfo()) {
		cl_event start_event;
		clEnqueueMarker(get_queue(), &start_event);
		clSetEventCallback(start_event, CL_COMPLETE, resetTimerOnComplete, &timer);
		clReleaseEvent(start_event);
	}
	for(int iter = 0; iter < get_parameters()->get_cgmax(); iter ++) {
		if(iter % get_parameters()->get_iter_refresh() == 0) {
			//rn = A*inout
			f(inout, clmem_rn_eoprec, gf);
			//rn = source - A*inout
			saxpy_eoprec_device(clmem_rn_eoprec, source, clmem_one, clmem_rn_eoprec);
			//p = rn
			copy_buffer_on_device(clmem_rn_eoprec, clmem_p_eoprec, get_parameters()->get_eo_sf_buf_size());
			//omega = (rn,rn)
			set_complex_to_scalar_product_eoprec_device(clmem_rn_eoprec, clmem_rn_eoprec, clmem_omega);
		} else {
			//update omega
			copy_buffer_on_device(clmem_rho_next, clmem_omega, sizeof(hmc_complex));
		}
		//v = A pn
		f(clmem_p_eoprec, clmem_v_eoprec, gf);
		//alpha = (rn, rn)/(pn, Apn) --> alpha = omega/rho
		set_complex_to_scalar_product_eoprec_device(clmem_p_eoprec, clmem_v_eoprec, clmem_rho);
		set_complex_to_ratio_device(clmem_omega, clmem_rho, clmem_alpha);
		set_complex_to_product_device(clmem_alpha, clmem_minusone, clmem_tmp1);

		//xn+1 = xn + alpha*p = xn - tmp1*p = xn - (-tmp1)*p
		saxpy_eoprec_device(clmem_p_eoprec, inout, clmem_tmp1, inout);
		//rn+1 = rn - alpha*v -> rhat
		saxpy_eoprec_device(clmem_v_eoprec, clmem_rn_eoprec, clmem_alpha, clmem_rn_eoprec);

		//calc residuum
		//NOTE: for beta one needs a complex number at the moment, therefore, this is done with "rho_next" instead of "resid"
		set_complex_to_scalar_product_eoprec_device(clmem_rn_eoprec, clmem_rn_eoprec, clmem_rho_next);
		hmc_float resid;
		get_buffer_from_device(clmem_rho_next, &resid, sizeof(hmc_float));
		//this is the orig. call
		//set_float_to_global_squarenorm_device(clmem_rn, clmem_resid);
		//get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));

		logger.debug() << "resid: " << resid;

		if(resid < prec) {
			// report on performance
			if(logger.beInfo()) {
				// we are always synchroneous here, as we had to recieve the residium from the device
				uint64_t duration = timer.getTime();

				// calculate flops
				unsigned refreshs = iter / get_parameters()->get_iter_refresh() + 1;
				cl_ulong mf_flops = f.get_Flops();

				cl_ulong total_flops = mf_flops + 3 * get_flop_size("scalar_product_eoprec") + 2 * get_flop_size("ratio") + 2 * get_flop_size("product") + 3 * get_flop_size("saxpy_eoprec");
				total_flops *= iter;

				total_flops += refreshs * (mf_flops + get_flop_size("saxpy_eoprec") + get_flop_size("scalar_product_eoprec"));

				// report performanc
				logger.info() << "CG completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f) << " Gflops.";
			}

			return iter;
		}

		//beta = (rn+1, rn+1)/(rn, rn) --> alpha = rho_next/omega
		set_complex_to_scalar_product_eoprec_device(clmem_rn_eoprec, clmem_rn_eoprec, clmem_rho_next);
		set_complex_to_ratio_device(clmem_rho_next, clmem_omega, clmem_beta);

		//pn+1 = rn+1 + beta*pn
		set_complex_to_product_device(clmem_beta, clmem_minusone, clmem_tmp2);
		saxpy_eoprec_device(clmem_p_eoprec, clmem_rn_eoprec, clmem_tmp2, clmem_p_eoprec);
	}
	return -1;
}


void Opencl_Module_Fermions::solver(const Matrix_Function & f, cl_mem inout, cl_mem source, cl_mem gf, usetimer * solvertimer)
{
	/** This solves the sparse-matrix system
	 *  A x = b
	 *  with  x == inout
	 *        A == f
	 *        b == source
	 * using a Krylov-Solver (BiCGStab or CG)
	 */
	int converged = -1;

	if(get_parameters()->get_profile_solver() ) (*solvertimer).reset();
	if(get_parameters()->get_use_eo() == true) {
		/**
		 * If even-odd-preconditioning is used, the inversion is split up
		 * into even and odd parts using Schur decomposition, assigning the
		 * non-trivial inversion to the even sites (see DeGran/DeTar p 174ff).
		 */
		//convert source and input-vector to eoprec-format
		convert_to_eoprec_device(clmem_source_even, clmem_source_odd, source);
		//prepare sources
		/**
		 * This changes the even source according to (with A = M + D):
		 *  b_e = b_e - D_eo M_inv b_o
		 */
		if(get_parameters()->get_fermact() == WILSON) {
			//in this case, the diagonal matrix is just 1 and falls away.
			M_tm_inverse_sitediagonal_device(clmem_source_odd, clmem_tmp_eoprec_1);
			dslash_eoprec_device(clmem_source_odd, clmem_tmp_eoprec_1, gf, EVEN);
			saxpy_eoprec_device(clmem_source_even, clmem_tmp_eoprec_1, clmem_one, clmem_source_even);
		} else if(get_parameters()->get_fermact() == TWISTEDMASS) {
			M_tm_inverse_sitediagonal_device(clmem_source_odd, clmem_tmp_eoprec_1);
			dslash_eoprec_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2, gf, EVEN);
			saxpy_eoprec_device(clmem_source_even, clmem_tmp_eoprec_2, clmem_one, clmem_source_even);
		}

		//Trial solution
		///@todo this should go into a more general function
		this->set_eoprec_spinorfield_cold_device(this->get_clmem_inout_eoprec());
		logger.debug() << "start eoprec-inversion";
		//even solution
		if(get_parameters()->get_use_cg() == true)
			converged = cg_eoprec(f, clmem_inout_eoprec, clmem_source_even, gf, get_parameters()->get_solver_prec());
		else
			converged = bicgstab_eoprec(f, this->get_clmem_inout_eoprec(), clmem_source_even, gf, get_parameters()->get_solver_prec());

		//odd solution
		/** The odd solution is obtained from the even one according to:
		 *  x_o = M_inv D x_e - M_inv b_o
		 */
		if(get_parameters()->get_fermact() == WILSON) {
			//in this case, the diagonal matrix is just 1 and falls away.
			dslash_eoprec_device(clmem_inout_eoprec, clmem_tmp_eoprec_1, gf, ODD);
			saxpy_eoprec_device(clmem_tmp_eoprec_1, clmem_source_odd, clmem_one, clmem_tmp_eoprec_1);
		} else if(get_parameters()->get_fermact() == TWISTEDMASS) {
			dslash_eoprec_device(clmem_inout_eoprec, clmem_tmp_eoprec_2, gf, ODD);
			M_tm_inverse_sitediagonal_device(clmem_tmp_eoprec_2, clmem_tmp_eoprec_1);
			M_tm_inverse_sitediagonal_device(clmem_source_odd, clmem_tmp_eoprec_2);
			saxpy_eoprec_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2, clmem_one, clmem_tmp_eoprec_1);
		}
		//CP: whole solution
		//CP: suppose the even sol is saved in inout_eoprec, the odd one in clmem_tmp_eoprec_1
		convert_from_eoprec_device(clmem_inout_eoprec, clmem_tmp_eoprec_1, inout);
	} else {
		//Trial solution
		///@todo this should go into a more general function
		this->set_spinorfield_cold_device(inout);

		if(get_parameters()->get_use_cg() == true)
			converged = cg(f, inout, source, gf, get_parameters()->get_solver_prec());
		else
			converged = bicgstab(f, inout, source, gf, get_parameters()->get_solver_prec());
	}
	if(get_parameters()->get_profile_solver() ) {
		clFinish(get_queue());
		(*solvertimer).add();
	}

	if (converged < 0) {
		if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
		else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
	} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";

	return;
}

cl_mem Opencl_Module_Fermions::get_clmem_inout()
{
	return clmem_inout;
}

cl_mem Opencl_Module_Fermions::get_clmem_source()
{
	return clmem_source;
}

cl_mem Opencl_Module_Fermions::get_clmem_tmp()
{
	return clmem_tmp;
}

cl_mem Opencl_Module_Fermions::get_clmem_inout_eoprec()
{
	return clmem_inout_eoprec;
}

cl_mem Opencl_Module_Fermions::get_clmem_tmp_eoprec_1()
{
	return clmem_tmp_eoprec_1;
}

cl_mem Opencl_Module_Fermions::get_clmem_tmp_eoprec_2()
{
	return clmem_tmp_eoprec_2;
}

cl_mem Opencl_Module_Fermions::get_clmem_source_even()
{
	return clmem_source_even;
}

cl_mem Opencl_Module_Fermions::get_clmem_source_odd()
{
	return clmem_source_odd;
}

cl_mem Opencl_Module_Fermions::get_clmem_minusone()
{
	return clmem_minusone;
}

cl_mem Opencl_Module_Fermions::get_clmem_one()
{
	return clmem_one;
}

hmc_float Opencl_Module_Fermions::print_info_inv_field(cl_mem in, bool eo, std::string msg)
{
	cl_mem clmem_sqnorm_tmp = create_rw_buffer(sizeof(hmc_float));
	hmc_float tmp;
	if(eo) set_float_to_global_squarenorm_eoprec_device(in, clmem_sqnorm_tmp);
	else set_float_to_global_squarenorm_device(in, clmem_sqnorm_tmp);
	get_buffer_from_device(clmem_sqnorm_tmp, &tmp, sizeof(hmc_float));
	cout.precision(10);
	logger.debug() << std::scientific << msg << tmp;
	int clerr = clReleaseMemObject(clmem_sqnorm_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	return tmp;
}

#ifdef _PROFILING_
usetimer* Opencl_Module_Fermions::get_timer(const char * in)
{
	usetimer *noop = NULL;
	noop = Opencl_Module_Spinors::get_timer(in);
	if(noop != NULL) return noop;

	if (strcmp(in, "M_wilson") == 0) {
		return &(this->timer_M_wilson);
	}
	if (strcmp(in, "gamma5") == 0) {
		return &this->timer_gamma5;
	}
	if (strcmp(in, "M_tm_plus") == 0) {
		return &this->timer_M_tm_plus;
	}
	if (strcmp(in, "M_tm_minus") == 0) {
		return &this->timer_M_tm_minus;
	}
	if (strcmp(in, "gamma5_eoprec") == 0) {
		return &this->timer_gamma5_eoprec;
	}
	if (strcmp(in, "dslash_eoprec") == 0) {
		return &this->timer_dslash_eoprec;
	}
	if (strcmp(in, "M_tm_sitediagonal") == 0) {
		return &this->timer_M_tm_sitediagonal;
	}
	if (strcmp(in, "M_tm_inverse_sitediagonal") == 0) {
		return &this->timer_M_tm_inverse_sitediagonal;
	}
	if (strcmp(in, "M_tm_sitediagonal_minus") == 0) {
		return &this->timer_M_tm_sitediagonal_minus;
	}
	if (strcmp(in, "M_tm_inverse_sitediagonal_minus") == 0) {
		return &this->timer_M_tm_inverse_sitediagonal_minus;
	}
	if (strcmp(in, "ps_correlator") == 0) {
		return &this->timer_ps_correlator;
	}

	//if the kernelname has not matched, return NULL
	else {
		return NULL;
	}
}
#endif

int Opencl_Module_Fermions::get_read_write_size(const char * in)
{
	int result = Opencl_Module_Spinors::get_read_write_size(in);
	if (result != 0) return result;
	//Depending on the compile-options, one has different sizes...
	int D = (*parameters).get_float_size();
	//this returns the number of entries in an su3-matrix
	int R = (*parameters).get_mat_size();
	int S = get_parameters()->get_spinorfieldsize();
	int Seo = get_parameters()->get_eoprec_spinorfieldsize();
	//factor for complex numbers
	int C = 2;
	//this is the same as in the function above
	//NOTE: 1 spinor has NC*NDIM = 12 complex entries
	if (strcmp(in, "M_wilson") == 0) {
		//this kernel reads 9 spinors, 8 su3matrices and writes 1 spinor:
		return (C * 12 * (9 + 1) + C * 8 * R) * D * S;
	}
	if (strcmp(in, "gamma5") == 0) {
		//this kernel reads 1 spinor and writes 1 spinor:
		return 2 * C * 12 * D * S;
	}
	if (strcmp(in, "M_tm_plus") == 0) {
		//this kernel reads 9 spinors, 8 su3matrices and writes 1 spinor:
		return (C * 12 * (9 + 1) + C * 8 * R) * D * S;
	}
	if (strcmp(in, "M_tm_minus") == 0) {
		//this kernel reads 9 spinors, 8 su3matrices and writes 1 spinor:
		return (C * 12 * (9 + 1) + C * 8 * R) * D * S;
	}
	if (strcmp(in, "gamma5_eoprec") == 0) {
		//this kernel reads 1 spinor and writes 1 spinor:
		return 48 * D * Seo;
	}
	if (strcmp(in, "M_tm_sitediagonal") == 0) {
		//this kernel reads 1 spinor and writes 1 spinor:
		return 48 * D * Seo;
	}
	if (strcmp(in, "M_tm_inverse_sitediagonal") == 0) {
		//this kernel reads 1 spinor and writes 1 spinor:
		return 48 * D * Seo;
	}
	if (strcmp(in, "M_tm_sitediagonal_minus") == 0) {
		//this kernel reads 1 spinor and writes 1 spinor:
		return 48 * D * Seo;
	}
	if (strcmp(in, "M_tm_inverse_sitediagonal_minus") == 0) {
		//this kernel reads 1 spinor and writes 1 spinor:
		return 48 * D * Seo;
	}
	if (strcmp(in, "dslash_eoprec") == 0) {
		//this kernel reads 8 spinors, 8 su3matrices and writes 1 spinor:
		return (C * 12 * (8 + 1) + C * 8 * R) * D * Seo;
	}
	return 0;
}

int flop_dslash_per_site(inputparameters * parameters)
{
	/** @NOTE: this is the "original" dslash without any simplifications, counting everything "full". this is a much too hight number!!
	   *  //this kernel performs for each eo site a 2*NDIM sum over (1 + gamma_mu) * su3matrix * spinor
	   *  //return  NDIM * 2 * ( get_parameters()->get_flop_su3_spinor() + get_parameters()->get_flop_gamma_spinor() ) ;
	  @NOTE: However, in 0911.3191 the dslash_eo is quoted to perform 1320 flop per site
	   *  If I count our implementation of the dslash-kernel, I get 1632 flop per site:
	   *  //the kernel performs 6 su3vec_acc, 2 su3_su3vec and 2 su3vec_complex in NDIM * 2 directions per site
	  */
	return NDIM * 2 * (parameters->get_flop_su3_su3vec() * 2 + 6 * NC * 2 + 2 * NC * parameters->get_flop_complex_mult() );
	/// I do not know for sure, but if one leaves out the 2 su3vec_complex operations in each directions one has almost 1320 flop per site. This would correspond to have the kappa in the diagonal matrix still.

}

int Opencl_Module_Fermions::get_flop_size(const char * in)
{
	int result = Opencl_Module_Spinors::get_flop_size(in);
	if (result != 0) return result;
	int S = get_parameters()->get_spinorfieldsize();
	int Seo = get_parameters()->get_eoprec_spinorfieldsize();
	if (strcmp(in, "M_wilson") == 0) {
		//this kernel performs one dslash on each site and adds this to a spinor
		return S * (flop_dslash_per_site(get_parameters()) + NC * NDIM * get_parameters()->get_flop_complex_mult() + NC * NDIM * 2 );
	}
	if (strcmp(in, "gamma5") == 0) {
		//this kernel performs ND*NC*2/2 real mults
		return S * NDIM * NC;
	}
	if (strcmp(in, "M_tm_plus") == 0) {
		//this kernel performs ND*NC complex mults and one dslash on each site and adds the results
		return S * (flop_dslash_per_site(get_parameters()) + NC * NDIM * get_parameters()->get_flop_complex_mult() + NC * NDIM * 2 );
	}
	if (strcmp(in, "M_tm_minus") == 0) {
		//this kernel performs ND*NC complex mults and one dslash on each site and adds the results
		return S * (flop_dslash_per_site(get_parameters()) + NC * NDIM * get_parameters()->get_flop_complex_mult() + NC * NDIM * 2 );
	}
	if (strcmp(in, "gamma5_eoprec") == 0) {
		//this kernel performs ND*NC*2/2 real mults
		return Seo * NDIM * NC;
	}
	if (strcmp(in, "M_tm_sitediagonal") == 0) {
		//this kernel performs ND*NC complex mults
		return Seo * ( NC * NDIM * get_parameters()->get_flop_complex_mult() );
	}
	if (strcmp(in, "M_tm_inverse_sitediagonal") == 0) {
		//this kernel performs ND*NC complex mults and ND*NC*2 real mults
		return Seo * ( NC * NDIM * get_parameters()->get_flop_complex_mult() + NC * NDIM * 2  );
	}
	if (strcmp(in, "M_tm_sitediagonal_minus") == 0) {
		//this kernel performs ND*NC complex mults
		return Seo * ( NC * NDIM * get_parameters()->get_flop_complex_mult() );
	}
	if (strcmp(in, "M_tm_inverse_sitediagonal_minus") == 0) {
		//this kernel performs ND*NC complex mults and ND*NC*2 real mults
		return Seo * ( NC * NDIM * get_parameters()->get_flop_complex_mult() + NC * NDIM * 2 );
	}
	if (strcmp(in, "dslash_eoprec") == 0) {
		return Seo * flop_dslash_per_site(get_parameters());
	}
	return 0;
}

#ifdef _PROFILING_
void Opencl_Module_Fermions::print_profiling(std::string filename, int number)
{
	Opencl_Module_Spinors::print_profiling(filename, number);
	const char * kernelName;
	kernelName = "M_wilson";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "gamma5";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "M_tm_plus";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "M_tm_minus";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "gamma5_eoprec";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "M_tm_sitediagonal";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "M_tm_inverse_sitediagonal";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "M_tm_sitediagonal_minus";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "M_tm_inverse_sitediagonal_minus";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "dslash_eoprec";
	Opencl_Module::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
}
#endif
