#include "opencl_module_fermions.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"
#include "meta/util.hpp"

using namespace std;

/**
 * What follows are functions that call opencl_fermions-class-functions.
 * This is needed to be able to pass different fermionmatrices as
 *  arguments to class-functions.
 */
void M::operator()(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa, hmc_float mubar) const
{
	that->M(in, out, gf, kappa, mubar);
}
cl_ulong M::get_Flops() const
{
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			return that->get_flop_size("M_wilson");
		case meta::Inputparameters::twistedmass:
			return that->get_flop_size("M_tm_plus");
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
}
cl_ulong M::get_Bytes() const
{
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			return that->get_read_write_size("M_wilson");
		case meta::Inputparameters::twistedmass:
			return that->get_read_write_size("M_tm_plus");
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
}

void Qplus::operator()(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa, hmc_float mubar) const
{
	that->Qplus(in, out, gf, kappa, mubar);
}
cl_ulong Qplus::get_Flops() const
{
	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = that->get_flop_size("M_wilson");
			break;
		case meta::Inputparameters::twistedmass:
			res = that->get_flop_size("M_tm_plus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += that->get_flop_size("gamma5");
	return res;
}
cl_ulong Qplus::get_Bytes() const
{
	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = that->get_read_write_size("M_wilson");
			break;
		case meta::Inputparameters::twistedmass:
			res = that->get_read_write_size("M_tm_plus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += that->get_read_write_size("gamma5");
	return res;
}

void Qminus::operator()(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa, hmc_float mubar) const
{
	that->Qminus(in, out, gf, kappa, mubar);
}
cl_ulong Qminus::get_Flops() const
{
	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = that->get_flop_size("M_wilson");
			break;
		case meta::Inputparameters::twistedmass:
			res = that->get_flop_size("M_tm_minus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += that->get_flop_size("gamma5");
	return res;
}
cl_ulong Qminus::get_Bytes() const
{
	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = that->get_read_write_size("M_wilson");
			break;
		case meta::Inputparameters::twistedmass:
			res = that->get_read_write_size("M_tm_minus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += that->get_read_write_size("gamma5");
	return res;
}

void QplusQminus::operator()(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa, hmc_float mubar) const
{
	that->QplusQminus(in, out, gf, kappa, mubar);
}
cl_ulong QplusQminus::get_Flops() const
{
	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_flop_size("M_wilson");
			break;
		case meta::Inputparameters::twistedmass:
			res = that->get_flop_size("M_tm_plus");
			res += that->get_flop_size("M_tm_minus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += 2 * that->get_flop_size("gamma5");
	return res;
}
cl_ulong QplusQminus::get_Bytes() const
{
	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_read_write_size("M_wilson");
			break;
		case meta::Inputparameters::twistedmass:
			res = that->get_read_write_size("M_tm_plus");
			res += that->get_read_write_size("M_tm_minus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += 2 * that->get_read_write_size("gamma5");
	return res;
}

void Aee::operator()(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa, hmc_float mubar) const
{
	that->Aee(in, out, gf, kappa, mubar);
}
cl_ulong Aee::get_Flops() const
{
	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_flop_size("dslash_eo");
			res += that->get_flop_size("saxpy_eoprec");
			break;
		case meta::Inputparameters::twistedmass:
			res = 2 * that->get_flop_size("dslash_eo");
			res += that->get_flop_size("M_tm_inverse_sitediagonal");
			res += that->get_flop_size("M_tm_sitediagonal");
			res += that->get_flop_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	return res;
}
cl_ulong Aee::get_Bytes() const
{
	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_read_write_size("dslash_eo");
			res += that->get_read_write_size("saxpy_eoprec");
			break;
		case meta::Inputparameters::twistedmass:
			res = 2 * that->get_read_write_size("dslash_eo");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal");
			res += that->get_read_write_size("M_tm_sitediagonal");
			res += that->get_read_write_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	return res;
}

void Qplus_eo::operator()(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa, hmc_float mubar) const
{
	that->Qplus_eo(in, out, gf, kappa, mubar);
}
cl_ulong Qplus_eo::get_Flops() const
{
	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_flop_size("dslash_eo");
			res += that->get_flop_size("saxpy_eoprec");
			break;
		case meta::Inputparameters::twistedmass:
			res = 2 * that->get_flop_size("dslash_eo");
			res += that->get_flop_size("M_tm_inverse_sitediagonal");
			res += that->get_flop_size("M_tm_sitediagonal");
			res += that->get_flop_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += that->get_flop_size("gamma5_eo");
	return res;
}
cl_ulong Qplus_eo::get_Bytes() const
{
	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_read_write_size("dslash_eo");
			res += that->get_read_write_size("saxpy_eoprec");
			break;
		case meta::Inputparameters::twistedmass:
			res = 2 * that->get_read_write_size("dslash_eo");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal");
			res += that->get_read_write_size("M_tm_sitediagonal");
			res += that->get_read_write_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += that->get_flop_size("gamma5_eo");
	return res;
}

void Qminus_eo::operator()(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa, hmc_float mubar) const
{
	that->Qminus_eo(in, out, gf, kappa, mubar);
}
cl_ulong Qminus_eo::get_Flops() const
{
	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_flop_size("dslash_eo");
			res += that->get_flop_size("saxpy_eoprec");
			break;
		case meta::Inputparameters::twistedmass:
			res = 2 * that->get_flop_size("dslash_eo");
			res += that->get_flop_size("M_tm_inverse_sitediagonal_minus");
			res += that->get_flop_size("M_tm_sitediagonal_minus");
			res += that->get_flop_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += that->get_flop_size("gamma5_eo");
	return res;
}
cl_ulong Qminus_eo::get_Bytes() const
{
	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_read_write_size("dslash_eo");
			res += that->get_read_write_size("saxpy_eoprec");
			break;
		case meta::Inputparameters::twistedmass:
			res = 2 * that->get_read_write_size("dslash_eo");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal_minus");
			res += that->get_read_write_size("M_tm_sitediagonal_minus");
			res += that->get_read_write_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += that->get_flop_size("gamma5_eo");
	return res;
}

void QplusQminus_eo::operator()(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa, hmc_float mubar) const
{
	that->QplusQminus_eo(in, out, gf, kappa, mubar);
}
cl_ulong QplusQminus_eo::get_Flops() const
{
	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_flop_size("dslash_eo");
			res += that->get_flop_size("saxpy_eoprec");
			res *= 2;
			break;
		case meta::Inputparameters::twistedmass:
			res = 4 * that->get_flop_size("dslash_eo");
			res += that->get_flop_size("M_tm_inverse_sitediagonal");
			res += that->get_flop_size("M_tm_sitediagonal");
			res += that->get_flop_size("M_tm_inverse_sitediagonal_minus");
			res += that->get_flop_size("M_tm_sitediagonal_minus");
			res += 2 * that->get_flop_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += 2 * that->get_flop_size("gamma5_eo");
	return res;
}
cl_ulong QplusQminus_eo::get_Bytes() const
{
	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_read_write_size("dslash_eo");
			res += that->get_read_write_size("saxpy_eoprec");
			res *= 2;
			break;
		case meta::Inputparameters::twistedmass:
			res = 4 * that->get_read_write_size("dslash_eo");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal");
			res += that->get_read_write_size("M_tm_sitediagonal");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal_minus");
			res += that->get_read_write_size("M_tm_sitediagonal_minus");
			res += 2 * that->get_read_write_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += 2 * that->get_flop_size("gamma5_eo");
	return res;
}


void Opencl_Module_Fermions::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module_Spinors::fill_collect_options(collect_options);

	switch (get_parameters().get_fermact()) {
		case meta::Inputparameters::twistedmass :
			*collect_options << " -D_TWISTEDMASS_";
			break;
		case meta::Inputparameters::clover :
			*collect_options << " -D_CLOVER_";
			break;
	}

	//CP: These are the BCs in spatial and temporal direction
	hmc_float tmp_spatial = (get_parameters().get_theta_fermion_spatial() * PI) / ( (hmc_float) get_parameters().get_nspace());
	hmc_float tmp_temporal = (get_parameters().get_theta_fermion_temporal() * PI) / ( (hmc_float) get_parameters().get_ntime());
	//BC: on the corners in each direction: exp(i theta) -> on each site exp(i theta*PI /LATEXTENSION) = cos(tmp2) + isin(tmp2)
	*collect_options << " -DSPATIAL_RE=" << cos(tmp_spatial);
	*collect_options << " -DMSPATIAL_RE=" << -cos(tmp_spatial);
	*collect_options << " -DSPATIAL_IM=" << sin(tmp_spatial);
	*collect_options << " -DMSPATIAL_IM=" << -sin(tmp_spatial);

	*collect_options << " -DTEMPORAL_RE=" << cos(tmp_temporal);
	*collect_options << " -DMTEMPORAL_RE=" << -cos(tmp_temporal);
	*collect_options << " -DTEMPORAL_IM=" << sin(tmp_temporal);
	*collect_options << " -DMTEMPORAL_IM=" << -sin(tmp_temporal);

	return;
}


void Opencl_Module_Fermions::fill_buffers()
{
	Opencl_Module_Spinors::fill_buffers();

	int clerr = CL_SUCCESS;

	int spinorfield_size = sizeof(spinor) * meta::get_spinorfieldsize(get_parameters());
	int eoprec_spinorfield_buffer_size = get_eoprec_spinorfield_buffer_size();
	int complex_size = sizeof(hmc_complex);
	int float_size = sizeof(hmc_float);

	logger.debug() << "init general spinorfield-buffers";
	clmem_inout = create_rw_buffer(spinorfield_size);
	clmem_tmp = create_rw_buffer(spinorfield_size);
	clmem_source = create_rw_buffer(spinorfield_size);

	if(get_parameters().get_use_eo() == true) {
		logger.debug() << "init general eoprec-spinorfield-buffers";
		clmem_inout_eo = create_rw_buffer(eoprec_spinorfield_buffer_size);
		clmem_source_even = create_rw_buffer(eoprec_spinorfield_buffer_size);
		clmem_source_odd = create_rw_buffer(eoprec_spinorfield_buffer_size);
		clmem_tmp_eo_1 = create_rw_buffer(eoprec_spinorfield_buffer_size);
		//this field is used only with twistedmass
		if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			clmem_tmp_eo_2 = create_rw_buffer(eoprec_spinorfield_buffer_size);
		} else
			clmem_tmp_eo_2 = 0;
	}

	clmem_one = create_rw_buffer(complex_size);
	clmem_minusone = create_rw_buffer(complex_size);
	//fill buffers with one and minusone
	hmc_complex one = hmc_complex_one;
	hmc_complex minusone = hmc_complex_minusone;
	clerr = clEnqueueWriteBuffer(get_queue(), clmem_one, CL_TRUE, 0, complex_size, &one, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);
	clerr = clEnqueueWriteBuffer(get_queue(), clmem_minusone, CL_TRUE, 0, complex_size, &minusone, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);

	//fill buffers for solver
	this->fill_solver_buffers();

}

void Opencl_Module_Fermions::fill_solver_buffers()
{
	int clerr = CL_SUCCESS;

	int spinorfield_size = sizeof(spinor) * get_spinorfieldsize(get_parameters());
	int eoprec_spinorfield_buffer_size = get_eoprec_spinorfield_buffer_size();
	int complex_size = sizeof(hmc_complex);
	int float_size = sizeof(hmc_float);

	///@todo some buffers can be saved here if only cg is used
	if(get_parameters().get_use_eo() == false) {
		logger.debug() << "init solver spinorfield-buffers";
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
		clmem_rn_eo = create_rw_buffer(eoprec_spinorfield_buffer_size);
		clmem_rhat_eo = create_rw_buffer(eoprec_spinorfield_buffer_size);
		clmem_v_eo = create_rw_buffer(eoprec_spinorfield_buffer_size);
		clmem_p_eo = create_rw_buffer(eoprec_spinorfield_buffer_size);
		clmem_s_eo = create_rw_buffer(eoprec_spinorfield_buffer_size);
		clmem_t_eo = create_rw_buffer(eoprec_spinorfield_buffer_size);
		clmem_aux_eo = create_rw_buffer(eoprec_spinorfield_buffer_size);
	}

	logger.debug() << "create buffers for complex and real numbers";
	clmem_rho = create_rw_buffer(complex_size);
	clmem_rho_next = create_rw_buffer(complex_size);
	clmem_alpha = create_rw_buffer(complex_size);
	clmem_omega = create_rw_buffer(complex_size);
	clmem_beta = create_rw_buffer(complex_size);
	clmem_tmp1 = create_rw_buffer(complex_size);
	clmem_tmp2 = create_rw_buffer(complex_size);
	clmem_resid = create_rw_buffer(float_size);
	clmem_trueresid = create_rw_buffer(float_size);

	logger.debug() << "write contents to some buffers";

	hmc_complex zero = hmc_complex_zero;
	clerr = clEnqueueWriteBuffer(get_queue(), clmem_resid, CL_TRUE, 0, float_size, &zero, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);
	clerr = clEnqueueWriteBuffer(get_queue(), clmem_trueresid, CL_TRUE, 0, float_size, &zero, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);

	logger.trace() << "init solver buffers done...";
}


void Opencl_Module_Fermions::fill_kernels()
{
	Opencl_Module_Spinors::fill_kernels();

	M_wilson = 0;
	M_tm_plus = 0;
	M_tm_minus = 0;

	logger.debug() << "Create fermion kernels...";
	if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
		M_wilson = createKernel("M_wilson") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_m.cl";
	} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
		M_tm_plus = createKernel("M_tm_plus") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_m_tm_plus.cl";
		M_tm_minus = createKernel("M_tm_minus") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_m_tm_minus.cl";
	} else if(get_parameters().get_fermact() == meta::Inputparameters::clover) {
		throw Print_Error_Message("no kernels for CLOVER-discretization implemented yet, aborting... ", __FILE__, __LINE__);
	} else {
		throw Print_Error_Message("there was a problem with which fermion-discretization to use, aborting... ", __FILE__, __LINE__);
	}

	gamma5 = createKernel("gamma5") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_gamma5.cl";


	//Kernels needed if eoprec is used
	if(get_parameters().get_use_eo() == true) {
		if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			M_tm_sitediagonal = createKernel("M_tm_sitediagonal") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m.cl";
			M_tm_inverse_sitediagonal = createKernel("M_tm_inverse_sitediagonal") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m.cl";
			M_tm_sitediagonal_minus = createKernel("M_tm_sitediagonal_minus") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m.cl";
			M_tm_inverse_sitediagonal_minus = createKernel("M_tm_inverse_sitediagonal_minus") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m.cl";
		}
		dslash_eo = createKernel("dslash_eo") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_dslash.cl";
		gamma5_eo = createKernel("gamma5_eo") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_eo_gamma5.cl";
		//merged kernels
		if (get_parameters().get_use_merge_kernels_fermion() == true) {
			dslash_AND_gamma5_eo = createKernel("dslash_AND_gamma5_eo") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_dslash_AND_gamma5.cl";
			dslash_AND_M_tm_inverse_sitediagonal_eo = createKernel("dslash_AND_M_tm_inverse_sitediagonal_eo") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_dslash_AND_M_tm_inverse_sitediagonal.cl";
			dslash_AND_M_tm_inverse_sitediagonal_minus_eo = createKernel("dslash_AND_M_tm_inverse_sitediagonal_minus_eo") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_dslash_AND_M_tm_inverse_sitediagonal_minus.cl";
			M_tm_sitediagonal_AND_gamma5_eo = createKernel("M_tm_sitediagonal_AND_gamma5_eo") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m_merged.cl";
			M_tm_sitediagonal_minus_AND_gamma5_eo = createKernel("M_tm_sitediagonal_minus_AND_gamma5_eo") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m_merged.cl";
		}

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
	if(gamma5) {
		clerr = clReleaseKernel(gamma5);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}

	if(get_parameters().get_use_eo()) {
		clerr = clReleaseKernel(dslash_eo);
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

	return;
}

void Opencl_Module_Fermions::clear_buffers()
{
	Opencl_Module_Spinors::clear_buffers();

	cl_uint clerr = CL_SUCCESS;

	logger.trace() << "clear general fermion buffer...";

	if(clmem_inout) {
		clerr = clReleaseMemObject(clmem_inout);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clmem_inout = 0;
	}
	if(clmem_source) {
		clerr = clReleaseMemObject(clmem_source);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clmem_source = 0;
	}
	if(clmem_tmp) {
		clerr = clReleaseMemObject(clmem_tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clmem_tmp = 0;
	}

	if(get_parameters().get_use_eo()) {
		if(clmem_inout_eo) {
			clerr = clReleaseMemObject(clmem_inout_eo);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_inout_eo = 0;
		}
		if(clmem_source_even) {
			clerr = clReleaseMemObject(clmem_source_even);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_source_even = 0;
		}
		if(clmem_source_odd) {
			clerr = clReleaseMemObject(clmem_source_odd);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_source_odd = 0;
		}
		if(clmem_tmp_eo_1) {
			clerr = clReleaseMemObject(clmem_tmp_eo_1);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_tmp_eo_1 = 0;
		}
		if(clmem_tmp_eo_2) {
			clerr = clReleaseMemObject(clmem_tmp_eo_2);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_tmp_eo_2 = 0;
		}
	}

	if(clmem_minusone) {
		clerr = clReleaseMemObject(clmem_minusone);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clmem_minusone = 0;
	}
	if(clmem_one) {
		clerr = clReleaseMemObject(clmem_one);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clmem_one = 0;
	}

	//clear solver buffers
	this->clear_solver_buffers();

	return;
}

void Opencl_Module_Fermions::clear_solver_buffers()
{
	cl_uint clerr = CL_SUCCESS;

	logger.trace() << "clear solver buffer...";

	if(get_parameters().get_use_eo()) {
		if(clmem_rn_eo) {
			clerr = clReleaseMemObject(clmem_rn_eo);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_rn_eo = 0;
		}
		if(clmem_rhat_eo) {
			clerr = clReleaseMemObject(clmem_rhat_eo);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_rhat_eo = 0;
		}
		if(clmem_v_eo) {
			clerr = clReleaseMemObject(clmem_v_eo);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_v_eo = 0;
		}
		if(clmem_p_eo) {
			clerr = clReleaseMemObject(clmem_p_eo);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_p_eo = 0;
		}
		if(clmem_s_eo) {
			clerr = clReleaseMemObject(clmem_s_eo);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_s_eo = 0;
		}
		if(clmem_t_eo) {
			clerr = clReleaseMemObject(clmem_t_eo);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_t_eo = 0;
		}
		if(clmem_aux_eo) {
			clerr = clReleaseMemObject(clmem_aux_eo);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_aux_eo = 0;
		}
	} else {
		if(clmem_rn) {
			clerr = clReleaseMemObject(clmem_rn);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_rn = 0;
		}
		if(clmem_rhat) {
			clerr = clReleaseMemObject(clmem_rhat);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_rhat = 0;
		}
		if(clmem_v) {
			clerr = clReleaseMemObject(clmem_v);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_v = 0;
		}
		if(clmem_p) {
			clerr = clReleaseMemObject(clmem_p);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_p = 0;
		}
		if(clmem_s) {
			clerr = clReleaseMemObject(clmem_s);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_s = 0;
		}
		if(clmem_t) {
			clerr = clReleaseMemObject(clmem_t);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_t = 0;
		}
		if(clmem_aux) {
			clerr = clReleaseMemObject(clmem_aux);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
			clmem_aux = 0;
		}
	}

	if(clmem_rho) {
		clerr = clReleaseMemObject(clmem_rho);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clmem_rho = 0;
	}
	if(clmem_rho_next) {
		clerr = clReleaseMemObject(clmem_rho_next);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clmem_rho_next = 0;
	}
	if(clmem_omega) {
		clerr = clReleaseMemObject(clmem_omega);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clmem_omega = 0;
	}
	if(clmem_beta) {
		clerr = clReleaseMemObject(clmem_beta);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clmem_beta = 0;
	}
	if(clmem_tmp2) {
		clerr = clReleaseMemObject(clmem_tmp2);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clmem_tmp2 = 0;
	}
	if(clmem_resid) {
		clerr = clReleaseMemObject(clmem_resid);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clmem_resid = 0;
	}
	if(clmem_trueresid) {
		clerr = clReleaseMemObject(clmem_trueresid);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clmem_trueresid = 0;
	}
	if(clmem_alpha) {
		clerr = clReleaseMemObject(clmem_alpha);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clmem_alpha = 0;
	}
	if(clmem_tmp1) {
		clerr = clReleaseMemObject(clmem_tmp1);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
		clmem_tmp1 = 0;
	}

	return;
}

void Opencl_Module_Fermions::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module_Spinors::get_work_sizes(kernel, ls, gs, num_groups);

	return;
}


//compound fermionmatrix-functions without eoprec
void Opencl_Module_Fermions::M(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa , hmc_float mubar )
{

	if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
		//in the pure Wilson case there is just one fermionmatrix
		M_wilson_device(in, out, gf, kappa);
	} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
		M_tm_plus_device(in, out, gf, kappa, mubar);
	}
}

void Opencl_Module_Fermions::Qplus(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa , hmc_float mubar )
{
	if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
		//in the pure Wilson case there is just one fermionmatrix
		M_wilson_device(in, out, gf, kappa);
	} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
		M_tm_plus_device(in, out, gf, kappa, mubar);
	}
	gamma5_device(out);
}

void Opencl_Module_Fermions::Qminus(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa , hmc_float mubar )
{
	if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
		//in the pure Wilson case there is just one fermionmatrix
		M_wilson_device(in, out, gf, kappa);
	} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
		M_tm_minus_device(in, out, gf, kappa, mubar);
	}
	gamma5_device(out);
}

void Opencl_Module_Fermions::QplusQminus(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa , hmc_float mubar )
{
	/** @todo one could save one field here if an additional copying would be included in the end...
	 * or the field should be created in here, local */
	Qminus(in, clmem_tmp, gf, kappa, mubar);
	Qplus(clmem_tmp, out, gf, kappa, mubar);
}

//explicit fermionmatrix-kernel calling functions
void Opencl_Module_Fermions::M_wilson_device(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa)
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
	int clerr = clSetKernelArg(M_wilson, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_wilson, 1, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_wilson, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_wilson, 3, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( M_wilson, gs2, ls2);
}

void Opencl_Module_Fermions::M_tm_plus_device(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa , hmc_float mubar )
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
	int clerr = clSetKernelArg(M_tm_plus, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_plus, 1, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_plus, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_plus, 3, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_plus, 4, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( M_tm_plus, gs2, ls2);
}

void Opencl_Module_Fermions::M_tm_minus_device(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa , hmc_float mubar )
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
	int clerr = clSetKernelArg(M_tm_minus, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_minus, 1, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_minus, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_minus, 3, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_minus, 4, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( M_tm_minus, gs2, ls2);
}

void Opencl_Module_Fermions::gamma5_device(cl_mem inout)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gamma5, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gamma5, 0, sizeof(cl_mem), &inout);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(gamma5 , gs2, ls2);
}

//compound fermionmatrix-functions with eoprec
void Opencl_Module_Fermions::Aee(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa , hmc_float mubar )
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
	if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
		//in this case, the diagonal matrix is just 1 and falls away.
		dslash_eo_device(in, clmem_tmp_eo_1, gf, odd, kappa);
		dslash_eo_device(clmem_tmp_eo_1, out, gf, even, kappa);
		saxpy_eoprec_device(out, in, clmem_one, out);
	} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
		dslash_eo_device(in, clmem_tmp_eo_1, gf, odd, kappa);
		M_tm_inverse_sitediagonal_device(clmem_tmp_eo_1, clmem_tmp_eo_2, mubar);
		dslash_eo_device(clmem_tmp_eo_2, out, gf, even, kappa);
		M_tm_sitediagonal_device(in, clmem_tmp_eo_1, mubar);
		saxpy_eoprec_device(out, clmem_tmp_eo_1, clmem_one, out);
	}
}

/**
 *  This is the equivalent of the above function, but for the lower
 *  flavour, which essentially means mu -> -mu in the tm-case and
 *  no changes in the meta::Inputparameters::wilson case.
 */
void Opencl_Module_Fermions::Aee_minus(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa , hmc_float mubar )
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
	if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
		//in this case, the diagonal matrix is just 1 and falls away.
		dslash_eo_device(in, clmem_tmp_eo_1, gf, odd, kappa);
		dslash_eo_device(clmem_tmp_eo_1, out, gf, even, kappa);
		saxpy_eoprec_device(out, in, clmem_one, out);
	} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
		dslash_eo_device(in, clmem_tmp_eo_1, gf, odd, kappa);
		M_tm_inverse_sitediagonal_minus_device(clmem_tmp_eo_1, clmem_tmp_eo_2, mubar);
		dslash_eo_device(clmem_tmp_eo_2, out, gf, even, kappa);
		M_tm_sitediagonal_minus_device(in, clmem_tmp_eo_1, mubar);
		saxpy_eoprec_device(out, clmem_tmp_eo_1, clmem_one, out);
	}
}

//merged fermionmatrix-functions with eoprec
void Opencl_Module_Fermions::Aee_AND_gamma5_eo(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa , hmc_float mubar )
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
	if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
		//in this case, the diagonal matrix is just 1 and falls away.
		//this case has not been adjusted for the merged kernels yet...
		logger.warn() << "No merged kernels implemented for pure Wilson case!";
		dslash_eo_device(in, clmem_tmp_eo_1, gf, odd, kappa);
		dslash_eo_device(clmem_tmp_eo_1, out, gf, even, kappa);
		saxpy_eoprec_device(out, in, clmem_one, out);
		gamma5_eo_device(out);
	} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
		/*
		dslash_eo_device(in, clmem_tmp_eo_1, gf, odd, kappa);
		M_tm_inverse_sitediagonal_device(clmem_tmp_eo_1, clmem_tmp_eo_2, mubar);
		*/
		dslash_AND_M_tm_inverse_sitediagonal_eo_device(in, clmem_tmp_eo_2, gf, odd, kappa, mubar);
		/*
		dslash_eo_device(clmem_tmp_eo_2, out, gf, even, kappa);
		gamma5_eo_device(out);
		*/
		dslash_AND_gamma5_eo_device(clmem_tmp_eo_2, out, gf, even, kappa);
		/*
		M_tm_sitediagonal_device(in, clmem_tmp_eo_1, mubar);
		gamma5_eo_device(clmem_tmp_eo_1);
		*/
		M_tm_sitediagonal_AND_gamma5_eo_device(in, clmem_tmp_eo_1, mubar);
		saxpy_eoprec_device(out, clmem_tmp_eo_1, clmem_one, out);
	}
}

/**
 *  This is the equivalent of the above function, but for the lower
 *  flavour, which essentially means mu -> -mu in the tm-case and
 *  no changes in the meta::Inputparameters::wilson case.
 */
void Opencl_Module_Fermions::Aee_minus_AND_gamma5_eo(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa , hmc_float mubar )
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
	if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
		//in this case, the diagonal matrix is just 1 and falls away.
		//this case has not been adjusted for the merged kernels yet...
		logger.warn() << "No merged kernels implemented for pure Wilson case!";
		dslash_eo_device(in, clmem_tmp_eo_1, gf, odd, kappa);
		dslash_eo_device(clmem_tmp_eo_1, out, gf, even, kappa);
		saxpy_eoprec_device(out, in, clmem_one, out);
		gamma5_eo_device(out);
	} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
		/*
		dslash_eo_device(in, clmem_tmp_eo_1, gf, odd, kappa);
		M_tm_inverse_sitediagonal_minus_device(clmem_tmp_eo_1, clmem_tmp_eo_2, mubar);
		*/
		dslash_AND_M_tm_inverse_sitediagonal_minus_eo_device(in, clmem_tmp_eo_2, gf, odd, kappa, mubar);
		/*
		dslash_eo_device(clmem_tmp_eo_2, out, gf, even, kappa);
		gamma5_eo_device(out);
		*/
		dslash_AND_gamma5_eo_device(clmem_tmp_eo_2, out, gf, even, kappa);
		/*
		M_tm_sitediagonal_minus_device(in, clmem_tmp_eo_1, mubar);
		gamma5_eo_device(clmem_tmp_eo_1);
		*/
		M_tm_sitediagonal_minus_AND_gamma5_eo_device(in, clmem_tmp_eo_1, mubar);
		saxpy_eoprec_device(out, clmem_tmp_eo_1, clmem_one, out);
	}
}


void Opencl_Module_Fermions::Qplus_eo(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa , hmc_float mubar )
{
	if(get_parameters().get_use_merge_kernels_fermion() == false) {
		Aee(in, out, gf, kappa, mubar);
		gamma5_eo_device(out);
	} else {
		Aee_AND_gamma5_eo(in, out, gf, kappa, mubar);
	}
	return;
}

void Opencl_Module_Fermions::Qminus_eo(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa , hmc_float mubar )
{
	if(get_parameters().get_use_merge_kernels_fermion() == false) {
		Aee_minus(in, out, gf, kappa, mubar);
		gamma5_eo_device(out);
	} else {
		Aee_minus_AND_gamma5_eo(in, out, gf, kappa, mubar);
	}
	return;
}

void Opencl_Module_Fermions::QplusQminus_eo(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa , hmc_float mubar )
{
	//CP: this should be an eoprec-sized field. However, this induces problems in the CG algorithm!!!
	//MB: This is because of padding, the eoprec buffer size shoulw always be queried from Opencl_Module_Spinor
	size_t spinorfield_size = get_eoprec_spinorfield_buffer_size();
	cl_mem sf_eo_tmp;
	sf_eo_tmp = create_rw_buffer(spinorfield_size);

	Qminus_eo(in, sf_eo_tmp, gf, kappa, mubar);
	Qplus_eo(sf_eo_tmp, out, gf, kappa, mubar);

	int clerr = clReleaseMemObject(sf_eo_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);

	return;
}

//explicit eoprec fermionmatrix functions
void Opencl_Module_Fermions::gamma5_eo_device(cl_mem inout)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gamma5_eo, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gamma5_eo, 0, sizeof(cl_mem), &inout);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( gamma5_eo, gs2, ls2);
}

void Opencl_Module_Fermions::dslash_eo_device(cl_mem in, cl_mem out, cl_mem gf, int evenodd, hmc_float kappa)
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
	int clerr = clSetKernelArg(dslash_eo, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_eo, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_eo, 2, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_eo, 3, sizeof(cl_int), &eo);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_eo, 4, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(dslash_eo , gs2, ls2);
}

void Opencl_Module_Fermions::dslash_AND_gamma5_eo_device(cl_mem in, cl_mem out, cl_mem gf, int evenodd, hmc_float kappa)
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
	int clerr = clSetKernelArg(dslash_AND_gamma5_eo, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_gamma5_eo, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_gamma5_eo, 2, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_gamma5_eo, 3, sizeof(cl_int), &eo);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_gamma5_eo, 4, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(dslash_AND_gamma5_eo , gs2, ls2);
}

void Opencl_Module_Fermions::dslash_AND_M_tm_inverse_sitediagonal_eo_device(cl_mem in, cl_mem out, cl_mem gf, int evenodd, hmc_float kappa, hmc_float mubar)
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
	int clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_eo, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_eo, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_eo, 2, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_eo, 3, sizeof(cl_int), &eo);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_eo, 4, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_eo, 5, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(dslash_AND_M_tm_inverse_sitediagonal_eo , gs2, ls2);
}

void Opencl_Module_Fermions::dslash_AND_M_tm_inverse_sitediagonal_minus_eo_device(cl_mem in, cl_mem out, cl_mem gf, int evenodd, hmc_float kappa, hmc_float mubar)
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
	int clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, 2, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, 3, sizeof(cl_int), &eo);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, 4, sizeof(hmc_float), &kappa_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(dslash_AND_M_tm_inverse_sitediagonal_minus_eo, 5, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(dslash_AND_M_tm_inverse_sitediagonal_minus_eo , gs2, ls2);
}

void Opencl_Module_Fermions::M_tm_inverse_sitediagonal_device(cl_mem in, cl_mem out, hmc_float mubar)
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
	int clerr = clSetKernelArg(M_tm_inverse_sitediagonal, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_inverse_sitediagonal, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_inverse_sitediagonal, 2, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( M_tm_inverse_sitediagonal, gs2, ls2);
}

void Opencl_Module_Fermions::M_tm_sitediagonal_device(cl_mem in, cl_mem out, hmc_float mubar)
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
	int clerr = clSetKernelArg(M_tm_sitediagonal, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal, 2, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(M_tm_sitediagonal , gs2, ls2);
}

void Opencl_Module_Fermions::M_tm_sitediagonal_AND_gamma5_eo_device(cl_mem in, cl_mem out, hmc_float mubar)
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
	int clerr = clSetKernelArg(M_tm_sitediagonal_AND_gamma5_eo, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal_AND_gamma5_eo, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal_AND_gamma5_eo, 2, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(M_tm_sitediagonal_AND_gamma5_eo , gs2, ls2);
}

void Opencl_Module_Fermions::M_tm_sitediagonal_minus_AND_gamma5_eo_device(cl_mem in, cl_mem out, hmc_float mubar)
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
	int clerr = clSetKernelArg(M_tm_sitediagonal_minus_AND_gamma5_eo, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal_minus_AND_gamma5_eo, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal_minus_AND_gamma5_eo, 2, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(M_tm_sitediagonal_minus_AND_gamma5_eo , gs2, ls2);
}

void Opencl_Module_Fermions::M_tm_inverse_sitediagonal_minus_device(cl_mem in, cl_mem out, hmc_float mubar)
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
	int clerr = clSetKernelArg(M_tm_inverse_sitediagonal_minus, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_inverse_sitediagonal_minus, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_inverse_sitediagonal_minus, 2, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( M_tm_inverse_sitediagonal_minus, gs2, ls2);
}

void Opencl_Module_Fermions::M_tm_sitediagonal_minus_device(cl_mem in, cl_mem out, hmc_float mubar)
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
	int clerr = clSetKernelArg(M_tm_sitediagonal_minus, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal_minus, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(M_tm_sitediagonal_minus, 2, sizeof(hmc_float), &mubar_tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(M_tm_sitediagonal_minus , gs2, ls2);
}

int Opencl_Module_Fermions::bicgstab(const Matrix_Function & f, cl_mem inout, cl_mem source, cl_mem gf, hmc_float prec, hmc_float kappa, hmc_float mubar)
{
	//"save" version, with comments. this is called if "bicgstab_save" is choosen.
	if (get_parameters().get_solver() == meta::Inputparameters::bicgstab_save) {
		hmc_float resid;
		for(int iter = 0; iter < get_parameters().get_cgmax(); iter++) {
			if(iter % get_parameters().get_iter_refresh() == 0) {
				set_zero_spinorfield_device(clmem_v);
				set_zero_spinorfield_device(clmem_p);
				//initial r_n
				f(inout, clmem_rn, gf, kappa, mubar);
				saxpy_device(clmem_rn, source, clmem_one, clmem_rn);
				//rhat = r_n
				copy_buffer_on_device(clmem_rn, clmem_rhat, meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor));
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
				//print the last residuum
				logger.fatal() << "\t\t\tsolver stuck at resid:\t" << resid;
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
			f(clmem_p, clmem_v, gf, kappa, mubar);
			//tmp1 = (rhat, v)
			set_complex_to_scalar_product_device(clmem_rhat, clmem_v, clmem_tmp1);
			//alpha = rho/tmp1 = (..)/(rhat, v)
			set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha);
			//s = - alpha * v - r_n
			saxpy_device(clmem_v, clmem_rn, clmem_alpha, clmem_s);
			//t = A s
			f(clmem_s, clmem_t, gf, kappa, mubar);
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
			get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));

			logger.debug() << "resid: " << resid;
			//test if resid is NAN
			if(resid != resid) {
				logger.fatal() << "\tNAN occured in bicgstab!";
				return -iter;
			}
			if(resid < prec) {
				//aux = A inout
				f(inout, clmem_aux, gf, kappa, mubar);
				//aux = -aux + source
				saxpy_device(clmem_aux, source, clmem_one, clmem_aux);
				//trueresid = (aux, aux)
				set_float_to_global_squarenorm_device(clmem_aux, clmem_trueresid);
				hmc_float trueresid;
				get_buffer_from_device(clmem_trueresid, &trueresid, sizeof(hmc_float));
				logger.debug() << "\tsolver converged! true resid:\t" << trueresid;
				if(trueresid < prec)
					return iter;
			}
		}
		return -1;
	}
	//version with different structure than "save" one, similar to tmlqcd. This should be the default bicgstab.
	//  In particular this version does not perform the check if the "real" residuum is sufficiently small!
	else if (get_parameters().get_solver() == meta::Inputparameters::bicgstab) {
		hmc_float resid;
		for(int iter = 0; iter < get_parameters().get_cgmax(); iter++) {
			if(iter % get_parameters().get_iter_refresh() == 0) {
				//initial r_n, saved in p
				f(inout, clmem_rn, gf, kappa, mubar);
				saxpy_device(clmem_rn, source, clmem_one, clmem_p);
				//rhat = p
				copy_buffer_on_device(clmem_p, clmem_rhat, meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor));
				//r_n = p
				copy_buffer_on_device(clmem_p, clmem_rn, meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor));
				//rho = (rhat, rn)
				set_complex_to_scalar_product_device(clmem_rhat, clmem_rn, clmem_rho);
			}
			//resid = (rn,rn)
			set_float_to_global_squarenorm_device(clmem_rn, clmem_resid);
			get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));
			//test if resid is NAN
			if(resid != resid) {
				logger.fatal() << "\tNAN occured in bicgstab!";
				return -iter;
			}
			if(resid < prec) {
				return iter;
			}
			//v = A*p
			f(clmem_p, clmem_v, gf, kappa, mubar);
			//tmp1 = (rhat, v)
			set_complex_to_scalar_product_device(clmem_rhat, clmem_v, clmem_tmp1);
			//alpha = rho/tmp1 = (rhat, rn)/(rhat, v)
			set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha);
			//s = - alpha * v - r_n
			saxpy_device(clmem_v, clmem_rn, clmem_alpha, clmem_s);
			//t = A s
			f(clmem_s, clmem_t, gf, kappa, mubar);
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
				//print the last residuum
				logger.fatal() << "\t\t\tsolver stuck at resid:\t" << resid;
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
	logger.fatal() << "This line should never be reached! Program logic is inconsistent!";
	throw std::logic_error("Somebody goofed up - the solver you selected is not implemented.");
}

int Opencl_Module_Fermions::bicgstab_eo(const Matrix_Function & f, cl_mem inout, cl_mem source, cl_mem gf, hmc_float prec, hmc_float kappa, hmc_float mubar)
{
	cl_int clerr = CL_SUCCESS;

	//"save" version, with comments. this is called if "bicgstab_save" is choosen.
	if (get_parameters().get_solver() == meta::Inputparameters::bicgstab_save) {
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
		int cgmax = get_parameters().get_cgmax();
		for(int iter = 0; iter < cgmax; iter++) {
			if(iter % get_parameters().get_iter_refresh() == 0) {
				set_zero_spinorfield_eoprec_device(clmem_v_eo);
				set_zero_spinorfield_eoprec_device(clmem_p_eo);

				f(inout, clmem_rn_eo, gf, kappa, mubar);

				saxpy_eoprec_device(clmem_rn_eo, source, clmem_one, clmem_rn_eo);

				copy_buffer_on_device(clmem_rn_eo, clmem_rhat_eo, get_eoprec_spinorfield_buffer_size());

				copy_buffer_on_device(clmem_one, clmem_alpha, sizeof(hmc_complex));
				copy_buffer_on_device(clmem_one, clmem_omega, sizeof(hmc_complex));
				copy_buffer_on_device(clmem_one, clmem_rho, sizeof(hmc_complex));
			}
			set_complex_to_scalar_product_eoprec_device(clmem_rhat_eo, clmem_rn_eo, clmem_rho_next);
			//check if algorithm is stuck
			hmc_complex check;
			get_buffer_from_device(clmem_rho_next, &check, sizeof(hmc_complex));
			if(abs(check.re) < 1e-25 && abs(check.im) < 1e-25 ) {
				//print the last residuum
				logger.fatal() << "\t\t\tsolver stuck at resid:\t" << resid;
				return -iter;
			}
			set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1);
			copy_buffer_on_device(clmem_rho_next, clmem_rho, sizeof(hmc_complex));
			set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2);
			set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta);

			set_complex_to_product_device(clmem_beta, clmem_omega, clmem_tmp1);
			set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2);
			saxsbypz_eoprec_device(clmem_p_eo, clmem_v_eo, clmem_rn_eo, clmem_beta, clmem_tmp2, clmem_p_eo);

			f(clmem_p_eo, clmem_v_eo, gf, kappa, mubar);

			set_complex_to_scalar_product_eoprec_device(clmem_rhat_eo, clmem_v_eo, clmem_tmp1);
			set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha);

			saxpy_eoprec_device(clmem_v_eo, clmem_rn_eo, clmem_alpha, clmem_s_eo);

			f(clmem_s_eo, clmem_t_eo, gf, kappa, mubar);

			set_complex_to_scalar_product_eoprec_device(clmem_t_eo, clmem_s_eo, clmem_tmp1);
			//!!CP: can this also be global_squarenorm??
			set_complex_to_scalar_product_eoprec_device(clmem_t_eo, clmem_t_eo, clmem_tmp2);
			set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega);

			saxpy_eoprec_device(clmem_t_eo, clmem_s_eo, clmem_omega, clmem_rn_eo);

			saxsbypz_eoprec_device(clmem_p_eo, clmem_s_eo, inout, clmem_alpha, clmem_omega, inout);

			set_float_to_global_squarenorm_eoprec_device(clmem_rn_eo, clmem_resid);
			get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));

			logger.debug() << "resid: " << resid;
			//test if resid is NAN
			if(resid != resid) {
				logger.fatal() << "\tNAN occured in bicgstab_eo!";
				return -iter;
			}
			if(resid < prec) {
				++retests;

				f(inout, clmem_aux_eo, gf, kappa, mubar);
				saxpy_eoprec_device(clmem_aux_eo, source, clmem_one, clmem_aux_eo);

				set_float_to_global_squarenorm_eoprec_device(clmem_aux_eo, clmem_trueresid);
				get_buffer_from_device(clmem_trueresid, &trueresid, sizeof(hmc_float));
				logger.debug() << "\ttrueresiduum:\t" << trueresid;
				if(trueresid < prec) {
					// report on performance
					if(logger.beInfo()) {
						// we are always synchroneous here, as we had to recieve the residium from the device
						uint64_t duration = timer.getTime();

						// calculate flops
						unsigned refreshs = iter / get_parameters().get_iter_refresh() + 1;
						cl_ulong mf_flops = f.get_Flops();

						cl_ulong total_flops = 4 * get_flop_size("scalar_product_eoprec") + 4 * get_flop_size("ratio") + 3 * get_flop_size("product") + 2 * get_flop_size("saxsbypz_eoprec") + 2 * mf_flops + 2 * get_flop_size("saxpy_eoprec") + get_flop_size("global_squarenorm_eoprec");
						total_flops *= iter;

						total_flops += refreshs * (mf_flops + get_flop_size("saxpy_eoprec"));

						total_flops += retests * (mf_flops + get_flop_size("saxpy_eoprec") + get_flop_size("global_squarenorm_eoprec"));

						// report performanc
						logger.info() << "BiCGstab_save completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f) << " Gflops. Performed " << iter << " iterations";
					}

					// we are done here
					return iter;
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
		for(int iter = 0; iter < get_parameters().get_cgmax(); iter++) {
			if(iter % get_parameters().get_iter_refresh() == 0) {
				//initial r_n, saved in p
				f(inout, clmem_rn_eo, gf, kappa, mubar);
				saxpy_eoprec_device(clmem_rn_eo, source, clmem_one, clmem_p_eo);
				//rhat = p
				copy_buffer_on_device(clmem_p_eo, clmem_rhat_eo, get_eoprec_spinorfield_buffer_size());
				//r_n = p
				copy_buffer_on_device(clmem_p_eo, clmem_rn_eo, get_eoprec_spinorfield_buffer_size());
				//rho = (rhat, rn)
				set_complex_to_scalar_product_eoprec_device(clmem_rhat_eo, clmem_rn_eo, clmem_rho);
			}
			//resid = (rn,rn)
			set_float_to_global_squarenorm_eoprec_device(clmem_rn_eo, clmem_resid);
			hmc_float resid;
			get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));

			logger.debug() << "resid: " << resid;
			//test if resid is NAN
			if(resid != resid) {
				logger.fatal() << "\tNAN occured in bicgstab_eo!";
				return -iter;
			}
			if(resid < prec) {
				// report on performance
				if(logger.beInfo()) {
					// we are always synchroneous here, as we had to recieve the residium from the device
					uint64_t duration = timer.getTime();

					// calculate flops
					unsigned refreshs = iter / get_parameters().get_iter_refresh() + 1;
					cl_ulong mf_flops = f.get_Flops();

					cl_ulong total_flops = get_flop_size("global_squarenorm_eoprec") + 2 * mf_flops + 4 * get_flop_size("scalar_product_eoprec") + 4 * get_flop_size("ratio") + 2 * get_flop_size("saxpy_eoprec") + 2 * get_flop_size("saxsbypz_eoprec") + 3 * get_flop_size("product");
					total_flops *= iter;

					total_flops += refreshs * (mf_flops + get_flop_size("saxpy_eoprec") + get_flop_size("scalar_product_eoprec"));

					// report performanc
					logger.info() << "BiCGstab completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f) << " Gflops. Performed " << iter << " iterations";
				}

				// we are done here
				return iter;
			}
			//v = A*p
			f(clmem_p_eo, clmem_v_eo, gf, kappa, mubar);
			clerr = clFlush(get_queue());
			if(clerr) {
				throw Opencl_Error(clerr, "Failed to flush command queue");
			}
			//tmp1 = (rhat, v)
			set_complex_to_scalar_product_eoprec_device(clmem_rhat_eo, clmem_v_eo, clmem_tmp1);
			//alpha = rho/tmp1 = (rhat, rn)/(rhat, v)
			set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha);
			//s = - alpha * v - r_n
			saxpy_eoprec_device(clmem_v_eo, clmem_rn_eo, clmem_alpha, clmem_s_eo);
			//t = A s
			f(clmem_s_eo, clmem_t_eo, gf, kappa, mubar);
			//tmp1 = (t, s)
			set_complex_to_scalar_product_eoprec_device(clmem_t_eo, clmem_s_eo, clmem_tmp1);
			//!!CP: this can also be global_squarenorm, but one needs a complex number here
			//tmp2 = (t,t)
			set_complex_to_scalar_product_eoprec_device(clmem_t_eo, clmem_t_eo, clmem_tmp2);
			//omega = tmp1/tmp2 = (t,s)/(t,t)
			set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega);
			//inout = alpha*p + omega * s + inout
			saxsbypz_eoprec_device(clmem_p_eo, clmem_s_eo, inout, clmem_alpha, clmem_omega, inout);
			//r_n = - omega*t - s
			saxpy_eoprec_device(clmem_t_eo, clmem_s_eo, clmem_omega, clmem_rn_eo);
			//rho_next = (rhat, rn)
			set_complex_to_scalar_product_eoprec_device(clmem_rhat_eo, clmem_rn_eo, clmem_rho_next);
			//check if algorithm is stuck
			hmc_complex check;
			get_buffer_from_device(clmem_rho_next, &check, sizeof(hmc_complex));
			if(abs(check.re) < 1e-25 && abs(check.im) < 1e-25 ) {
				//print the last residuum
				logger.fatal() << "\t\t\tsolver stuck at resid:\t" << resid;
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
			saxsbypz_eoprec_device(clmem_p_eo, clmem_v_eo, clmem_rn_eo, clmem_beta, clmem_tmp2, clmem_p_eo);
			clerr = clFlush(get_queue());
			if(clerr) {
				throw Opencl_Error(clerr, "Failed to flush command queue");
			}
			//rho_next = rho
			copy_buffer_on_device(clmem_rho_next, clmem_rho, sizeof(hmc_complex));
		}
		return -1;
	}
}

int Opencl_Module_Fermions::cg(const Matrix_Function & f, cl_mem inout, cl_mem source, cl_mem gf, hmc_float prec, hmc_float kappa, hmc_float mubar)
{
	//CP: here I do not use clmem_rnhat anymore and saved one scalar_product (omega)
	//NOTE: here, most of the complex numbers may also be just hmc_floats. However, for this one would need some add. functions...
	for(int iter = 0; iter < get_parameters().get_cgmax(); iter ++) {
		if(iter % get_parameters().get_iter_refresh() == 0) {
			//rn = A*inout
			f(inout, clmem_rn, gf, kappa, mubar);
			//rn = source - A*inout
			saxpy_device(clmem_rn, source, clmem_one, clmem_rn);
			//p = rn
			copy_buffer_on_device(clmem_rn, clmem_p, meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor));
			//omega = (rn,rn)
			set_complex_to_scalar_product_device(clmem_rn, clmem_rn, clmem_omega);
		} else {
			//update omega
			copy_buffer_on_device(clmem_rho_next, clmem_omega, sizeof(hmc_complex));
		}
		//v = A pn
		f(clmem_p, clmem_v, gf, kappa, mubar);
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
		//test if resid is NAN
		if(resid != resid) {
			logger.fatal() << "\tNAN occured in cg!";
			return -iter;
		}
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

int Opencl_Module_Fermions::cg_eo(const Matrix_Function & f, cl_mem inout, cl_mem source, cl_mem gf, hmc_float prec, hmc_float kappa, hmc_float mubar)
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
	for(int iter = 0; iter < get_parameters().get_cgmax(); iter ++) {
		if(iter % get_parameters().get_iter_refresh() == 0) {
			//rn = A*inout
			f(inout, clmem_rn_eo, gf, kappa, mubar);
			//rn = source - A*inout
			saxpy_eoprec_device(clmem_rn_eo, source, clmem_one, clmem_rn_eo);
			//p = rn
			copy_buffer_on_device(clmem_rn_eo, clmem_p_eo, get_eoprec_spinorfield_buffer_size());
			//omega = (rn,rn)
			set_complex_to_scalar_product_eoprec_device(clmem_rn_eo, clmem_rn_eo, clmem_omega);
		} else {
			//update omega
			copy_buffer_on_device(clmem_rho_next, clmem_omega, sizeof(hmc_complex));
		}
		//v = A pn
		f(clmem_p_eo, clmem_v_eo, gf, kappa, mubar);

		//alpha = (rn, rn)/(pn, Apn) --> alpha = omega/rho
		set_complex_to_scalar_product_eoprec_device(clmem_p_eo, clmem_v_eo, clmem_rho);
		set_complex_to_ratio_device(clmem_omega, clmem_rho, clmem_alpha);
		set_complex_to_product_device(clmem_alpha, clmem_minusone, clmem_tmp1);

		//xn+1 = xn + alpha*p = xn - tmp1*p = xn - (-tmp1)*p
		saxpy_eoprec_device(clmem_p_eo, inout, clmem_tmp1, inout);
		//switch between original version and kernel merged one
		if(get_parameters().get_use_merge_kernels_spinor() == false) {
			//rn+1 = rn - alpha*v -> rhat
			saxpy_eoprec_device(clmem_v_eo, clmem_rn_eo, clmem_alpha, clmem_rn_eo);

			//calc residuum
			//NOTE: for beta one needs a complex number at the moment, therefore, this is done with "rho_next" instead of "resid"
			set_complex_to_scalar_product_eoprec_device(clmem_rn_eo, clmem_rn_eo, clmem_rho_next);
		} else {
			//merge two calls:
			//rn+1 = rn - alpha*v -> rhat
			//and
			//rho_next = |rhat|^2
			//rho_next is a complex number, set its imag to zero
			hmc_complex zero = hmc_complex_zero;
			int clerr = clEnqueueWriteBuffer(get_queue(), clmem_rho_next, CL_TRUE, 0, sizeof(hmc_complex), &zero, 0, 0, NULL);
			if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);
			saxpy_AND_squarenorm_eo_device(clmem_v_eo, clmem_rn_eo, clmem_alpha, clmem_rn_eo, clmem_rho_next);
		}
		hmc_float resid;
		get_buffer_from_device(clmem_rho_next, &resid, sizeof(hmc_float));
		//this is the orig. call
		//set_float_to_global_squarenorm_device(clmem_rn, clmem_resid);
		//get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));

		logger.debug() << "resid: " << resid;
		//test if resid is NAN
		if(resid != resid) {
			logger.fatal() << "\tNAN occured in cg_eo!";
			return -iter;
		}
		if(resid < prec) {
			// report on performance
			if(logger.beInfo()) {
				// we are always synchroneous here, as we had to recieve the residium from the device
				uint64_t duration = timer.getTime();

				// calculate flops
				unsigned refreshs = iter / get_parameters().get_iter_refresh() + 1;
				cl_ulong mf_flops = f.get_Flops();

				cl_ulong total_flops = mf_flops + 3 * get_flop_size("scalar_product_eoprec") + 2 * get_flop_size("ratio") + 2 * get_flop_size("product") + 3 * get_flop_size("saxpy_eoprec");
				total_flops *= iter;

				total_flops += refreshs * (mf_flops + get_flop_size("saxpy_eoprec") + get_flop_size("scalar_product_eoprec"));

				// report performanc
				logger.info() << "CG completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f) << " Gflops. Performed " << iter << " iterations";
			}

			return iter;
		}

		//beta = (rn+1, rn+1)/(rn, rn) --> alpha = rho_next/omega
		set_complex_to_ratio_device(clmem_rho_next, clmem_omega, clmem_beta);

		//pn+1 = rn+1 + beta*pn
		set_complex_to_product_device(clmem_beta, clmem_minusone, clmem_tmp2);
		saxpy_eoprec_device(clmem_p_eo, clmem_rn_eo, clmem_tmp2, clmem_p_eo);
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

	if(get_parameters().get_profile_solver() ) (*solvertimer).reset();
	if(get_parameters().get_use_eo() == true) {
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
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			//in this case, the diagonal matrix is just 1 and falls away.
			dslash_eo_device(clmem_source_odd, clmem_tmp_eo_1, gf, EVEN);
			saxpy_eoprec_device(clmem_source_even, clmem_tmp_eo_1, clmem_one, clmem_source_even);
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			M_tm_inverse_sitediagonal_device(clmem_source_odd, clmem_tmp_eo_1);
			dslash_eo_device(clmem_tmp_eo_1, clmem_tmp_eo_2, gf, EVEN);
			saxpy_eoprec_device(clmem_source_even, clmem_tmp_eo_2, clmem_one, clmem_source_even);
		}

		//Trial solution
		///@todo this should go into a more general function
		this->set_eoprec_spinorfield_cold_device(this->get_clmem_inout_eo());
		logger.debug() << "start eoprec-inversion";
		//make sure buffer are initialised
		this->fill_solver_buffers();
		//even solution
		if(get_parameters().get_solver() == meta::Inputparameters::cg)
			converged = cg_eo(f, clmem_inout_eo, clmem_source_even, gf, get_parameters().get_solver_prec());
		else
			converged = bicgstab_eo(f, this->get_clmem_inout_eo(), clmem_source_even, gf, get_parameters().get_solver_prec());

		//release solver buffers
		this->clear_solver_buffers();

		//odd solution
		/** The odd solution is obtained from the even one according to:
		 *  x_o = M_inv D x_e - M_inv b_o
		 */
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			//in this case, the diagonal matrix is just 1 and falls away.
			dslash_eo_device(clmem_inout_eo, clmem_tmp_eo_1, gf, ODD);
			saxpy_eoprec_device(clmem_tmp_eo_1, clmem_source_odd, clmem_one, clmem_tmp_eo_1);
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			dslash_eo_device(clmem_inout_eo, clmem_tmp_eo_2, gf, ODD);
			M_tm_inverse_sitediagonal_device(clmem_tmp_eo_2, clmem_tmp_eo_1);
			M_tm_inverse_sitediagonal_device(clmem_source_odd, clmem_tmp_eo_2);
			saxpy_eoprec_device(clmem_tmp_eo_1, clmem_tmp_eo_2, clmem_one, clmem_tmp_eo_1);
		}
		//CP: whole solution
		//CP: suppose the even sol is saved in inout_eoprec, the odd one in clmem_tmp_eo_1
		convert_from_eoprec_device(clmem_inout_eo, clmem_tmp_eo_1, inout);
	} else {
		//Trial solution
		///@todo this should go into a more general function
		this->set_spinorfield_cold_device(inout);

		if(get_parameters().get_solver() == meta::Inputparameters::cg)
			converged = cg(f, inout, source, gf, get_parameters().get_solver_prec());
		else
			converged = bicgstab(f, inout, source, gf, get_parameters().get_solver_prec());
	}
	if(get_parameters().get_profile_solver() ) {
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

cl_mem Opencl_Module_Fermions::get_clmem_inout_eo()
{
	return clmem_inout_eo;
}

cl_mem Opencl_Module_Fermions::get_clmem_tmp_eo_1()
{
	return clmem_tmp_eo_1;
}

cl_mem Opencl_Module_Fermions::get_clmem_tmp_eo_2()
{
	return clmem_tmp_eo_2;
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

size_t Opencl_Module_Fermions::get_read_write_size(const std::string& in) const
{
	size_t result = Opencl_Module_Spinors::get_read_write_size(in);
	if (result != 0) return result;
	//Depending on the compile-options, one has different sizes...
	size_t D = meta::get_float_size(parameters);
	//this returns the number of entries in an su3-matrix
	size_t R = meta::get_mat_size(parameters);
	//this is the number of spinors in the system (or number of sites)
	size_t S = meta::get_spinorfieldsize(get_parameters());
	size_t Seo = meta::get_eoprec_spinorfieldsize(get_parameters());
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

uint64_t Opencl_Module_Fermions::get_flop_size(const std::string& in) const
{
	uint64_t result = Opencl_Module_Spinors::get_flop_size(in);
	if (result != 0) return result;
	size_t S = meta::get_spinorfieldsize(get_parameters());
	size_t Seo = meta::get_eoprec_spinorfieldsize(get_parameters());
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

void Opencl_Module_Fermions::print_profiling(const std::string& filename, int number)
{
	Opencl_Module_Spinors::print_profiling(filename, number);
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
