#include "fermionmatrix.hpp"

/**
 * Implementations of fermion matrices
 */

bool physics::fermionmatrix::Fermionmatrix_basic::get_is_hermitian() const noexcept
{
	return is_hermitian;
}

void physics::fermionmatrix::M::operator()(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf) const
{
	that->M(in, out, gf, kappa, mubar);
}
cl_ulong physics::fermionmatrix::M::get_Flops() const
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
cl_ulong physics::fermionmatrix::M::get_Bytes() const
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

void physics::fermionmatrix::Qplus::operator()(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf) const
{
	that->Qplus(in, out, gf, kappa, mubar);
}
cl_ulong physics::fermionmatrix::Qplus::get_Flops() const
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
cl_ulong physics::fermionmatrix::Qplus::get_Bytes() const
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

void physics::fermionmatrix::Qminus::operator()(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf) const
{
	that->Qminus(in, out, gf, kappa, mubar);
}
cl_ulong physics::fermionmatrix::Qminus::get_Flops() const
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
cl_ulong physics::fermionmatrix::Qminus::get_Bytes() const
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

void physics::fermionmatrix::QplusQminus::operator()(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf) const
{
	that->QplusQminus(in, out, gf, kappa, mubar);
}
cl_ulong physics::fermionmatrix::QplusQminus::get_Flops() const
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
cl_ulong physics::fermionmatrix::QplusQminus::get_Bytes() const
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

void physics::fermionmatrix::Aee::operator()(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf) const
{
	that->Aee(in, out, gf, kappa, mubar);
}
cl_ulong physics::fermionmatrix::Aee::get_Flops() const
{
	auto spinor_code = that->get_device()->get_spinor_code();

	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_flop_size("dslash_eo");
			res += spinor_code->get_flop_size("saxpy_eoprec");
			break;
		case meta::Inputparameters::twistedmass:
			res = 2 * that->get_flop_size("dslash_eo");
			res += that->get_flop_size("M_tm_inverse_sitediagonal");
			res += that->get_flop_size("M_tm_sitediagonal");
			res += spinor_code->get_flop_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	return res;
}
cl_ulong physics::fermionmatrix::Aee::get_Bytes() const
{
	auto spinor_code = that->get_device()->get_spinor_code();

	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_read_write_size("dslash_eo");
			res += spinor_code->get_read_write_size("saxpy_eoprec");
			break;
		case meta::Inputparameters::twistedmass:
			res = 2 * that->get_read_write_size("dslash_eo");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal");
			res += that->get_read_write_size("M_tm_sitediagonal");
			res += spinor_code->get_read_write_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	return res;
}

void physics::fermionmatrix::Qplus_eo::operator()(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf) const
{
	that->Qplus_eo(in, out, gf, kappa, mubar);
}
cl_ulong physics::fermionmatrix::Qplus_eo::get_Flops() const
{
	auto spinor_code = that->get_device()->get_spinor_code();

	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_flop_size("dslash_eo");
			res += spinor_code->get_flop_size("saxpy_eoprec");
			break;
		case meta::Inputparameters::twistedmass:
			res = 2 * that->get_flop_size("dslash_eo");
			res += that->get_flop_size("M_tm_inverse_sitediagonal");
			res += that->get_flop_size("M_tm_sitediagonal");
			res += spinor_code->get_flop_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += that->get_flop_size("gamma5_eo");
	return res;
}
cl_ulong physics::fermionmatrix::Qplus_eo::get_Bytes() const
{
	auto spinor_code = that->get_device()->get_spinor_code();

	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_read_write_size("dslash_eo");
			res += spinor_code->get_read_write_size("saxpy_eoprec");
			break;
		case meta::Inputparameters::twistedmass:
			res = 2 * that->get_read_write_size("dslash_eo");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal");
			res += that->get_read_write_size("M_tm_sitediagonal");
			res += spinor_code->get_read_write_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += that->get_flop_size("gamma5_eo");
	return res;
}

void physics::fermionmatrix::Qminus_eo::operator()(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf) const
{
	that->Qminus_eo(in, out, gf, kappa, mubar);
}
cl_ulong physics::fermionmatrix::Qminus_eo::get_Flops() const
{
	auto spinor_code = that->get_device()->get_spinor_code();

	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_flop_size("dslash_eo");
			res += spinor_code->get_flop_size("saxpy_eoprec");
			break;
		case meta::Inputparameters::twistedmass:
			res = 2 * that->get_flop_size("dslash_eo");
			res += that->get_flop_size("M_tm_inverse_sitediagonal_minus");
			res += that->get_flop_size("M_tm_sitediagonal_minus");
			res += spinor_code->get_flop_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += that->get_flop_size("gamma5_eo");
	return res;
}
cl_ulong physics::fermionmatrix::Qminus_eo::get_Bytes() const
{
	auto spinor_code = that->get_device()->get_spinor_code();

	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_read_write_size("dslash_eo");
			res += spinor_code->get_read_write_size("saxpy_eoprec");
			break;
		case meta::Inputparameters::twistedmass:
			res = 2 * that->get_read_write_size("dslash_eo");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal_minus");
			res += that->get_read_write_size("M_tm_sitediagonal_minus");
			res += spinor_code->get_read_write_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += that->get_flop_size("gamma5_eo");
	return res;
}

void physics::fermionmatrix::QplusQminus_eo::operator()(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf) const
{
	that->QplusQminus_eo(in, out, gf, kappa, mubar);
}
cl_ulong physics::fermionmatrix::QplusQminus_eo::get_Flops() const
{
	auto spinor_code = that->get_device()->get_spinor_code();

	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_flop_size("dslash_eo");
			res += spinor_code->get_flop_size("saxpy_eoprec");
			res *= 2;
			break;
		case meta::Inputparameters::twistedmass:
			res = 4 * that->get_flop_size("dslash_eo");
			res += that->get_flop_size("M_tm_inverse_sitediagonal");
			res += that->get_flop_size("M_tm_sitediagonal");
			res += that->get_flop_size("M_tm_inverse_sitediagonal_minus");
			res += that->get_flop_size("M_tm_sitediagonal_minus");
			res += 2 * spinor_code->get_flop_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += 2 * that->get_flop_size("gamma5_eo");
	return res;
}
cl_ulong physics::fermionmatrix::QplusQminus_eo::get_Bytes() const
{
	auto spinor_code = that->get_device()->get_spinor_code();

	cl_ulong res;
	switch(that->get_parameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			res = 2 * that->get_read_write_size("dslash_eo");
			res += spinor_code->get_read_write_size("saxpy_eoprec");
			res *= 2;
			break;
		case meta::Inputparameters::twistedmass:
			res = 4 * that->get_read_write_size("dslash_eo");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal");
			res += that->get_read_write_size("M_tm_sitediagonal");
			res += that->get_read_write_size("M_tm_inverse_sitediagonal_minus");
			res += that->get_read_write_size("M_tm_sitediagonal_minus");
			res += 2 * spinor_code->get_read_write_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or meta::Inputparameters::twistedmass", that->get_parameters().get_fermact());
	}
	res += 2 * that->get_flop_size("gamma5_eo");
	return res;
}
