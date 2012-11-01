#include "fermionmatrix.hpp"

bool physics::fermionmatrix::Fermionmatrix_basic::get_is_hermitian() const noexcept
{
  return is_hermitian;
}

void physics::fermionmatrix::M::operator()(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa, hmc_float mubar) const
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

