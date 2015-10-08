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

#include "fermionmatrix.hpp"

#include "../../hardware/code/spinors.hpp"

/**
* Implementations of fermion matrices
*/

bool physics::fermionmatrix::Fermionmatrix_basic::is_hermitian() const noexcept
{
	return _is_hermitian;
}

hmc_float physics::fermionmatrix::Fermionmatrix_basic::get_kappa() const noexcept
{
	return kappa;
}

hmc_float physics::fermionmatrix::Fermionmatrix_basic::get_mubar() const noexcept
{
	return mubar;
}

const hardware::System& physics::fermionmatrix::Fermionmatrix_basic::get_system() const noexcept
{
	return system;
}

void physics::fermionmatrix::M::operator()(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in) const
{
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			//in the pure Wilson case there is just one fermionmatrix
			M_wilson(out, gf, in, get_kappa());
			break;
		case meta::action::twistedmass:
			M_tm_plus(out, gf, in, get_kappa(), get_mubar());
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
}
cl_ulong physics::fermionmatrix::M::get_flops() const
{

	auto devices = get_system().get_devices();
	auto fermion_code = devices[0]->get_fermion_code();
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			return fermion_code->get_flop_size("M_wilson");
		case meta::action::twistedmass:
			return fermion_code->get_flop_size("M_tm_plus");
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
}
cl_ulong physics::fermionmatrix::M::get_read_write_size() const
{

	auto devices = get_system().get_devices();
	auto fermion_code = devices[0]->get_fermion_code();
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			return fermion_code->get_read_write_size("M_wilson");
		case meta::action::twistedmass:
			return fermion_code->get_read_write_size("M_tm_plus");
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
}

void physics::fermionmatrix::Qplus::operator()(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in) const
{
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			//in the pure Wilson case there is just one fermionmatrix
			M_wilson(out, gf, in, get_kappa());
			break;
		case meta::action::twistedmass:
			M_tm_plus(out, gf, in, get_kappa(), get_mubar());
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
	out->gamma5();
}
cl_ulong physics::fermionmatrix::Qplus::get_flops() const
{
	auto devices = get_system().get_devices();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res;
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			res = fermion_code->get_flop_size("M_wilson");
			break;
		case meta::action::twistedmass:
			res = fermion_code->get_flop_size("M_tm_plus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
	res += fermion_code->get_flop_size("gamma5");
	return res;
}
cl_ulong physics::fermionmatrix::Qplus::get_read_write_size() const
{
	auto devices = get_system().get_devices();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res;
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			res = fermion_code->get_read_write_size("M_wilson");
			break;
		case meta::action::twistedmass:
			res = fermion_code->get_read_write_size("M_tm_plus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
	res += fermion_code->get_read_write_size("gamma5");
	return res;
}
void physics::fermionmatrix::Qminus::operator()(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in) const
{
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			//in the pure Wilson case there is just one fermionmatrix
			M_wilson(out, gf, in, get_kappa());
			break;
		case meta::action::twistedmass:
			M_tm_minus(out, gf, in, get_kappa(), get_mubar());
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
	out->gamma5();
}
cl_ulong physics::fermionmatrix::Qminus::get_flops() const
{
	auto devices = get_system().get_devices();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res;
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			res = fermion_code->get_flop_size("M_wilson");
			break;
		case meta::action::twistedmass:
			res = fermion_code->get_flop_size("M_tm_minus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
	res += fermion_code->get_flop_size("gamma5");
	return res;
}
cl_ulong physics::fermionmatrix::Qminus::get_read_write_size() const
{
	auto devices = get_system().get_devices();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res;
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			res = fermion_code->get_read_write_size("M_wilson");
			break;
		case meta::action::twistedmass:
			res = fermion_code->get_read_write_size("M_tm_minus");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
	res += fermion_code->get_read_write_size("gamma5");
	return res;
}
void physics::fermionmatrix::QplusQminus::operator()(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in) const
{
	q_minus(&tmp, gf, in);
	q_plus(out, gf, tmp);
}
cl_ulong physics::fermionmatrix::QplusQminus::get_flops() const
{
	return q_minus.get_flops() + q_plus.get_flops();
}
cl_ulong physics::fermionmatrix::QplusQminus::get_read_write_size() const
{
	return q_minus.get_read_write_size() + q_plus.get_read_write_size();
}
void physics::fermionmatrix::Aee::operator()(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const
{
	/**
	 * This is the even-odd preconditioned fermion matrix with the
	 * non-trivial inversion on the even sites (see DeGran/DeTar p. 174).
	 * If one has fermionmatrix
	 *  M = R + D,
	 * then Aee is:
	 * Aee = R_e - D_eo R_o_inv D_oe
	 */

	/** @todo The local creation of the temporary field is known to cause performance problems... */


	hmc_float kappa = get_kappa();
	hmc_float mubar = get_mubar();

	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			//in this case, the diagonal matrix is just 1 and falls away.
			dslash(&tmp, gf, in, ODD, kappa);
			dslash(out, gf, tmp, EVEN, kappa);
			saxpy(out, {1., 0.}, *out, in);
			break;
		case meta::action::twistedmass:
			dslash(&tmp, gf, in, ODD, kappa);
			M_tm_inverse_sitediagonal(&tmp2, tmp, mubar);
			dslash(out, gf, tmp2, EVEN, kappa);
			M_tm_sitediagonal(&tmp, in, mubar);
			saxpy(out, {1., 0.}, *out, tmp);
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
}
cl_ulong physics::fermionmatrix::Aee::get_flops() const
{
	auto devices = get_system().get_devices();
	auto spinor_code = devices[0]->get_spinor_code();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res;
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			res = 2 * fermion_code->get_flop_size("dslash_eo");
			res += spinor_code->get_flop_size("saxpy_eoprec");
			break;
		case meta::action::twistedmass:
			res = 2 * fermion_code->get_flop_size("dslash_eo");
			res += fermion_code->get_flop_size("M_tm_inverse_sitediagonal");
			res += fermion_code->get_flop_size("M_tm_sitediagonal");
			res += spinor_code->get_flop_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
	logger.trace() << "Aee flops: " << res;
	return res;
}
cl_ulong physics::fermionmatrix::Aee::get_read_write_size() const
{
	auto devices = get_system().get_devices();
	auto spinor_code = devices[0]->get_spinor_code();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res;
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			res = 2 * fermion_code->get_read_write_size("dslash_eo");
			res += spinor_code->get_read_write_size("saxpy_eoprec");
			break;
		case meta::action::twistedmass:
		        res = 2 * fermion_code->get_read_write_size("dslash_eo");
			res += fermion_code->get_read_write_size("M_tm_inverse_sitediagonal");
			res += fermion_code->get_read_write_size("M_tm_sitediagonal");
			res += spinor_code->get_read_write_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
	logger.trace() << "Aee read-write size: " << res;
	return res;
}
void physics::fermionmatrix::Aee_AND_gamma5_eo::operator()(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const
{
	/**
	 * This is the even-odd preconditioned fermion matrix with the
	 * non-trivial inversion on the even sites (see DeGran/DeTar p. 174).
	 * If one has fermionmatrix
	 *  M = R + D,
	 * then Aee is:
	 * Aee = R_e - D_eo R_o_inv D_oe
	 */

	/** @todo The local creation of the temporary field is known to cause performance problems... */
	hmc_float kappa = get_kappa();
	hmc_float mubar = get_mubar();

	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			//in this case, the diagonal matrix is just 1 and falls away.
			dslash(&tmp, gf, in, ODD, kappa);
			dslash(out, gf, tmp, EVEN, kappa);
			saxpy_AND_gamma5_eo(out, {1., 0.}, *out, in);
			break;
		case meta::action::twistedmass:
			dslash(&tmp, gf, in, ODD, kappa);
			M_tm_inverse_sitediagonal(&tmp2, tmp, mubar);
			dslash(out, gf, tmp2, EVEN, kappa);
			M_tm_sitediagonal(&tmp, in, mubar);
			saxpy_AND_gamma5_eo(out, {1., 0.}, *out, tmp);
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
}
cl_ulong physics::fermionmatrix::Aee_AND_gamma5_eo::get_flops() const
{
	auto devices = get_system().get_devices();
	auto spinor_code = devices[0]->get_spinor_code();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res;
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			res = 2 * fermion_code->get_flop_size("dslash_eo");
			res += spinor_code->get_flop_size("saxpy_AND_gamma5_eo");
			break;
		case meta::action::twistedmass:
			res = 2 * fermion_code->get_flop_size("dslash_eo");
			res += fermion_code->get_flop_size("M_tm_inverse_sitediagonal");
			res += fermion_code->get_flop_size("M_tm_sitediagonal");
			res += spinor_code->get_flop_size("saxpy_AND_gamma5_eo");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
	logger.trace() << "Aee_AND_gamma5_eo flops: " << res;
	return res;
}
cl_ulong physics::fermionmatrix::Aee_AND_gamma5_eo::get_read_write_size() const
{
	auto devices = get_system().get_devices();
	auto spinor_code = devices[0]->get_spinor_code();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res;
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			res = 2 * fermion_code->get_read_write_size("dslash_eo");
			res += spinor_code->get_read_write_size("saxpy_AND_gamma5_eo");
			break;
		case meta::action::twistedmass:
		        res = 2 * fermion_code->get_read_write_size("dslash_eo");
			res += fermion_code->get_read_write_size("M_tm_inverse_sitediagonal");
			res += fermion_code->get_read_write_size("M_tm_sitediagonal");
			res += spinor_code->get_read_write_size("saxpy_AND_gamma5_eo");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
	logger.trace() << "Aee_AND_gamma5_eo read-write size: " << res;
	return res;
}
void physics::fermionmatrix::Aee_minus::operator()(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const
{
	/**
	 * This is the even-odd preconditioned fermion matrix with the
	 * non-trivial inversion on the even sites (see DeGran/DeTar p. 174).
	 * If one has fermionmatrix
	 *  M = R + D,
	 * then Aee is:
	 * Aee = R_e - D_eo R_o_inv D_oe
	 */

	/** @todo The local creation of the temporary field is known to cause performance problems... */
	hmc_float kappa = get_kappa();
	hmc_float mubar = get_mubar();

	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			//in this case, the diagonal matrix is just 1 and falls away.
			dslash(&tmp, gf, in, ODD, kappa);
			dslash(out, gf, tmp, EVEN, kappa);
			saxpy(out, {1., 0.}, *out, in);
			break;
		case meta::action::twistedmass:
			dslash(&tmp, gf, in, ODD, kappa);
			M_tm_inverse_sitediagonal_minus(&tmp2, tmp, mubar);
			dslash(out, gf, tmp2, EVEN, kappa);
			M_tm_sitediagonal_minus(&tmp, in, mubar);
			saxpy(out, {1., 0.}, *out, tmp);
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
}
cl_ulong physics::fermionmatrix::Aee_minus::get_flops() const
{
	auto devices = get_system().get_devices();
	auto spinor_code = devices[0]->get_spinor_code();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res;
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			res = 2 * fermion_code->get_flop_size("dslash_eo");
			res += spinor_code->get_flop_size("saxpy_eoprec");
			break;
		case meta::action::twistedmass:
			res = 2 * fermion_code->get_flop_size("dslash_eo");
			res += fermion_code->get_flop_size("M_tm_inverse_sitediagonal_minus");
			res += fermion_code->get_flop_size("M_tm_sitediagonal_minus");
			res += spinor_code->get_flop_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
	logger.trace() << "Aee_minus flops: " << res;
	return res;
}
cl_ulong physics::fermionmatrix::Aee_minus::get_read_write_size() const
{
	auto devices = get_system().get_devices();
	auto spinor_code = devices[0]->get_spinor_code();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res;
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			res = 2 * fermion_code->get_read_write_size("dslash_eo");
			res += spinor_code->get_read_write_size("saxpy_eoprec");
			break;
		case meta::action::twistedmass:
		        res = 2 * fermion_code->get_read_write_size("dslash_eo");
			res += fermion_code->get_read_write_size("M_tm_inverse_sitediagonal_minus");
			res += fermion_code->get_read_write_size("M_tm_sitediagonal_minus");
			res += spinor_code->get_read_write_size("saxpy_eoprec");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
	logger.trace() << "Aee_minus read-write size: " << res;
	return res;
}
void physics::fermionmatrix::Aee_minus_AND_gamma5_eo::operator()(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const
{
	/**
	 * This is the even-odd preconditioned fermion matrix with the
	 * non-trivial inversion on the even sites (see DeGran/DeTar p. 174).
	 * If one has fermionmatrix
	 *  M = R + D,
	 * then Aee is:
	 * Aee = R_e - D_eo R_o_inv D_oe
	 */

	/** @todo The local creation of the temporary field is known to cause performance problems... */
	hmc_float kappa = get_kappa();
	hmc_float mubar = get_mubar();

	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			//in this case, the diagonal matrix is just 1 and falls away.
			dslash(&tmp, gf, in, ODD, kappa);
			dslash(out, gf, tmp, EVEN, kappa);
			saxpy_AND_gamma5_eo(out, {1., 0.}, *out, in);
			break;
		case meta::action::twistedmass:
			dslash(&tmp, gf, in, ODD, kappa);
			M_tm_inverse_sitediagonal_minus(&tmp2, tmp, mubar);
			dslash(out, gf, tmp2, EVEN, kappa);
			M_tm_sitediagonal_minus(&tmp, in, mubar);
			saxpy_AND_gamma5_eo(out, {1., 0.}, *out, tmp);
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
}
cl_ulong physics::fermionmatrix::Aee_minus_AND_gamma5_eo::get_flops() const
{
	auto devices = get_system().get_devices();
	auto spinor_code = devices[0]->get_spinor_code();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res;
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			res = 2 * fermion_code->get_flop_size("dslash_eo");
			res += spinor_code->get_flop_size("saxpy_AND_gamma5_eo");
			break;
		case meta::action::twistedmass:
			res = 2 * fermion_code->get_flop_size("dslash_eo");
			res += fermion_code->get_flop_size("M_tm_inverse_sitediagonal_minus");
			res += fermion_code->get_flop_size("M_tm_sitediagonal_minus");
			res += spinor_code->get_flop_size("saxpy_AND_gamma5_eo");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
	logger.trace() << "Aee_minus_AND_gamma5_eo flops: " << res;
	return res;
}
cl_ulong physics::fermionmatrix::Aee_minus_AND_gamma5_eo::get_read_write_size() const
{
	auto devices = get_system().get_devices();
	auto spinor_code = devices[0]->get_spinor_code();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res;
	switch(fermionmatrixParametersInterface->getFermionicActionType()) {
		case meta::action::wilson:
			res = 2 * fermion_code->get_read_write_size("dslash_eo");
			res += spinor_code->get_read_write_size("saxpy_AND_gamma5_eo");
			break;
		case meta::action::twistedmass:
		        res = 2 * fermion_code->get_read_write_size("dslash_eo");
			res += fermion_code->get_read_write_size("M_tm_inverse_sitediagonal_minus");
			res += fermion_code->get_read_write_size("M_tm_sitediagonal_minus");
			res += spinor_code->get_read_write_size("saxpy_AND_gamma5_eo");
			break;
		default:
			throw Invalid_Parameters("Unkown fermion action!", "wilson or twistedmass", fermionmatrixParametersInterface->getFermionicActionType());
	}
	logger.trace() << "Aee_minus_AND_gamma5_eo read-write size: " << res;
	return res;
}
void physics::fermionmatrix::Qplus_eo::operator()(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const
{
	if(fermionmatrixParametersInterface->useMergedFermionicKernels() == false) {
		aee(out, gf, in);
		out->gamma5();
	} else {
		aee_AND_gamma5_eo(out, gf, in);
	}
}
cl_ulong physics::fermionmatrix::Qplus_eo::get_flops() const
{
	auto devices = get_system().get_devices();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res = aee.get_flops();
	res += fermion_code->get_flop_size("gamma5_eo");
	logger.trace() << "Qplus_eo flops: " << res;
	return res;
}
void physics::fermionmatrix::Qminus_eo::operator()(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const
{
	if(fermionmatrixParametersInterface->useMergedFermionicKernels() == false) {
		aee_minus(out, gf, in);
		out->gamma5();
	} else {
		aee_minus_AND_gamma5_eo(out, gf, in);
	}
}
cl_ulong physics::fermionmatrix::Qminus_eo::get_flops() const
{
	auto devices = get_system().get_devices();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res = aee_minus.get_flops();
	res += fermion_code->get_flop_size("gamma5_eo");
	logger.trace() << "Qminus_eo flops: " << res;
	return res;
}
cl_ulong physics::fermionmatrix::Qplus_eo::get_read_write_size() const
{
	auto devices = get_system().get_devices();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res = aee.get_read_write_size();
	res += fermion_code->get_read_write_size("gamma5_eo");
	logger.trace() << "Qplus_eo read-write size: " << res;
	return res;
}
cl_ulong physics::fermionmatrix::Qminus_eo::get_read_write_size() const
{
	auto devices = get_system().get_devices();
	auto fermion_code = devices[0]->get_fermion_code();

	cl_ulong res = aee_minus.get_read_write_size();
	res += fermion_code->get_read_write_size("gamma5_eo");
	logger.trace() << "Qminus_eo read-write size: " << res;
	return res;
}
void physics::fermionmatrix::QplusQminus_eo::operator()(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const
{
	q_minus(&tmp, gf, in);
	q_plus(out, gf, tmp);
}
cl_ulong physics::fermionmatrix::QplusQminus_eo::get_flops() const
{
	cl_ulong res = q_minus.get_flops() + q_plus.get_flops();
	logger.trace() << "QplusQminus_eo flops: " << res;
	return res;

}
cl_ulong physics::fermionmatrix::QplusQminus_eo::get_read_write_size() const
{
	cl_ulong res = q_minus.get_read_write_size() + q_plus.get_read_write_size();
	logger.trace() << "QplusQminus_eo read-write size: " << res;
	return res;

}
