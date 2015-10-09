/** @file
 * Implementation of staggered fermionmatrix classes methods.
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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

#include "fermionmatrix_stagg.hpp"
#include "../../hardware/code/spinors_staggered.hpp"

//Generic basic class
bool physics::fermionmatrix::Fermionmatrix_stagg_basic::is_hermitian() const noexcept
{
	return _is_hermitian;
}

//Even if get_mass is a pure virtual method, the base class can still have an implementation
//of such a method that can be explicitly called with the scope resolution operator.
//In this way, derived classes that must return the value of the mass can call this default
//implementation, those that do not depend on the mass will throw an exception.
//REMARK: = 0 means derived classes must provide an implementation,
//        not that the base class can not provide an implementation.
hmc_float physics::fermionmatrix::Fermionmatrix_stagg_basic::get_mass() const
{
	return mass;
}

const hardware::System& physics::fermionmatrix::Fermionmatrix_stagg_basic::get_system() const noexcept
{
	return system;
}

//Class D_KS_eo
void physics::fermionmatrix::D_KS_eo::operator()(const physics::lattices::Staggeredfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Staggeredfield_eo& in) const
{
	DKS_eo(out, gf, in, evenodd);
}

cl_ulong physics::fermionmatrix::D_KS_eo::get_flops() const
{
	const hardware::System& system = get_system();
	auto devices = system.get_devices();
	auto fermion_code = devices[0]->get_fermion_staggered_code();
	
	return fermion_code->get_flop_size("D_KS_eo");
}

hmc_float physics::fermionmatrix::D_KS_eo::get_mass() const
{
	throw Print_Error_Message("Unable to recover right fermions mass from fermionmatrix::D_KS_eo object.", __FILE__, __LINE__);
}

//Class MdagM_eo
void physics::fermionmatrix::MdagM_eo::operator()(const physics::lattices::Staggeredfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Staggeredfield_eo& in) const
{
	hmc_float mass = get_mass();
	if(upper_left==EVEN){
		//mass**2 - Deo*Doe
		DKS_eo(&tmp, gf, in, ODD);
		DKS_eo(out, gf, tmp, EVEN);
	} else {
		//mass**2 - Doe*Deo
		DKS_eo(&tmp, gf, in, EVEN);
		DKS_eo(out, gf, tmp, ODD);
	}
	saxpby(out, {mass*mass, 0.}, in, {-1., 0.}, *out);
}

cl_ulong physics::fermionmatrix::MdagM_eo::get_flops() const
{
	const hardware::System& system = get_system();
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	auto fermion_code = devices[0]->get_fermion_staggered_code();
	cl_ulong res;
	res = 2*fermion_code->get_flop_size("D_KS_eo");
	res += spinor_code->get_flop_size("saxpby_cplx_staggered_eoprec");
	
	return res;
}

bool physics::fermionmatrix::MdagM_eo::get_upper_left() const
{
	return upper_left;
}

hmc_float physics::fermionmatrix::MdagM_eo::get_mass() const noexcept
{
	return Fermionmatrix_stagg_basic::get_mass();
}

