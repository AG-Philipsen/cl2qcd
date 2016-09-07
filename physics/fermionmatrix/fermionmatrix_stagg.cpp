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
bool physics::fermionmatrix::Fermionmatrix_stagg_basic::isHermitian() const noexcept
{
	return isMatrixHermitian;
}

bool physics::fermionmatrix::Fermionmatrix_stagg_basic::hasMinimumEigenvalueThreshold() const noexcept
{
    return hasMatrixMinimumEigenvalueThreshold;
}

const hardware::System& physics::fermionmatrix::Fermionmatrix_stagg_basic::get_system() const noexcept
{
	return system;
}

hmc_float physics::fermionmatrix::Fermionmatrix_stagg_basic::getThresholdForMinimumEigenvalue(hmc_float) const
{
    throw Print_Error_Message("Threshold for minimum eigenvalue not existing or not implemented!");
}

//Class D_KS_eo
void physics::fermionmatrix::D_KS_eo::operator()(const physics::lattices::Staggeredfield_eo * out, const physics::lattices::Gaugefield& gf,
                                                 const physics::lattices::Staggeredfield_eo& in, const physics::AdditionalParameters* additionalParameters) const
{
    if(additionalParameters == NULL)
        DKS_eo(out, gf, in, evenodd);
    else
        throw Print_Error_Message("D_KS_eo operator applied passing to it some additional parameters! This should not happen!");
}

cl_ulong physics::fermionmatrix::D_KS_eo::get_flops() const
{
	const hardware::System& system = get_system();
	auto devices = system.get_devices();
	auto fermion_code = devices[0]->getFermionStaggeredCode();
	
	return fermion_code->get_flop_size("D_KS_eo");
}

//Class MdagM_eo
void physics::fermionmatrix::MdagM_eo::operator()(const physics::lattices::Staggeredfield_eo * out, const physics::lattices::Gaugefield& gf,
                                                  const physics::lattices::Staggeredfield_eo& in, const physics::AdditionalParameters* additionalParameters) const
{
    if(additionalParameters != NULL){
        if(upper_left==EVEN){
            //mass**2 - Deo*Doe
            DKS_eo(&tmp, gf, in, ODD);
            DKS_eo(out, gf, tmp, EVEN);
        } else {
            //mass**2 - Doe*Deo
            DKS_eo(&tmp, gf, in, EVEN);
            DKS_eo(out, gf, tmp, ODD);
        }
        hmc_float mass = additionalParameters->getMass();
        saxpby(out, {mass*mass, 0.}, in, {-1., 0.}, *out);
    } else
        throw Print_Error_Message("MdagM_eo operator applied without the MANDATORY additional parameters! Aborting...");
}

cl_ulong physics::fermionmatrix::MdagM_eo::get_flops() const
{
	const hardware::System& system = get_system();
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->getSpinorStaggeredCode();
	auto fermion_code = devices[0]->getFermionStaggeredCode();
	cl_ulong res;
	res = 2*fermion_code->get_flop_size("D_KS_eo");
	res += spinor_code->get_flop_size("saxpby_cplx_staggered_eoprec");
	
	return res;
}

bool physics::fermionmatrix::MdagM_eo::get_upper_left() const
{
	return upper_left;
}

hmc_float physics::fermionmatrix::MdagM_eo::getThresholdForMinimumEigenvalue(hmc_float mass) const
{
    return mass * mass;
}


