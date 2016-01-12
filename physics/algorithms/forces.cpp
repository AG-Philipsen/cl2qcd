/** @file
 * Implementation of the force algorithms
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2012-2013 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
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

#include "forces.hpp"

#include "../../meta/util.hpp"
#include "../../hardware/device.hpp"
#include "../../hardware/code/molecular_dynamics.hpp"

void physics::algorithms::gauge_force(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield& gf)
{
    auto gm_bufs = gm->get_buffers();
    auto gf_bufs = gf.get_buffers();
    size_t num_bufs = gm_bufs.size();
    if(num_bufs != gf_bufs.size()) {
        throw Print_Error_Message("Input buffers seem to use different devices.", __FILE__, __LINE__);
    }

    for (size_t i = 0; i < num_bufs; ++i) {
        auto gm_buf = gm_bufs[i];
        auto gf_buf = gf_bufs[i];
        auto code = gm_buf->get_device()->getMolecularDynamicsCode();
        code->gauge_force_device(gf_buf, gm_buf);
    }
    gm->update_halo();
}

void physics::algorithms::gauge_force_tlsym(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield& gf)
{
    auto gm_bufs = gm->get_buffers();
    auto gf_bufs = gf.get_buffers();
    size_t num_bufs = gm_bufs.size();
    if(num_bufs != gf_bufs.size()) {
        throw Print_Error_Message("Input buffers seem to use different devices.", __FILE__, __LINE__);
    }

    for (size_t i = 0; i < num_bufs; ++i) {
        auto gm_buf = gm_bufs[i];
        auto gf_buf = gf_bufs[i];
        auto code = gm_buf->get_device()->getMolecularDynamicsCode();
        code->gauge_force_tlsym_device(gf_buf, gm_buf);
    }
    gm->update_halo();
}

void physics::algorithms::calc_gauge_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf, physics::InterfacesHandler& interfacesHandler)
{
    gauge_force(gm, gf);
    if(interfacesHandler.getForcesParametersInterface().getUseRectangles()) {
        gauge_force_tlsym(gm, gf);
    }
}

template<class SPINORFIELD> static void calc_total_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
        const SPINORFIELD& phi, const hardware::System& system, physics::InterfacesHandler& interfacesHandler, const physics::AdditionalParameters& additionalParameters)
{
    using namespace physics::algorithms;

    force->zero();
    if(!interfacesHandler.getForcesParametersInterface().getUseGaugeOnly())
        calc_fermion_forces(force, gf, phi, system, interfacesHandler, additionalParameters);
    calc_gauge_force(force, gf, interfacesHandler);
}

void physics::algorithms::calc_total_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
        const physics::lattices::Spinorfield& phi, const hardware::System& system, physics::InterfacesHandler& interfacesHandler,
        const physics::AdditionalParameters& additionalParameters)
{
    ::calc_total_force(force, gf, phi, system, interfacesHandler, additionalParameters);
}
void physics::algorithms::calc_total_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
        const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfacesHandler,
        const physics::AdditionalParameters& additionalParameters)
{
    ::calc_total_force(force, gf, phi, system, interfacesHandler, additionalParameters);
}

void physics::algorithms::calc_total_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
        const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler,
        const physics::AdditionalParameters& additionalParameters)
{
    ::calc_total_force(force, gf, phi, system, interfaceHandler, additionalParameters);
}

//Here we do not need the last argument mubar and than we do not use the template above
//void physics::algorithms::calc_total_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
//        const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler,
//        const hmc_float mass)
//{
//    using namespace physics::algorithms;
//    force->zero();
//    if(!interfaceHandler.getForcesParametersInterface().getUseGaugeOnly())
//        calc_fermion_forces(force, gf, phi, system, interfaceHandler, mass);
//    calc_gauge_force(force, gf, interfaceHandler);
//}

