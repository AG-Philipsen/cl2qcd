/** @file
 * Implementation of the fermion_force functions
 *
 * Copyright (c) 2013, 2017 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
 * Copyright (c) 2017 Francesca Cuteri <cuteri@th.physik.uni-frankfurt.de>
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

#include "fermion_force_staggered.hpp"

#include "solver_shifted.hpp"
#include "../../hardware/code/molecular_dynamics.hpp"

/**
 * This function reconstructs the fermionic contribution to the force (in the RHMC). Now, here
 * it is particularly easy to get lost because of minus signs. What is called force is somehow
 * arbitrary. For sure, instead, there is no doubts about the definition of the time derivative
 * of the momentum field conjugated to the gaugefield, it is -dS/dq. Now, if we refer to the
 * Gattringer (page 197) notation (as done in the Wilson code), we choose to call force F_\mu(n)
 * only dS/dq and later, in the update of the gaugemomentum, we will take into account the minus sign.
 * Indicating the gaugemomentum field as H_\mu(n), we have that F_\mu(n) = - Hdot_\mu(n). To be
 * coherent with the Wilson code, we have to add here a minus sign to obtain F_\mu(n) from Hdot_\mu(n).
 * Starting from the field Hdot, we have
 * @code
 * Hdot_\mu(n) = -i * [U_\mu(n)*\sum_{i=1}^k c_i Q^i_\mu(n)]_TA
 * @endcode
 * where k is the order of rational approximation, c_i are the numerators and
 * @code
 *               | +eta_\mu(n) (D_oe X_e^i)_{n+\mu} (X_e^i^\dag)_n     if evenodd = EVEN 
 *  Q^i_\mu(n) = | 
 *               | -eta_\mu(n) (X_e^i)_{n+\mu} ((D_oe X_e^i)^\dag)_n   if evenodd = ODD 
 * 
 * @endcode
 * If we put U_\mu(n) into Q^i_\mu(n) we have
 * @code
 * Hdot_\mu(n) = -i * [\sum_{i=1}^k c_i QQ^i_\mu(n)]_TA
 *             = -i * \sum_{i=1}^k [c_i QQ^i_\mu(n)]_TA
 *             = \sum_{i=1}^k {c_i * (-i) * [QQ^i_\mu(n)]_TA}
 * @endcode
 * with
 * @code
 *                | +eta_\mu(n) U_\mu(n) * (D_oe X_e^i)_{n+\mu} (X_e^i^\dag)_n     if evenodd = EVEN 
 *  QQ^i_\mu(n) = | 
 *                | -eta_\mu(n) U_\mu(n) * (X_e^i)_{n+\mu} ((D_oe X_e^i)^\dag)_n   if evenodd = ODD 
 * 
 * @endcode
 * Now, (-i) * [QQ^i_\mu(n)]_TA is exactly what the function fermion_force calculates, given the fields
 * (D_oe X_e^i) and (X_e^i). So, basically, here we have to calculate them with the
 * multi-shifted inverter:
 * @code
 *  X^i_e = (M^\dagM + p_i)^{-1} * phi_e  ==>  D_oe X_e^i
 * @endcode
 * and then reconstruct the force, using the Rational Coefficients:
 * @code
 * Hdot_\mu(n) = \sum_{i=1}^k {c_i * out_fermion_force} ==> F_\mu(n) = \sum_{i=1}^k {-c_i * out_fermion_force}
 * @endcode
 * where we add a minus sign to pass from Hdot_\mu(n) to F_\mu(n).
 * 
 * @note To perform the sum above, the saxpy operation of the gaugemomenta is used.
 * 
 * @warning Remember that this function add to the Gaugemomenta field "force" the fermionic
 *          contribution. Therefore such a field must be properly initialized.
 * 
 * @attention If an imaginary chemical potential is used, this function is not modified,
 *            because chem_pot_im is included in the kernel. See force_staggered_fermion_eo.cl
 *            file documentation for further information.
 */
void physics::algorithms::calc_fermion_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
        const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfacesHandler,
        const physics::AdditionalParameters& additionalParameters)
{
    using physics::lattices::Staggeredfield_eo;
    using namespace physics::algorithms::solvers;
    using namespace physics::fermionmatrix;

    const physics::algorithms::ForcesParametersInterface & parametersInterface = interfacesHandler.getForcesParametersInterface();
    logger.debug() << "\t\tcalc_fermion_force...";

    logger.debug() << "\t\t\tstart solver";
    std::vector<std::shared_ptr<Staggeredfield_eo> > X;
    std::vector<std::shared_ptr<Staggeredfield_eo> > Y;
    for(unsigned int i = 0; i < phi.getOrder(); i++) {
        X.emplace_back(std::make_shared<Staggeredfield_eo>(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>()));
        Y.emplace_back(std::make_shared<Staggeredfield_eo>(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>()));
    }
    const MdagM_eo fm(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>());
    cg_m(X, fm, gf, phi.get_b(), *phi[0], system, interfacesHandler, parametersInterface.getForcePreconditioning(), additionalParameters);
    logger.debug() << "\t\t\t  end solver";

    //Now that I have X^i I can calculate Y^i = D_oe X_e^i and in the same for loop
    //reconstruct the force. I will use a temporary Gaugemomenta to calculate the
    //partial force (on the whole lattice) that will be later added to "force"
    const D_KS_eo Doe(system, interfacesHandler.getInterface<physics::fermionmatrix::D_KS_eo>(), ODD);   //with ODD it is the Doe operator
    physics::lattices::Gaugemomenta tmp(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

    for(unsigned int i = 0; i < phi.getOrder(); i++) {
        Doe(Y[i].get(), gf, *X[i]);
        tmp.zero();
        fermion_force(&tmp, *Y[i], *X[i], gf, EVEN);
        fermion_force(&tmp, *X[i], *Y[i], gf, ODD);
        physics::lattices::saxpy(force, -1. * (phi.get_a())[i], tmp);
    }

    logger.debug() << "\t\t...end calc_fermion_force!";
}

void physics::algorithms::calc_fermion_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
        const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler,
        const physics::AdditionalParameters& additionalParameters)
{
    using physics::lattices::Gaugefield;
    using namespace physics::algorithms;

    calc_fermion_force(force, gf, phi, system, interfaceHandler, additionalParameters);

    const physics::algorithms::ForcesParametersInterface & parametersInterface = interfaceHandler.getForcesParametersInterface();
    if(parametersInterface.getUseSmearing() == true) {
        throw Print_Error_Message("Smeared Gaugefield force is not implemented.", __FILE__, __LINE__);
    }
}

void physics::algorithms::fermion_force(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Staggeredfield_eo& A,
        const physics::lattices::Staggeredfield_eo& B, const physics::lattices::Gaugefield& gf, const int evenodd)
{
    auto gm_bufs = gm->get_buffers();
    auto A_bufs = A.get_buffers();
    auto B_bufs = B.get_buffers();
    auto gf_bufs = gf.get_buffers();
    size_t num_bufs = gm_bufs.size();
    if(num_bufs != A_bufs.size() || num_bufs != B_bufs.size() || num_bufs != gf_bufs.size()) {
        throw Print_Error_Message(std::string(__func__) + " is only implemented for a single device.", __FILE__, __LINE__);
    }

    for (size_t i = 0; i < num_bufs; ++i) {
        auto gm_buf = gm_bufs[i];
        auto A_buf = A_bufs[i];
        auto B_buf = B_bufs[i];
        auto gf_buf = gf_bufs[i];
        auto code = gm_buf->get_device()->getMolecularDynamicsCode();
        code->fermion_staggered_partial_force_device(gf_buf, A_buf, B_buf, gm_buf, evenodd);
    }
    gm->update_halo();
}

