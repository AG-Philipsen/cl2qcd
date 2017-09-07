/** @file
 * Implementation of the molecular dynamics algorithms
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2012-2013 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
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

#include "molecular_dynamics.hpp"

#include <stdexcept>
#include <type_traits>
#include "../fermionmatrix/fermionmatrix.hpp"
#include "../../meta/util.hpp"
#include "solvers/solvers.hpp"
#include "solver_shifted.hpp"
#include "forces.hpp"
#include "../../hardware/code/molecular_dynamics.hpp"
#include "../lattices/util.hpp"

void physics::algorithms::md_update_gaugefield(const physics::lattices::Gaugefield * const gf, const physics::lattices::Gaugemomenta& gm, const hmc_float eps)
{
    auto gf_bufs = gf->get_buffers();
    auto gm_bufs = gm.get_buffers();
    size_t num_bufs = gf_bufs.size();
    if(num_bufs != gm_bufs.size()) {
        throw std::invalid_argument("The gaugefield and the gaugemomenta use different devices.");
    }

    for (size_t i = 0; i < num_bufs; ++i) {
        auto gf_buf = gf_bufs[i];
        auto gm_buf = gm_bufs[i];
        auto code = gf_buf->get_device()->getMolecularDynamicsCode();
        code->md_update_gaugefield_device(gm_buf, gf_buf, eps);
    }

    logger.debug() << "\tHMC [UP]:\tupdate GF [" << eps << "]";
}

void physics::algorithms::md_update_spinorfield(const physics::lattices::Spinorfield * const out, const physics::lattices::Gaugefield& gf,
                                                const physics::lattices::Spinorfield& orig, const hardware::System& system,
                                                physics::InterfacesHandler & interfacesHandler, const physics::AdditionalParameters& additionalParameters)
{
    logger.debug() << "\tHMC [UP]:\tupdate SF";
    physics::fermionmatrix::Qplus qplus(system, interfacesHandler.getInterface<physics::fermionmatrix::Qplus>());
    qplus(out, gf, orig, additionalParameters);
    log_squarenorm("Spinorfield after update", *out);
}

void physics::algorithms::md_update_spinorfield(const physics::lattices::Spinorfield_eo * const out, const physics::lattices::Gaugefield& gf,
                                                const physics::lattices::Spinorfield_eo& orig, const hardware::System& system,
                                                physics::InterfacesHandler & interfacesHandler, const physics::AdditionalParameters& additionalParameters)
{
    logger.debug() << "\tHMC [UP]:\tupdate SF";
    physics::fermionmatrix::Qplus_eo qplus(system, interfacesHandler.getInterface<physics::fermionmatrix::Qplus_eo>());
    qplus(out, gf, orig, additionalParameters);
    log_squarenorm("Spinorfield after update", *out);
}

/**
 * This is the update of the Staggeredfield_eo in the RHMC. Actually, it is a tool to obtain
 * the new Staggeredfield_eo at the beginning of each iteration, given a gaussian drawn field.
 * In fact, this function makes the operator (Mdag*M) to some rational power act on the
 * Staggeredfield_eo "orig". To be clearer, here the following operation is implemented:
 * @code
 *   out = (a_0 + \sum_{i=1}^k a_i * (Mdag*M + b_i)^{-1} ) * orig 
 *       = a_0 * orig + \sum_{i=1}^k a_i * [(Mdag*M + b_i)^{-1} * orig]
 * @endcode
 * where k is the order of the rational approximation, a_0, a_i and b_i are the coefficients.
 * 
 * @note The coefficients of the approximation are stored in out, since this is the field that
 *       appears in the perform_RHMC_step function.
 */
void physics::algorithms::md_update_spinorfield(const physics::lattices::Rooted_Staggeredfield_eo * out, const physics::lattices::Gaugefield& gf,
                                                const physics::lattices::Rooted_Staggeredfield_eo& orig, const hardware::System& system,
                                                physics::InterfacesHandler & interfacesHandler, const physics::AdditionalParameters& additionalParameters)
{
    logger.debug() << "\tRHMC [UP]:\tupdate SF";
    const physics::algorithms::MolecularDynamicsInterface & parametersInterface = interfacesHandler.getMolecularDynamicsInterface();
    const physics::fermionmatrix::MdagM_eo fm(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>());

    const unsigned int numberOfPseudofermions = interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>().getNumberOfPseudofermions();
    for(unsigned int j=0; j<numberOfPseudofermions; j++){
        logger.trace() << "\t\tstart solver...";
        //Temporary fields for shifted inverter
        std::vector<std::shared_ptr<physics::lattices::Staggeredfield_eo> > X;
        for (unsigned int i = 0; i < out->getOrder(); i++)
            X.emplace_back(std::make_shared<physics::lattices::Staggeredfield_eo>(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>()));
        //Here the inversion must be performed with high precision, because it'll be used for Metropolis test
        const int iterations = physics::algorithms::solvers::cg_m(X, fm, gf, out->get_b(), *orig[j], system, interfacesHandler, parametersInterface.getSolverPrec(), additionalParameters);
        logger.trace() << "\t\t...end solver in " << iterations << " iterations";

        const physics::lattices::Staggeredfield_eo& pseudofermionInOut = *(*out)[j];
        physics::lattices::sax(&pseudofermionInOut, { out->get_a0(), 0. }, *orig[j]);
        for (unsigned int i = 0; i < out->getOrder(); i++)
            physics::lattices::saxpy(&pseudofermionInOut, { (out->get_a())[i], 0. }, *X[i], pseudofermionInOut);

        log_squarenorm("Staggeredfield_eo after update", pseudofermionInOut);
    }
}

/**
 * template for md_update_spinorfield_mp
 * this needs 3 fermionmatrices in order to use cg as default solver (because for the cg one needs a hermitian matrix)
 */
template<class FERMIONMATRIX, class FERMIONMATRIX_CONJ, class FERMIONMATRIX_HERM, class SPINORFIELD> static void md_update_spinorfield_mp(
        const SPINORFIELD * const out, const physics::lattices::Gaugefield& gf, const SPINORFIELD& orig, const hardware::System& system,
        physics::InterfacesHandler& interfacesHandler)
{
    SPINORFIELD temporarySpinorfield(system, interfacesHandler.getInterface<SPINORFIELD>());
    FERMIONMATRIX qplus(system, interfacesHandler.getInterface<FERMIONMATRIX>());
    const physics::algorithms::MolecularDynamicsInterface & parametersInterface = interfacesHandler.getMolecularDynamicsInterface();
    const physics::AdditionalParameters& additionalParameters = interfacesHandler.getAdditionalParameters<SPINORFIELD>(false);
    const physics::AdditionalParameters& additionalParametersMp = interfacesHandler.getAdditionalParameters<SPINORFIELD>(true);

    log_squarenorm("Spinorfield before update: ", orig);

    qplus(&temporarySpinorfield, gf, orig, additionalParameters);

    /**
     * Now one needs ( Qplus )^-1 (heavy_mass) using tmp as source to get phi_mp
     */
    try {
        /**
         * Use BiCGStab as default here
         * an exception will be thrown if the solver cannot solve
         */
        logger.debug() << "\t\t\tstart solver";

        /**
         * @todo at the moment, we can only put in a cold spinorfield
         * or a point-source spinorfield as trial-solution
         */
        out->zero();
        out->gamma5();

        FERMIONMATRIX qplusMp(system, interfacesHandler.getInterface<FERMIONMATRIX>());
        physics::algorithms::solvers::bicgstab(out, qplusMp, gf, temporarySpinorfield, system, interfacesHandler, parametersInterface.getSolverPrec(), additionalParametersMp);
    }   //try
    catch (physics::algorithms::solvers::SolverException& e) {
        logger.fatal() << e.what();
        logger.info() << "Retry with CG...";
        SPINORFIELD tmp2(system, interfacesHandler.getInterface<SPINORFIELD>());
        /**
         * @todo at the moment, we can only put in a cold spinorfield
         * or a point-source spinorfield as trial-solution
         */
        tmp2.zero();
        tmp2.gamma5();

        FERMIONMATRIX_HERM fm_herm(system, interfacesHandler.getInterface<FERMIONMATRIX_HERM>());
        physics::algorithms::solvers::cg(&tmp2, fm_herm, gf, temporarySpinorfield, system, interfacesHandler, parametersInterface.getSolverPrec(), additionalParametersMp);
        FERMIONMATRIX_CONJ fmConj(system, interfacesHandler.getInterface<FERMIONMATRIX_CONJ>());
        fmConj(out, gf, tmp2, additionalParametersMp);
    }
}

void physics::algorithms::md_update_spinorfield_mp(const physics::lattices::Spinorfield * const out, const physics::lattices::Gaugefield& gf,
                                                   const physics::lattices::Spinorfield& orig, const hardware::System& system,
                                                   physics::InterfacesHandler& interfacesHandler)
{
    logger.debug() << "\tHMC [UP]:\tupdate SF_MP";
    using physics::fermionmatrix::Qplus;
    using physics::fermionmatrix::Qminus;
    using physics::fermionmatrix::QplusQminus;

    ::md_update_spinorfield_mp<Qplus, Qminus, QplusQminus>(out, gf, orig, system, interfacesHandler);
}

void physics::algorithms::md_update_spinorfield_mp(const physics::lattices::Spinorfield_eo * const out, const physics::lattices::Gaugefield& gf,
                                                   const physics::lattices::Spinorfield_eo& orig, const hardware::System& system,
                                                   physics::InterfacesHandler& interfacesHandler)
{
    logger.debug() << "\tHMC [UP]:\tupdate SF_MP";
    using physics::fermionmatrix::Qplus_eo;
    using physics::fermionmatrix::Qminus_eo;
    using physics::fermionmatrix::QplusQminus_eo;

    ::md_update_spinorfield_mp<Qplus_eo, Qminus_eo, QplusQminus_eo>(out, gf, orig, system, interfacesHandler);
}

template<class SPINORFIELD> static void md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps,
                                                                const physics::lattices::Gaugefield& gf, const SPINORFIELD& phi, const hardware::System& system,
                                                                physics::InterfacesHandler& interfacesHandler, const physics::AdditionalParameters& additionalParameters)
{
    using namespace physics::algorithms;

    physics::lattices::Gaugemomenta delta_p(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
    delta_p.zero();
    calc_total_force(&delta_p, gf, phi, system, interfacesHandler, additionalParameters);

    logger.debug() << "\t[R]HMC " << "[UP]:\tupdate GM [" << eps << "]";
    physics::lattices::saxpy(inout, -1. * eps, delta_p);
}
void physics::algorithms::md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf,
                                                  const physics::lattices::Spinorfield& phi, const hardware::System& system,
                                                  physics::InterfacesHandler& interfaceHandler, const physics::AdditionalParameters& additionalParameters)
{
    ::md_update_gaugemomentum(inout, eps, gf, phi, system, interfaceHandler, additionalParameters);
}
void physics::algorithms::md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf,
                                                  const physics::lattices::Spinorfield_eo& phi, const hardware::System& system,
                                                  physics::InterfacesHandler& interfaceHandler, const physics::AdditionalParameters& additionalParameters)
{
    ::md_update_gaugemomentum(inout, eps, gf, phi, system, interfaceHandler, additionalParameters);
}

void physics::algorithms::md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf,
                                                  const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system,
                                                  physics::InterfacesHandler& interfaceHandler, const physics::AdditionalParameters& additionalParameters)
{
    ::md_update_gaugemomentum(inout, eps, gf, phi, system, interfaceHandler, additionalParameters);
}

//With staggered fermions we do not need the last argument mubar and than we do not use the template above
//void physics::algorithms::md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf,
//                                                  const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system,
//                                                  physics::InterfacesHandler& interfaceHandler, const physics::AdditionalParameters& additionalParameters)
//{
//    using namespace physics::algorithms;
//
//    physics::lattices::Gaugemomenta delta_p(system, interfaceHandler.getInterface<physics::lattices::Gaugemomenta>());
//    delta_p.zero();
//    calc_total_force(&delta_p, gf, phi, system, interfaceHandler, mass);
//
//    logger.debug() << "\tRHMC " << "[UP]:\tupdate GM [" << eps << "]";
//    physics::lattices::saxpy(inout, -1. * eps, delta_p);
//}

void physics::algorithms::md_update_gaugemomentum_gauge(const physics::lattices::Gaugemomenta * const gm, const hmc_float eps,
                                                        const physics::lattices::Gaugefield& gf, const hardware::System& system,
                                                        physics::InterfacesHandler& interfaceHandler)
{
    const physics::lattices::Gaugemomenta force(system, interfaceHandler.getInterface<physics::lattices::Gaugemomenta>());
    force.zero();
    calc_gauge_force(&force, gf, interfaceHandler);
    log_squarenorm("\tHMC [UP]:\tFORCE [GAUGE]:\t", force);

    logger.debug() << "\tHMC [UP]:\tupdate GM [" << eps << "]";
    physics::lattices::saxpy(gm, -1. * eps, force);
}

template<class SPINORFIELD> static void md_update_gaugemomentum_fermion(const physics::lattices::Gaugemomenta * const inout, hmc_float eps,
                                                                        const physics::lattices::Gaugefield& gf, const SPINORFIELD& phi,
                                                                        const hardware::System& system, physics::InterfacesHandler& interfacesHandler,
                                                                        const physics::AdditionalParameters& additionalParameters)
{
    using namespace physics::algorithms;

    const physics::lattices::Gaugemomenta force(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
    force.zero();
    calc_fermion_force(&force, gf, phi, system, interfacesHandler, additionalParameters);
    log_squarenorm("\t[R]HMC [UP]:\tFORCE [DET]:\t", force);
    logger.debug() << "\t[R]HMC " << "[UP]:\tupdate GM [" << eps << "]";
    physics::lattices::saxpy(inout, -1. * eps, force);
}
void physics::algorithms::md_update_gaugemomentum_fermion(const physics::lattices::Gaugemomenta * const inout, hmc_float eps,
                                                          const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi,
                                                          const hardware::System& system, physics::InterfacesHandler& interfaceHandler,
                                                          const physics::AdditionalParameters& additionalParameters)
{
    ::md_update_gaugemomentum_fermion(inout, eps, gf, phi, system, interfaceHandler, additionalParameters);
}
void physics::algorithms::md_update_gaugemomentum_fermion(const physics::lattices::Gaugemomenta * const inout, hmc_float eps,
                                                          const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi,
                                                          const hardware::System& system, physics::InterfacesHandler& interfaceHandler,
                                                          const physics::AdditionalParameters& additionalParameters)
{
    ::md_update_gaugemomentum_fermion(inout, eps, gf, phi, system, interfaceHandler, additionalParameters);
}

void physics::algorithms::md_update_gaugemomentum_fermion(const physics::lattices::Gaugemomenta * const inout, hmc_float eps,
                                                          const physics::lattices::Gaugefield& gf, const physics::lattices::Rooted_Staggeredfield_eo& phi,
                                                          const hardware::System& system, physics::InterfacesHandler& interfaceHandler,
                                                          const physics::AdditionalParameters& additionalParameters)
{
    ::md_update_gaugemomentum_fermion(inout, eps, gf, phi, system, interfaceHandler, additionalParameters);
}

//With staggered fermions we do not need the last argument mubar and than we do not use the template above
//void physics::algorithms::md_update_gaugemomentum_fermion(const physics::lattices::Gaugemomenta * const inout, hmc_float eps,
//                                                          const physics::lattices::Gaugefield& gf, const physics::lattices::Rooted_Staggeredfield_eo& phi,
//                                                          const hardware::System& system, physics::InterfacesHandler& interfaceHandler,
//                                                          const physics::AdditionalParameters& additionalParameters)
//{
//    using namespace physics::algorithms;
//
//    const physics::lattices::Gaugemomenta force(system, interfaceHandler.getInterface<physics::lattices::Gaugemomenta>());
//    force.zero();
//    calc_fermion_force(&force, gf, phi, system, interfaceHandler, mass);
//    log_squarenorm("\tRHMC [UP]:\tFORCE [DET]:\t", force);
//    logger.debug() << "\tRHMC " << "[UP]:\tupdate GM [" << eps << "]";
//    physics::lattices::saxpy(inout, -1. * eps, force);
//}

template<class SPINORFIELD> static void md_update_gaugemomentum_detratio(const physics::lattices::Gaugemomenta * const inout, hmc_float eps,
                                                                         const physics::lattices::Gaugefield& gf, const SPINORFIELD& phi_mp,
                                                                         const hardware::System& system, physics::InterfacesHandler& interfacesHandler)
{
    using namespace physics::algorithms;

    const physics::lattices::Gaugemomenta force(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
    force.zero();
    calc_detratio_forces(&force, gf, phi_mp, system, interfacesHandler);
    log_squarenorm("\tHMC [UP]:\tFORCE [DETRAT]:\t", force);

    logger.debug() << "\tHMC [UP]:\tupdate GM [" << eps << "]";
    physics::lattices::saxpy(inout, -1. * eps, force);
}
void physics::algorithms::md_update_gaugemomentum_detratio(const physics::lattices::Gaugemomenta * const inout, hmc_float eps,
                                                           const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi,
                                                           const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    ::md_update_gaugemomentum_detratio(inout, eps, gf, phi, system, interfaceHandler);
}
void physics::algorithms::md_update_gaugemomentum_detratio(const physics::lattices::Gaugemomenta * const inout, hmc_float eps,
                                                           const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi,
                                                           const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    ::md_update_gaugemomentum_detratio(inout, eps, gf, phi, system, interfaceHandler);
}
