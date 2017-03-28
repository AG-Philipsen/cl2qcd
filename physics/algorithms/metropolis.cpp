/** @file
 * Implementation of the metropolis algorithm
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2012-2014 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
 * Copyright (c) 2013, 2017 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
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

#include "metropolis.hpp"

#include "solvers/solvers.hpp"
#include "solver_shifted.hpp"
#include "../lattices/util.hpp"
#include "../../meta/util.hpp"
#include <cmath>
#include "../observables/gaugeObservables.hpp"

static void print_info_debug(physics::InterfacesHandler& interfacesHandler, std::string metropolis_part, hmc_float value, bool info = true);

hmc_float physics::algorithms::calc_s_fermion(const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system,
                                              physics::InterfacesHandler& interfacesHandler, const physics::AdditionalParameters& additionalParameters)
{
    using physics::lattices::Spinorfield;
    using namespace physics::algorithms::solvers;
    using namespace physics::fermionmatrix;

    const physics::algorithms::MetropolisParametersInterface & parametersInterface = interfacesHandler.getMetropolisParametersInterface();

    const Spinorfield phi_inv(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

    logger.debug() << "\t\t\tstart solver";

    /**
     * @todo at the moment, we can only put in a cold spinorfield
     * or a point-source spinorfield as trial-solution
     */
    const Spinorfield solution(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
    solution.cold();
    int iterations = 0;

    if(parametersInterface.getSolver() == common::cg) {
        const QplusQminus fm(system, interfacesHandler.getInterface<physics::fermionmatrix::QplusQminus>());
        iterations = cg(&solution, fm, gf, phi, system, interfacesHandler, parametersInterface.getSolverPrec(), additionalParameters);
        const Qminus qminus(system, interfacesHandler.getInterface<physics::fermionmatrix::Qminus>());
        qminus(&phi_inv, gf, solution, additionalParameters);

    } else {
        const Qplus fm(system, interfacesHandler.getInterface<physics::fermionmatrix::Qplus>());
        iterations = bicgstab(&solution, fm, gf, phi, system, interfacesHandler, parametersInterface.getSolverPrec(), additionalParameters);
        copyData(&phi_inv, solution);
    }
    logger.trace() << "Calulated S_fermion, solver took " << iterations << " iterations.";
    return squarenorm(phi_inv);
}

hmc_float physics::algorithms::calc_s_fermion(const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system,
                                              physics::InterfacesHandler& interfacesHandler, const physics::AdditionalParameters& additionalParameters)
{
    using physics::lattices::Spinorfield_eo;
    using namespace physics::algorithms::solvers;
    using namespace physics::fermionmatrix;

    const physics::algorithms::MetropolisParametersInterface & parametersInterface = interfacesHandler.getMetropolisParametersInterface();

    const Spinorfield_eo phi_inv(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

    logger.debug() << "\t\t\tstart solver";

    /** @todo at the moment, we can only put in a cold spinorfield
     * or a point-source spinorfield as trial-solution
     */
    const Spinorfield_eo solution(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
    int iterations = 0;
    logger.debug() << "\t\t\tstart solver";

    //the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
    if(parametersInterface.getSolver() == common::cg) {
        solution.cold();

        const QplusQminus_eo fm(system, interfacesHandler.getInterface<physics::fermionmatrix::QplusQminus_eo>());
        iterations = cg(&solution, fm, gf, phi, system, interfacesHandler, parametersInterface.getSolverPrec(), additionalParameters);
        const Qminus_eo qminus(system, interfacesHandler.getInterface<physics::fermionmatrix::Qminus_eo>());
        qminus(&phi_inv, gf, solution, additionalParameters);
    } else {
        solution.zero();
        solution.gamma5();
        const Qplus_eo fm(system, interfacesHandler.getInterface<physics::fermionmatrix::Qplus_eo>());
        iterations = bicgstab(&solution, fm, gf, phi, system, interfacesHandler, parametersInterface.getSolverPrec(), additionalParameters);
        copyData(&phi_inv, solution);
    }
    logger.trace() << "Calulated S_fermion, solver took " << iterations << " iterations.";
    return squarenorm(phi_inv);
}

hmc_float physics::algorithms::calc_s_fermion_mp(const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi,
                                                 const hardware::System& system, physics::InterfacesHandler& interfacesHandler)
{
    //this function essentially performs the same steps as in the non mass-prec case, however, one has to apply one more matrix multiplication
    //  therefore, comments are deleted here...
    //  Furthermore, in the bicgstab-case, the second inversions are not needed

    using physics::lattices::Spinorfield;
    using namespace physics::algorithms::solvers;
    using namespace physics::fermionmatrix;

    logger.trace() << "\tHMC [DH]:\tcalc final fermion energy...";

    const Spinorfield phi_inv(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

    const physics::algorithms::MetropolisParametersInterface & parametersInterface = interfacesHandler.getMetropolisParametersInterface();
    const physics::AdditionalParameters& additionalParameters = interfacesHandler.getAdditionalParameters<physics::lattices::Spinorfield>();
    const physics::AdditionalParameters& additionalParametersMp = interfacesHandler.getAdditionalParameters<physics::lattices::Spinorfield>(true);

    const Spinorfield tmp(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
    const Qplus qplus_mp(system, interfacesHandler.getInterface<physics::fermionmatrix::Qplus>());
    qplus_mp(&tmp, gf, phi, additionalParametersMp);

    const Spinorfield solution(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
    solution.cold();

    logger.debug() << "\t\t\tstart solver";
    int iterations = 0;

    if(parametersInterface.getSolver() == common::cg) {
        const QplusQminus fm(system, interfacesHandler.getInterface<physics::fermionmatrix::QplusQminus>());
        iterations = cg(&solution, fm, gf, tmp, system, interfacesHandler, parametersInterface.getSolverPrec(), additionalParameters);
        const Qminus qminus(system, interfacesHandler.getInterface<physics::fermionmatrix::Qminus>());
        qminus(&phi_inv, gf, solution, additionalParameters);
    } else {
        const Qplus fm(system, interfacesHandler.getInterface<physics::fermionmatrix::Qplus>());
        iterations = bicgstab(&solution, fm, gf, tmp, system, interfacesHandler, parametersInterface.getSolverPrec(), additionalParameters);
        copyData(&phi_inv, solution);
    }
    logger.trace() << "Calulated S_fermion, solver took " << iterations << " iterations.";
    return squarenorm(phi_inv);
}

hmc_float physics::algorithms::calc_s_fermion_mp(const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi,
                                                 const hardware::System& system, physics::InterfacesHandler& interfacesHandler)
{
    //this function essentially performs the same steps as in the non mass-prec case, however, one has to apply one more matrix multiplication
    //  therefore, comments are deleted here...
    //  Furthermore, in the bicgstab-case, the second inversions are not needed
    using physics::lattices::Spinorfield_eo;
    using namespace physics::algorithms::solvers;
    using namespace physics::fermionmatrix;

    logger.trace() << "\tHMC [DH]:\tcalc final fermion energy...";

    const physics::algorithms::MetropolisParametersInterface & parametersInterface = interfacesHandler.getMetropolisParametersInterface();
    const physics::AdditionalParameters& additionalParameters = interfacesHandler.getAdditionalParameters<physics::lattices::Spinorfield_eo>();
    const physics::AdditionalParameters& additionalParametersMp = interfacesHandler.getAdditionalParameters<physics::lattices::Spinorfield_eo>(true);

    const Spinorfield_eo phi_inv(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

    logger.debug() << "\t\t\tstart solver";
    int iterations = 0;

    //sf_tmp = Qplus(light_mass) phi_mp
    const Spinorfield_eo tmp(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
    const Qplus_eo qplus_mp(system, interfacesHandler.getInterface<physics::fermionmatrix::Qplus_eo>());
    qplus_mp(&tmp, gf, phi, additionalParametersMp);

    const Spinorfield_eo solution(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

    logger.debug() << "\t\t\tstart solver";

    if(parametersInterface.getSolver() == common::cg) {
        solution.cold();
        const QplusQminus_eo fm(system, interfacesHandler.getInterface<physics::fermionmatrix::QplusQminus_eo>());
        iterations = cg(&solution, fm, gf, tmp, system, interfacesHandler, parametersInterface.getSolverPrec(), additionalParameters);
        const Qminus_eo qminus(system, interfacesHandler.getInterface<physics::fermionmatrix::Qminus_eo>());
        qminus(&phi_inv, gf, solution, additionalParameters);
    } else {
        solution.zero();
        solution.gamma5();

        const Qplus_eo fm(system, interfacesHandler.getInterface<physics::fermionmatrix::Qplus_eo>());
        iterations = bicgstab(&solution, fm, gf, tmp, system, interfacesHandler, parametersInterface.getSolverPrec(), additionalParameters);

        copyData(&phi_inv, solution);
    }
    logger.trace() << "Calulated S_fermion, solver took " << iterations << " iterations.";
    return squarenorm(phi_inv);
}

hmc_float physics::algorithms::calc_s_fermion_mp(const physics::lattices::Gaugefield&, const physics::lattices::Rooted_Staggeredfield_eo&,
                                                 const hardware::System&, physics::InterfacesHandler&)
{
    throw std::runtime_error("Not implemented!");
}

/**
 * This function returns the value of the fermionic part of the action for the RHMC, i.e.
 * \f$ \phi^*\,(M^\dag\,M)^{-\frac{N_f}{4}}\,\phi \f$.
 * 
 * Here, we use the coefficients of the rational expansion of \f$ x^{-\frac{N_f}{4}} \f$ 
 * included in the Rooted_Staggeredfield_eo object to calculate \f$ (M^\dag\,M)^{-\frac{N_f}{4}}\,\phi \f$
 * (using the multi-shifted inverter). Then a scalar product give the returning value.
 * 
 */
hmc_float physics::algorithms::calc_s_fermion(const physics::lattices::Gaugefield& gf, const physics::lattices::Rooted_Staggeredfield_eo& phi,
                                              const hardware::System& system, physics::InterfacesHandler& interfacesHandler,
                                              const physics::AdditionalParameters& additionalParameters)
{
    using physics::lattices::Rooted_Staggeredfield_eo;
    using namespace physics::algorithms::solvers;
    using namespace physics::fermionmatrix;

    logger.trace() << "\tRHMC [DH]:\tcalc final fermion energy...";

    const physics::algorithms::MetropolisParametersInterface & parametersInterface = interfacesHandler.getMetropolisParametersInterface();

    const physics::fermionmatrix::MdagM_eo fm(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>());
    int iterations = 0;

    //Temporary fields for shifted inverter
    logger.debug() << "\t\tstart solver...";
    std::vector<std::shared_ptr<physics::lattices::Staggeredfield_eo> > X;
    for (int i = 0; i < phi.getOrder(); i++)
        X.emplace_back(std::make_shared<physics::lattices::Staggeredfield_eo>(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>()));
    //Here the inversion must be performed with high precision, because it'll be used for Metropolis test
    iterations = physics::algorithms::solvers::cg_m(X, fm, gf, phi.get_b(), phi[0], system, interfacesHandler, parametersInterface.getSolverPrec(), additionalParameters);
    logger.debug() << "\t\t...end solver in " << iterations << " iterations";

    //this is to reconstruct (MdagM)^{-\frac{N_f}{4}}\,\phi
    physics::lattices::Staggeredfield_eo tmp(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
    sax(&tmp, { phi.get_a0(), 0. }, phi[0]);
    for (int i = 0; i < phi.getOrder(); i++) {
        saxpy(&tmp, { (phi.get_a())[i], 0. }, *X[i], tmp);
    }

    logger.trace() << "Calulated S_fermion, solver took " << iterations << " iterations.";
    return scalar_product(phi[0], tmp).re;

}

template<class SPINORFIELD> static hmc_observables metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf,
        const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p, const physics::lattices::Gaugemomenta& new_p,
        const SPINORFIELD& phi, const hmc_float spinor_energy_init, const SPINORFIELD * const phi_mp, const hmc_float spinor_energy_mp_init,
        const hardware::System& system, physics::InterfacesHandler& interfacesHandler)
{
    using namespace physics::algorithms;

    const physics::algorithms::MetropolisParametersInterface& parametersInterface = interfacesHandler.getMetropolisParametersInterface();
    const physics::observables::GaugeObservablesParametersInterface& gaugeobservablesParameters = interfacesHandler.getGaugeObservablesParametersInterface();

    //Calc Hamiltonian
    print_info_debug(interfacesHandler, "[DH]:\tCalculate Hamiltonian", sqrt(-1.), false);
    hmc_float deltaH = 0.;
    hmc_float s_old = 0.;
    hmc_float s_new = 0.;

    //Gauge-Part
    hmc_float plaq = physics::observables::measurePlaquetteWithoutNormalization(&gf, gaugeobservablesParameters);
    hmc_float plaq_new = physics::observables::measurePlaquetteWithoutNormalization(&new_u, gaugeobservablesParameters);

    hmc_float rect_new = 0.;
    hmc_float rect = 0.;

    if(parametersInterface.getUseRectangles() == true) {
        rect = physics::observables::measureRectangles(&gf, gaugeobservablesParameters);
        rect_new = physics::observables::measureRectangles(&new_u, gaugeobservablesParameters);
        hmc_float c0 = parametersInterface.getC0();
        hmc_float c1 = parametersInterface.getC1();
        deltaH = -beta * (c0 * (plaq - plaq_new) + c1 * (rect - rect_new));
        s_old = -beta * (c0 * (plaq) + c1 * (rect));
        s_new = -beta * (c0 * (plaq_new) + c1 * (rect_new));
    } else {
        /** NOTE: the minus here is introduced to fit tmlqcd!!! */
        deltaH = -(plaq - plaq_new) * beta;
        s_old = -(plaq) * beta;
        s_new = -(plaq_new) * beta;
    }

    print_info_debug(interfacesHandler, "[DH]:\tS[GF]_0:\t", s_old, false);
    print_info_debug(interfacesHandler, "[DH]:\tS[GF]_1:\t", s_new, false);
    print_info_debug(interfacesHandler, "[DH]:\tdS[GF]: \t", deltaH);
    //check on NANs
    if(s_old != s_old || s_new != s_new || deltaH != deltaH) {
        throw Print_Error_Message("NAN occured in Metropolis! Aborting!", __FILE__, __LINE__);
    }

    //Gaugemomentum-Part
    hmc_float p2 = squarenorm(p);
    hmc_float new_p2 = squarenorm(new_p);
    if(parametersInterface.getFermact() != common::action::rooted_stagg) {
        //the energy is half the squarenorm
        deltaH += 0.5 * (p2 - new_p2);
    } else {
        //The gaugemomenta part of the action in the staggered code is:
        // \sum_{n,\mu} 0.5 * Tr[P_\mu(n)P_\mu(n)] = 0.25 \sum_{n,\mu,A} ||P^A_\mu(n)||^2
        p2 *= 0.5;
        new_p2 *= 0.5;   //I multiply here by 0.5 so that in the following cout we have the correct numbers
        deltaH += 0.5 * (p2 - new_p2);
    }

    print_info_debug(interfacesHandler, "[DH]:\tS[GM]_0:\t", 0.5 * p2, false);
    print_info_debug(interfacesHandler, "[DH]:\tS[GM]_1:\t", 0.5 * new_p2, false);
    print_info_debug(interfacesHandler, "[DH]:\tdS[GM]: \t", 0.5 * (p2 - new_p2));
    //check on NANs
    if(p2 != p2 || new_p2 != new_p2 || deltaH != deltaH) {
        throw Print_Error_Message("NAN occured in Metropolis! Aborting!", __FILE__, __LINE__);
    }

    //Fermion-Part:
    if(!parametersInterface.getUseGaugeOnly()) {
        if(parametersInterface.getUseMp()) {
            if(parametersInterface.getFermact() == common::action::rooted_stagg) {
                throw Invalid_Parameters("Mass preconditioning not implemented for staggered fermions!", "NOT rooted_stagg", parametersInterface.getFermact());
            }
            //in this case one has contributions from det(m_light/m_heavy) and det(m_heavy)
            // det(m_heavy)
            hmc_float s_fermion_final;
            //initial energy has been computed in the beginning...
            s_fermion_final = calc_s_fermion(new_u, phi, system, interfacesHandler, interfacesHandler.getAdditionalParameters<SPINORFIELD>(true));
            deltaH += spinor_energy_init - s_fermion_final;

            print_info_debug(interfacesHandler, "[DH]:\tS[DET]_0:\t", spinor_energy_init, false);
            print_info_debug(interfacesHandler, "[DH]:\tS[DET]_1:\t", s_fermion_final, false);
            print_info_debug(interfacesHandler, "[DH]:\tdS[DET]:\t", spinor_energy_init - s_fermion_final);
            //check on NANs
            if(spinor_energy_init != spinor_energy_init || s_fermion_final != s_fermion_final || deltaH != deltaH) {
                throw Print_Error_Message("NAN occured in Metropolis! Aborting!", __FILE__, __LINE__);
            }

            // det(m_light/m_heavy)
            //initial energy has been computed in the beginning...
            hmc_float s_fermion_mp_final = calc_s_fermion_mp(new_u, *phi_mp, system, interfacesHandler);
            deltaH += spinor_energy_mp_init - s_fermion_mp_final;

            print_info_debug(interfacesHandler, "[DH]:\tS[DETRAT]_0:\t", spinor_energy_mp_init, false);
            print_info_debug(interfacesHandler, "[DH]:\tS[DETRAT]_1:\t", s_fermion_mp_final, false);
            print_info_debug(interfacesHandler, "[DH]:\tdS[DETRAT]:\t", spinor_energy_mp_init - s_fermion_mp_final);
            //check on NANs
            if(spinor_energy_mp_init != spinor_energy_mp_init || s_fermion_mp_final != s_fermion_mp_final || deltaH != deltaH) {
                throw Print_Error_Message("NAN occured in Metropolis! Aborting!", __FILE__, __LINE__);
            }
        } else {
            hmc_float s_fermion_final = calc_s_fermion(new_u, phi, system, interfacesHandler, interfacesHandler.getAdditionalParameters<SPINORFIELD>(false));
            deltaH += spinor_energy_init - s_fermion_final;

            print_info_debug(interfacesHandler, "[DH]:\tS[DET]_0:\t", spinor_energy_init, false);
            print_info_debug(interfacesHandler, "[DH]:\tS[DET]_1:\t", s_fermion_final, false);
            print_info_debug(interfacesHandler, "[DH]:\tdS[DET]: \t", spinor_energy_init - s_fermion_final);
            //check on NANs
            if(spinor_energy_init != spinor_energy_init || s_fermion_final != s_fermion_final || deltaH != deltaH) {
                throw Print_Error_Message("NAN occured in Metropolis! Aborting!", __FILE__, __LINE__);
            }
        }
    }
    //Metropolis-Part
    hmc_float compare_prob;
    if(deltaH < 0) {
        compare_prob = std::exp(deltaH);
    } else {
        compare_prob = 1.0;
    }
    print_info_debug(interfacesHandler, "[DH]:\tdH:\t\t", deltaH);
    print_info_debug(interfacesHandler, "[MET]:\tAcc-Prop:\t", compare_prob);

    //calc gaugeobservables
    //todo: calc only of final configuration
    auto plaqs = physics::observables::measureAllPlaquettes(&gf, gaugeobservablesParameters);
    plaq = plaqs.plaquette;
    hmc_float splaq = plaqs.spatialPlaquette;
    hmc_float tplaq = plaqs.temporalPlaquette;

    plaqs = physics::observables::measureAllPlaquettes(&new_u, gaugeobservablesParameters);
    plaq_new = plaqs.plaquette;
    hmc_float splaq_new = plaqs.spatialPlaquette;
    hmc_float tplaq_new = plaqs.temporalPlaquette;

    hmc_complex poly = physics::observables::measurePolyakovloop(&gf, gaugeobservablesParameters);
    hmc_complex poly_new = physics::observables::measurePolyakovloop(&new_u, gaugeobservablesParameters);

    hmc_observables tmp;
    if(rnd <= compare_prob) {
        tmp.accept = 1;
        tmp.plaq = plaq_new;
        tmp.tplaq = tplaq_new;
        tmp.splaq = splaq_new;
        tmp.poly = poly_new;
        tmp.deltaH = deltaH;
        tmp.prob = compare_prob;
        if(parametersInterface.getUseRectangles())
            tmp.rectangles = rect_new / parametersInterface.getRectanglesNormalization();
    } else {
        tmp.accept = 0;
        tmp.plaq = plaq;
        tmp.tplaq = tplaq;
        tmp.splaq = splaq;
        tmp.poly = poly;
        tmp.deltaH = deltaH;
        tmp.prob = compare_prob;
        if(parametersInterface.getUseRectangles())
            tmp.rectangles = rect / parametersInterface.getRectanglesNormalization();
    }

    return tmp;
}

hmc_observables physics::algorithms::metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf,
        const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p, const physics::lattices::Gaugemomenta& new_p,
        const physics::lattices::Spinorfield& phi, const hmc_float spinor_energy_init, const physics::lattices::Spinorfield * const phi_mp,
        const hmc_float spinor_energy_mp_init, const hardware::System& system, physics::InterfacesHandler& interfacesHandler)
{
    return ::metropolis(rnd, beta, gf, new_u, p, new_p, phi, spinor_energy_init, phi_mp, spinor_energy_mp_init, system, interfacesHandler);
}

hmc_observables physics::algorithms::metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf,
        const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p, const physics::lattices::Gaugemomenta& new_p,
        const physics::lattices::Spinorfield_eo& phi, const hmc_float spinor_energy_init, const physics::lattices::Spinorfield_eo * const phi_mp,
        const hmc_float spinor_energy_mp_init, const hardware::System& system, physics::InterfacesHandler& interfacesHandler)
{
    return ::metropolis(rnd, beta, gf, new_u, p, new_p, phi, spinor_energy_init, phi_mp, spinor_energy_mp_init, system, interfacesHandler);
}

hmc_observables physics::algorithms::metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf,
        const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p, const physics::lattices::Gaugemomenta& new_p,
        const physics::lattices::Rooted_Staggeredfield_eo& phi, const hmc_float spinor_energy_init,
        const physics::lattices::Rooted_Staggeredfield_eo* const phi_mp, const hmc_float spinor_energy_mp_init, const hardware::System& system,
        physics::InterfacesHandler& interfacesHandler)
{
    return ::metropolis(rnd, beta, gf, new_u, p, new_p, phi, spinor_energy_init, phi_mp, spinor_energy_mp_init, system, interfacesHandler);
}

static void print_info_debug(physics::InterfacesHandler& interfacesHandler, std::string metropolis_part, hmc_float value, bool info)
{

    if(info == false) {
        if(logger.beDebug()) {
            if(interfacesHandler.getMetropolisParametersInterface().getFermact() != common::action::rooted_stagg) {
                if(value == value)
                    logger.debug() << "\tHMC " << metropolis_part << std::setprecision(10) << value;
                else
                    logger.debug() << "\tHMC " << metropolis_part;
            } else {
                if(value == value)
                    logger.debug() << "\tRHMC " << metropolis_part << std::setprecision(10) << value;
                else
                    logger.debug() << "\tRHMC " << metropolis_part;
            }
        }
    } else {
        if(interfacesHandler.getMetropolisParametersInterface().getFermact() != common::action::rooted_stagg) {
            if(value == value)
                logger.info() << "\tHMC " << metropolis_part << std::setprecision(10) << value;
            else
                logger.info() << "\tHMC " << metropolis_part;
        } else {
            if(value == value)
                logger.info() << "\tRHMC " << metropolis_part << std::setprecision(10) << value;
            else
                logger.info() << "\tRHMC " << metropolis_part;
        }
    }

}
