/** @file
 * Implementation of the bicgstab algorithm
 *
 * Copyright (c) 2012-2014 Christopher Pinke <pinke@th.uni-frankfurt.de>
 * Copyright (c) 2012-2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#include "bicgstab.hpp"

/**
 * A "save" version of the bicgstab algorithm.
 * It is explicitely checked if the true residuum also fullfills the break condition.
 * This is chosen if "bicgstab_save" is selected.
 */
static int bicgstab_save(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& A, const physics::lattices::Gaugefield& gf,
                         const physics::lattices::Spinorfield& b, const hardware::System& system, physics::InterfacesHandler& interfacesHandler,
                         hmc_float prec, const physics::AdditionalParameters& additionalParameters);
static int bicgstab_save(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& A, const physics::lattices::Gaugefield& gf,
                         const physics::lattices::Spinorfield_eo& b, const hardware::System& system, physics::InterfacesHandler& interfacesHandler,
                         hmc_float prec, const physics::AdditionalParameters& additionalParameters);

/**
 * BICGstab implementation with a different structure than "save" one, similar to tmlqcd. This should be the default bicgstab.
 * In particular this version does not perform the check if the "real" residuum is sufficiently small!
 */
static int bicgstab_fast(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& A, const physics::lattices::Gaugefield& gf,
                         const physics::lattices::Spinorfield& b, const hardware::System& system, physics::InterfacesHandler& interfacesHandler,
                         hmc_float prec, const physics::AdditionalParameters& additionalParameters);
static int bicgstab_fast(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& A, const physics::lattices::Gaugefield& gf,
                         const physics::lattices::Spinorfield_eo& b, const hardware::System& system, physics::InterfacesHandler& interfacesHandler,
                         hmc_float prec, const physics::AdditionalParameters& additionalParameters);

static std::string create_log_prefix_solver(std::string name, int number) noexcept;
static std::string create_log_prefix_bicgstab(int number) noexcept;

int physics::algorithms::solvers::bicgstab(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& A,
                                           const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system,
                                           physics::InterfacesHandler& interfacesHandler, hmc_float prec, const physics::AdditionalParameters& additionalParameters)
{
    const physics::algorithms::SolversParametersInterface& parametersInterface = interfacesHandler.getSolversParametersInterface();

    if(parametersInterface.getSolver() == common::bicgstab_save) {
        return bicgstab_save(x, A, gf, b, system, interfacesHandler, prec, additionalParameters);
    } else {
        return bicgstab_fast(x, A, gf, b, system, interfacesHandler, prec, additionalParameters);
    }
}

static int bicgstab_save(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& f, const physics::lattices::Gaugefield& gf,
                         const physics::lattices::Spinorfield& b, const hardware::System& system, physics::InterfacesHandler& interfacesHandler,
                         hmc_float prec, const physics::AdditionalParameters& additionalParameters)
{
    // @todo this function often contains -1 in the comment but 1 in the code...
    using physics::lattices::Spinorfield;
    using physics::algorithms::solvers::SolverStuck;
    using physics::algorithms::solvers::SolverDidNotSolve;

    const physics::algorithms::SolversParametersInterface& parametersInterface = interfacesHandler.getSolversParametersInterface();

    /// @todo start timer synchronized with device(s)
    klepsydra::Monotonic timer;

    hmc_float resid;
    hmc_complex alpha, omega, rho { std::nan(""), std::nan("") };

    const Spinorfield v(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
    const Spinorfield p(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
    const Spinorfield rn(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
    const Spinorfield rhat(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
    const Spinorfield s(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
    const Spinorfield t(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
    const Spinorfield aux(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

    unsigned int iter = 0;
    log_squarenorm(create_log_prefix_bicgstab(iter) + "b: ", b);
    log_squarenorm(create_log_prefix_bicgstab(iter) + "x (initial): ", *x);

    for (iter = 0; iter < parametersInterface.getCgMax(); iter++) {
        if(iter % parametersInterface.getIterRefresh() == 0) {
            v.zero();
            p.zero();

            //initial r_n
            f(&rn, gf, *x, additionalParameters);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "rn: ", rn);
            saxpy(&rn, { 1., 0 }, rn, b);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "rn: ", rn);

            //rhat = r_n
            copyData(&rhat, rn);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "rhat: ", rhat);

            //set some constants to 1
            alpha = {1., 0.};
            omega = {1., 0.};
            rho = {1., 0.};
        }
        //rho_next = (rhat, rn)
        hmc_complex rho_next = scalar_product(rhat, rn);

        //check if algorithm is stuck
        //if rho is too small the algorithm will get stuck and will never converge!!
        if(std::abs(rho_next.re) < 1e-25 && std::abs(rho_next.im) < 1e-25 ) {
            //print the last residuum
            logger.fatal() << create_log_prefix_bicgstab(iter) << "Solver stuck at resid:\t" << resid;
            throw SolverStuck(iter, __FILE__, __LINE__);
        }

        //tmp1 = rho_next/rho = (rhat, rn)/..
        hmc_complex tmp1 = complexdivide(rho_next, rho);
        //rho_next = rho
        rho = rho_next;
        //tmp2 = alpha/omega = ...
        hmc_complex tmp2 = complexdivide(alpha, omega);
        //beta = tmp1*tmp2
        hmc_complex beta = complexmult(tmp1, tmp2);
        //tmp1 = beta*omega
        tmp1 = complexmult(beta, omega);
        //tmp2 = -tmp1
        tmp2 = complexsubtract( {0., 0.}, tmp1);
        //p = beta*p + tmp2*v + r_n = beta*p - beta*omega*v + r_n
        saxsbypz(&p, beta, p, tmp2, v, rn);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "p: ", p);

        //v = A*p
        f(&v, gf, p, additionalParameters);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "v: ", v);

        //tmp1 = (rhat, v)
        tmp1 = scalar_product(rhat, v);
        //alpha = rho/tmp1 = (..)/(rhat, v)
        alpha = complexdivide(rho, tmp1);
        //s = - alpha * v - r_n
        saxpy(&s, alpha, v, rn);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "s: ", s);

        //t = A s
        f(&t, gf, s, additionalParameters);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "t: ", t);

        //tmp1 = (t, s)
        tmp1 = scalar_product(t, s);
        //!!CP: this can also be global_squarenorm, but one needs a complex number here
        //tmp2 = (t,t)
        tmp2 = scalar_product(t, t);
        //omega = tmp1/tmp2 = (t,s)/(t,t)
        omega = complexdivide(tmp1, tmp2);

        //r_n = - omega*t - s
        saxpy(&rn, omega, t, s);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "rn: ", rn);

        //inout = alpha*p + omega * s + inout
        saxsbypz(x, alpha, p, omega, s, *x);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "x: ", *x);

        //resid = (rn,rn)
        resid = squarenorm(rn);
        logger.debug() << create_log_prefix_bicgstab(iter) << "resid: " << resid;

        //test if resid is NAN
        if(resid != resid) {
            logger.fatal() << create_log_prefix_bicgstab(iter) << "\tNAN occured!";
            throw SolverStuck(iter, __FILE__, __LINE__);
        }

        if(resid < prec) {
            //aux = A inout
            f(&aux, gf, *x, additionalParameters);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "aux: ", aux);

            //aux = -aux + source
            saxpy(&aux, {1., 0.}, aux, b);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "aux: ", aux);

            //trueresid = (aux, aux)
            hmc_float trueresid = squarenorm(aux);
            if(trueresid < prec) {
                logger.debug() << create_log_prefix_bicgstab(iter) << "Solver converged in " << iter << " iterations! true resid:\t" << trueresid;

                // report on performance
                if(logger.beInfo()) {
                    // we are always synchroneous here, as we had to recieve the residium from the device
                    const uint64_t duration = timer.getTime();

                    // calculate flops
                    /**
                     * @note this is not implemented since for non-eo this should not be of interest
                     */
                    // report performance
                    logger.info() << create_log_prefix_bicgstab(iter) << "Solver completed in " << duration / 1000 << " ms. Performed " << iter << " iterations";
                }

                log_squarenorm(create_log_prefix_bicgstab(iter) + "x (final): ", *x);
                return iter;
            }
        }
    }

    logger.fatal() << create_log_prefix_bicgstab(iter) << "Solver did not solve in " << parametersInterface.getCgMax() << " iterations. Last resid: " << resid;
    throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}

static int bicgstab_fast(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& f, const physics::lattices::Gaugefield& gf,
                         const physics::lattices::Spinorfield& b, const hardware::System& system, physics::InterfacesHandler& interfacesHandler,
                         hmc_float prec, const physics::AdditionalParameters& additionalParameters)
{
    using physics::lattices::Spinorfield;
    using physics::algorithms::solvers::SolverStuck;
    using physics::algorithms::solvers::SolverDidNotSolve;

    const physics::algorithms::SolversParametersInterface& parametersInterface = interfacesHandler.getSolversParametersInterface();

    /// @todo start timer synchronized with device(s)
    klepsydra::Monotonic timer;

    hmc_float resid;
    hmc_complex rho;

    const Spinorfield p(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
    const Spinorfield rn(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
    const Spinorfield rhat(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
    const Spinorfield v(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
    const Spinorfield s(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
    const Spinorfield t(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

    unsigned int iter = 0;
    log_squarenorm(create_log_prefix_bicgstab(iter) + "b: ", b);
    log_squarenorm(create_log_prefix_bicgstab(iter) + "x (initial): ", *x);

    for (iter = 0; iter < parametersInterface.getCgMax() ; iter++) {
        if(iter % parametersInterface.getIterRefresh() == 0) {
            //initial r_n, saved in p
            f(&rn, gf, *x, additionalParameters);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "rn: ", rn);
            saxpy(&p, { 1.0, 0 }, rn, b);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "p: ", p);

            //rhat = p
            copyData(&rhat, p);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "rhat: ", rhat);

            //r_n = p
            copyData(&rn, p);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "rn: ", rn);

            //rho = (rhat, rn)
            rho = scalar_product(rhat, rn);
        }

        //resid = (rn,rn)
        resid = squarenorm(rn);
        logger.debug() << create_log_prefix_bicgstab(iter) << "resid: " << resid;

        //test if resid is NAN
        if(resid != resid) {
            logger.fatal() << create_log_prefix_bicgstab(iter) << "NAN occured!";
            throw SolverStuck(iter, __FILE__, __LINE__);
        }
        if(resid < prec) {
            logger.debug() << create_log_prefix_bicgstab(iter) << "Solver converged in " << iter << " iterations! resid:\t" << resid;
            log_squarenorm(create_log_prefix_bicgstab(iter) + "x (final): ", *x);

            // report on performance
            if(logger.beInfo()) {
                // we are always synchroneous here, as we had to recieve the residium from the device
                const uint64_t duration = timer.getTime();

                // calculate flops
                /**
                 * @note this is not implemented since for non-eo this should not be of interest
                 */
                // report performance
                logger.info() << create_log_prefix_bicgstab(iter) << "Solver completed in " << duration / 1000 << " ms. Performed " << iter << " iterations";
            }

            return iter;
        }

        //v = A*p
        f(&v, gf, p, additionalParameters);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "v: ", v);

        //tmp1 = (rhat, v)
        hmc_complex tmp1 = scalar_product(rhat, v);
        //alpha = rho/tmp1 = (rhat, rn)/(rhat, v)
        hmc_complex alpha = complexdivide(rho, tmp1);

        //s = - alpha * v - r_n
        saxpy(&s, alpha, v, rn);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "s: ", s);

        //t = A s
        f(&t, gf, s, additionalParameters);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "t: ", t);

        //tmp1 = (t, s)
        tmp1 = scalar_product(t, s);
        //!!CP: this can also be global_squarenorm, but one needs a complex number here
        //tmp2 = (t,t)
        hmc_complex tmp2 = scalar_product(t, t);
        //omega = tmp1/tmp2 = (t,s)/(t,t)
        hmc_complex omega = complexdivide(tmp1, tmp2);

        //inout = alpha*p + omega * s + inout
        saxsbypz(x, alpha, p, omega, s, *x);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "x: ", *x);

        //r_n = - omega*t - s
        saxpy(&rn, omega, t, s);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "rn: ", rn);

        //rho_next = (rhat, rn)
        hmc_complex rho_next = scalar_product(rhat, rn);

        //check if algorithm is stuck
        //if rho is too small the algorithm will get stuck and will never converge!!
        if(std::abs(rho_next.re) < 1e-25 && std::abs(rho_next.im) < 1e-25) {
            //print the last residuum
            logger.fatal() << create_log_prefix_bicgstab(iter) << "Solver stuck at resid:" << resid;
            throw SolverStuck(iter, __FILE__, __LINE__);
        }

        //tmp1 = rho_next/rho = (rhat, rn)/..
        tmp1 = complexdivide(rho_next, rho);
        //tmp2 = alpha/omega = ...
        tmp2 = complexdivide(alpha, omega);
        //beta = tmp1*tmp2 = alpha*rho_next / (omega*rho)
        hmc_complex beta = complexmult(tmp1, tmp2);
        //tmp1 = beta*omega = alpha* rho_next / rho
        tmp1 = complexmult(beta, omega);
        //tmp2 = -tmp1
        tmp2 = complexsubtract( { 0., 0. }, tmp1);

        //p = beta*p + tmp2*v + r_n = beta*p - beta*omega*v + r_n
        saxsbypz(&p, beta, p, tmp2, v, rn);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "p: ", p);

        //rho_next = rho
        rho = rho_next;
    }

    logger.fatal() << create_log_prefix_bicgstab(iter) << "Solver did not solve in " << parametersInterface.getCgMax() << " iterations. Last resid: " << resid;
    throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}

int physics::algorithms::solvers::bicgstab(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& A,
                                           const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system,
                                           physics::InterfacesHandler& interfacesHandler, hmc_float prec, const physics::AdditionalParameters& additionalParameters)
{
    const physics::algorithms::SolversParametersInterface& parametersInterface = interfacesHandler.getSolversParametersInterface();

    if(parametersInterface.getSolver() == common::bicgstab_save) {
        return bicgstab_save(x, A, gf, b, system, interfacesHandler, prec, additionalParameters);
    } else {
        return bicgstab_fast(x, A, gf, b, system, interfacesHandler, prec, additionalParameters);
    }
}

static int bicgstab_save(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& f, const physics::lattices::Gaugefield& gf,
                         const physics::lattices::Spinorfield_eo& b, const hardware::System& system, physics::InterfacesHandler & interfacesHandler,
                         hmc_float prec, const physics::AdditionalParameters& additionalParameters)
{
    using physics::algorithms::solvers::SolverStuck;
    using physics::algorithms::solvers::SolverDidNotSolve;
    using namespace physics::lattices;

    // TODO start timer synchronized with device(s)
    klepsydra::Monotonic timer;

    const physics::algorithms::SolversParametersInterface& parametersInterface = interfacesHandler.getSolversParametersInterface();

    const Spinorfield_eo s(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
    const Spinorfield_eo t(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
    const Spinorfield_eo v(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
    const Spinorfield_eo p(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
    const Spinorfield_eo rn(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
    const Spinorfield_eo rhat(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
    const Spinorfield_eo aux(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

    unsigned retests = 0;
    int cgmax = parametersInterface.getCgMax();

    const Scalar<hmc_complex> alpha(system);
    const Scalar<hmc_complex> beta(system);
    const Scalar<hmc_complex> omega(system);
    const Scalar<hmc_complex> rho(system);
    const Scalar<hmc_complex> rho_next(system);
    const Scalar<hmc_complex> tmp1(system);
    const Scalar<hmc_complex> tmp2(system);
    const Scalar<hmc_complex> one(system);
    one.store(hmc_complex_one);
    const Scalar<hmc_complex> minus_one(system);
    minus_one.store(hmc_complex_minusone);
    hmc_float resid = 1.;
    int iter = 0;

    // report source and initial solution
    log_squarenorm(create_log_prefix_bicgstab(iter) + "b (initial): ", b);
    log_squarenorm(create_log_prefix_bicgstab(iter) + "x (initial): ", *x);

    // comments correspond to the bicgstab_fast version
    for (iter = 0; iter < cgmax; iter++) {
        if(iter % parametersInterface.getIterRefresh() == 0) {
            v.zero();
            p.zero();

            f(&rn, gf, *x, additionalParameters);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "rn: ", rn);
            saxpy(&rn, one, rn, b);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "rn: ", rn);

            copyData(&rhat, rn);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "rhat: ", rhat);

            copyData(&alpha, one);
            copyData(&omega, one);
            copyData(&rho, one);
        }
        scalar_product(&rho_next, rhat, rn);

        //check if algorithm is stuck
        const hmc_complex test = rho_next.get();
        if(std::abs(test.re) < 1e-25 && std::abs(test.im) < 1e-25) {
            //print the last residuum
            logger.fatal() << create_log_prefix_bicgstab(iter) << "Solver stuck at resid:\t" << resid;
            throw SolverStuck(iter, __FILE__, __LINE__);
        }

        divide(&tmp1, rho_next, rho);
        copyData(&rho, rho_next);
        divide(&tmp2, alpha, omega);
        multiply(&beta, tmp1, tmp2);

        multiply(&tmp1, beta, omega);
        multiply(&tmp2, minus_one, tmp1);
        saxsbypz(&p, beta, p, tmp2, v, rn);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "p: ", p);

        f(&v, gf, p, additionalParameters);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "v: ", v);

        scalar_product(&tmp1, rhat, v);
        divide(&alpha, rho, tmp1);

        saxpy(&s, alpha, v, rn);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "s: ", s);

        f(&t, gf, s, additionalParameters);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "t: ", t);

        scalar_product(&tmp1, t, s);
        //!!CP: can this also be global_squarenorm??
        scalar_product(&tmp2, t, t);
        divide(&omega, tmp1, tmp2);

        saxpy(&rn, omega, t, s);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "rn: ", rn);

        saxsbypz(x, alpha, p, omega, s, *x);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "x: ", *x);

        resid = squarenorm(rn);
        logger.debug() << create_log_prefix_bicgstab(iter) << "resid: " << resid;

        //test if resid is NAN
        if(resid != resid) {
            logger.fatal() << create_log_prefix_bicgstab(iter) << "NAN occured in bicgstab_eo!";
            throw SolverStuck(iter, __FILE__, __LINE__);
        }

        if(resid < prec) {
            ++retests;

            f(&aux, gf, *x, additionalParameters);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "aux: ", aux);
            saxpy(&aux, one, aux, b);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "aux: ", aux);

            hmc_float trueresid = squarenorm(aux);
            if(trueresid < prec) {
                logger.debug() << create_log_prefix_bicgstab(iter) << "Solver converged in " << iter << " iterations! true resid:\t" << trueresid;

                // report on performance
                if(logger.beInfo()) {
                    // we are always synchroneous here, as we had to recieve the residium from the device
                    const uint64_t duration = timer.getTime();

                    // calculate flops
                    const unsigned refreshs = iter / parametersInterface.getIterRefresh() + 1;
                    const size_t mf_flops = f.get_flops();

                    cl_ulong total_flops = 4 * get_flops<Spinorfield_eo, scalar_product>(system) + 4 * get_flops<hmc_complex, complexdivide>()
                            + 3 * get_flops<hmc_complex, complexmult>() + 2 * get_flops<Spinorfield_eo, saxsbypz>(system) + 2 * mf_flops
                            + 2 * get_flops<Spinorfield_eo, saxpy>(system) + get_flops<Spinorfield_eo, squarenorm>(system);
                    total_flops *= iter;

                    total_flops += refreshs * (mf_flops + get_flops<Spinorfield_eo, saxpy>(system));
                    total_flops += retests * (mf_flops + get_flops<Spinorfield_eo, saxpy>(system) + get_flops<Spinorfield_eo, squarenorm>(system));

                    // report performanc
                    logger.info() << create_log_prefix_bicgstab(iter) << "BiCGstab_save completed in " << duration / 1000 << " ms @ "
                            << (total_flops / duration / 1000.f) << " Gflops. Performed " << iter << " iterations";
                }

                log_squarenorm(create_log_prefix_bicgstab(iter) + "x (final): ", *x);
                return iter;
            }
        }
    }

    logger.fatal() << create_log_prefix_bicgstab(iter) << "Solver did not solve in " << parametersInterface.getCgMax() << " iterations. Last resid: " << resid;
    throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}

static int bicgstab_fast(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& f, const physics::lattices::Gaugefield& gf,
                         const physics::lattices::Spinorfield_eo& b, const hardware::System& system, physics::InterfacesHandler & interfacesHandler,
                         hmc_float prec, const physics::AdditionalParameters& additionalParameters)
{
    using namespace physics::lattices;
    using physics::algorithms::solvers::SolverStuck;
    using physics::algorithms::solvers::SolverDidNotSolve;

    // TODO start timer synchronized with device(s)
    klepsydra::Monotonic timer;

    const physics::algorithms::SolversParametersInterface& parametersInterface = interfacesHandler.getSolversParametersInterface();

    const Spinorfield_eo p(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
    const Spinorfield_eo rn(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
    const Spinorfield_eo rhat(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
    const Spinorfield_eo s(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
    const Spinorfield_eo t(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
    const Spinorfield_eo v(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

    const Scalar<hmc_complex> alpha(system);
    const Scalar<hmc_complex> beta(system);
    const Scalar<hmc_complex> omega(system);
    const Scalar<hmc_complex> rho(system);
    const Scalar<hmc_complex> rho_next(system);
    const Scalar<hmc_complex> tmp1(system);
    const Scalar<hmc_complex> tmp2(system);
    const Scalar<hmc_complex> one(system);
    one.store(hmc_complex_one);
    const Scalar<hmc_complex> minus_one(system);
    minus_one.store(hmc_complex_minusone);

    hmc_float resid;
    unsigned int iter = 0;

    // report source and initial solution
    log_squarenorm(create_log_prefix_bicgstab(iter) + "b (initial): ", b);
    log_squarenorm(create_log_prefix_bicgstab(iter) + "x (initial): ", *x);

    for (iter = 0; iter < parametersInterface.getCgMax(); iter++) {
        if(iter % parametersInterface.getIterRefresh() == 0) {
            //initial r_n, saved in p
            f(&rn, gf, *x, additionalParameters);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "rn: ", rn);
            saxpy(&p, one, rn, b);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "p: ", p);

            //rhat = p
            copyData(&rhat, p);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "rhat: ", rhat);

            //r_n = p
            copyData(&rn, p);
            log_squarenorm(create_log_prefix_bicgstab(iter) + "rn: ", rn);

            //rho = (rhat, rn)
            scalar_product(&rho, rhat, rn);
        }
        //resid = (rn,rn)
        resid = squarenorm(rn);
        logger.debug() << create_log_prefix_bicgstab(iter) << "resid: " << resid;

        //test if resid is NAN
        if(resid != resid) {
            logger.fatal() << create_log_prefix_bicgstab(iter) << "\tNAN occured!";
            throw SolverStuck(iter, __FILE__, __LINE__);
        }

        if(resid < prec) {
            logger.debug() << create_log_prefix_bicgstab(iter) << "Solver converged in " << iter << " iterations! resid:\t" << resid;

            // report on performance
            if(logger.beInfo()) {
                // we are always synchroneous here, as we had to recieve the residium from the device
                const uint64_t duration = timer.getTime();

                // calculate flops
                const unsigned refreshs = iter / parametersInterface.getIterRefresh() + 1;
                const cl_ulong mf_flops = f.get_flops();

                cl_ulong total_flops = get_flops<Spinorfield_eo, squarenorm>(system) + 2 * mf_flops + 4 * get_flops<Spinorfield_eo, scalar_product>(system)
                        + 4 * get_flops<hmc_complex, complexdivide>() + 2 * get_flops<Spinorfield_eo, saxpy>(system)
                        + 2 * get_flops<Spinorfield_eo, saxsbypz>(system) + 3 * get_flops<hmc_complex, complexmult>();
                total_flops *= iter;

                total_flops += refreshs * (mf_flops + get_flops<Spinorfield_eo, saxpy>(system) + get_flops<Spinorfield_eo, scalar_product>(system));

                // report performanc
                logger.info() << create_log_prefix_bicgstab(iter) << "Solver completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f)
                        << " Gflops. Performed " << iter << " iterations";
            }

            log_squarenorm(create_log_prefix_bicgstab(iter) + "x (final): ", *x);
            return iter;
        }
        //v = A*p
        f(&v, gf, p, additionalParameters);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "v: ", v);

        //tmp1 = (rhat, v)
        scalar_product(&tmp1, rhat, v);
        //alpha = rho/tmp1 = (rhat, rn)/(rhat, v)
        divide(&alpha, rho, tmp1);
        //s = - alpha * v - r_n
        saxpy(&s, alpha, v, rn);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "s: ", s);

        //t = A s
        f(&t, gf, s, additionalParameters);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "t: ", t);

        //tmp1 = (t, s)
        scalar_product(&tmp1, t, s);
        //!!CP: this can also be global_squarenorm, but one needs a complex number here
        //tmp2 = (t,t)
        scalar_product(&tmp2, t, t);
        //omega = tmp1/tmp2 = (t,s)/(t,t)
        divide(&omega, tmp1, tmp2);

        //inout = alpha*p + omega * s + inout
        saxsbypz(x, alpha, p, omega, s, *x);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "x: ", *x);

        //r_n = - omega*t - s
        saxpy(&rn, omega, t, s);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "rn: ", rn);

        //rho_next = (rhat, rn)
        scalar_product(&rho_next, rhat, rn);
        const hmc_complex test = rho_next.get();
        //check if algorithm is stuck
        if(std::abs(test.re) < 1e-25 && std::abs(test.im) < 1e-25) {
            //print the last residuum
            logger.fatal() << create_log_prefix_bicgstab(iter) << "Solver stuck at resid:\t" << resid;
            throw SolverStuck(iter, __FILE__, __LINE__);
        }

        //tmp1 = rho_next/rho = (rhat, rn)/..
        divide(&tmp1, rho_next, rho);
        //tmp2 = alpha/omega = ...
        divide(&tmp2, alpha, omega);
        //beta = tmp1*tmp2 = alpha*rho_next / (omega*rho)
        multiply(&beta, tmp1, tmp2);
        //tmp1 = beta*omega = alpha* rho_next / rho
        multiply(&tmp1, beta, omega);
        //tmp2 = -tmp1
        multiply(&tmp2, minus_one, tmp1);

        //p = beta*p + tmp2*v + r_n = beta*p - beta*omega*v + r_n
        saxsbypz(&p, beta, p, tmp2, v, rn);
        log_squarenorm(create_log_prefix_bicgstab(iter) + "p: ", p);

        //rho_next = rho
        copyData(&rho, rho_next);
    }

    logger.fatal() << create_log_prefix_bicgstab(iter) << "Solver did not solve in " << parametersInterface.getCgMax() << " iterations. Last resid: " << resid;
    throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}

static std::string create_log_prefix_solver(std::string name, int number) noexcept
{
    using namespace std;
    string separator_big = "\t";
    string separator_small = " ";
    string label = "SOLVER";

    stringstream strnumber;
    strnumber.fill('0');
    /// @todo this should be length(cgmax)
    strnumber.width(6);
    strnumber << right << number;
    stringstream outfilename;
    outfilename << separator_big << label << separator_small << "[" << name << "]" << separator_small << "[" << strnumber.str() << "]:" << separator_big;
    string outputfile = outfilename.str();
    return outputfile;
}

static std::string create_log_prefix_bicgstab(int number) noexcept
{
    return create_log_prefix_solver("BICGSTAB", number);
}
