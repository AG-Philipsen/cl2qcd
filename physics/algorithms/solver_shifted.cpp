/** @file
 * Implementation of the shifted solver algorithms
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.uni-frankfurt.de>
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

#include "solvers/solvers.hpp"
#include "solver_shifted.hpp"

#include "../../host_functionality/logger.hpp"
#include "../lattices/scalar_complex.hpp"
#include "../lattices/algebra_real.hpp"
#include "../lattices/staggeredfield_eo.hpp"
#include <sstream>
#include <vector>
#include <numeric>

int physics::algorithms::solvers::cg_m(std::vector<std::shared_ptr<physics::lattices::Staggeredfield_eo> > x,
                                       const physics::fermionmatrix::Fermionmatrix_stagg_eo& A,
                                       const physics::lattices::Gaugefield& gf, const std::vector<hmc_float> sigma,
                                       const physics::lattices::Staggeredfield_eo& b, const hardware::System& system,
                                       physics::InterfacesHandler& interfacesHandler, hmc_float prec, const physics::AdditionalParameters& additionalParameters)
{
    physics::algorithms::solvers::SolverShifted<physics::lattices::Staggeredfield_eo, physics::fermionmatrix::Fermionmatrix_stagg_eo>
    solverShifted(x, A, gf, sigma, b, system, interfacesHandler, prec, additionalParameters);
    x = solverShifted.solve();
    return solverShifted.getNumberOfIterationsDone();
}


template<typename FERMIONFIELD, typename FERMIONMATRIX>
physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::SolverShifted(const std::vector<std::shared_ptr<FERMIONFIELD> > xIn, const FERMIONMATRIX& AIn,
                                                                                        const physics::lattices::Gaugefield& gfIn, const std::vector<hmc_float> sigmaIn,
                                                                                        const FERMIONFIELD& bIn, const hardware::System& systemIn,
                                                                                        physics::InterfacesHandler& interfacesHandlerIn, hmc_float prec,
                                                                                        const physics::AdditionalParameters& additionalParametersIn)
     : x(xIn), A(AIn), gf(gfIn), sigma(sigmaIn), b(bIn), system(systemIn), solverPrecision(prec), additionalParameters(additionalParametersIn),
       parametersInterface(interfacesHandlerIn.getSolversParametersInterface()),
       hasSystemBeSolved(false), numberOfEquations(sigmaIn.size()), iterationNumber(0), residuumValue(NAN),
       r(FERMIONFIELD{system, interfacesHandlerIn.getInterface<FERMIONFIELD>()}),
       p(FERMIONFIELD{system, interfacesHandlerIn.getInterface<FERMIONFIELD>()}),
       ps(numberOfEquations),
       zeta_prev(physics::lattices::Vector<hmc_float>{numberOfEquations, system}),
       zeta(physics::lattices::Vector<hmc_float>{numberOfEquations, system}),
       zeta_foll(physics::lattices::Vector<hmc_float>{numberOfEquations, system}),
       alpha_vec(physics::lattices::Vector<hmc_float>{numberOfEquations, system}),
       beta_vec(physics::lattices::Vector<hmc_float>{numberOfEquations, system}),
       shift(physics::lattices::Vector<hmc_float>{numberOfEquations, system}),
       single_system_converged(numberOfEquations, false),
       single_system_iter(),
       resultSquarenorm(numberOfEquations, 0),
       alpha_scalar_prev(physics::lattices::Scalar<hmc_float>{system}),
       alpha_scalar(physics::lattices::Scalar<hmc_float>{system}),
       beta_scalar_prev(physics::lattices::Scalar<hmc_float>{system}),
       beta_scalar(physics::lattices::Scalar<hmc_float>{system}),
       zero(physics::lattices::Scalar<hmc_float>{system}),
       v(FERMIONFIELD(system, interfacesHandlerIn.getInterface<FERMIONFIELD>())),
       tmp1(physics::lattices::Scalar<hmc_float>{system}),
       tmp2(physics::lattices::Scalar<hmc_float>{system}),
       tmp3(physics::lattices::Scalar<hmc_float>{system}),
       single_eq_resid{nullptr},
       single_eq_resid_host{nullptr}
{
    for(auto& psElement: ps){
        psElement = std::make_shared<FERMIONFIELD>(system, interfacesHandlerIn.getInterface<FERMIONFIELD>());
    }
    if(sigma.size() != x.size())
        throw std::invalid_argument("Wrong size of multi-shifted inverter parameters!");
    if(parametersInterface.getUseMergeKernelsSpinor() == true) {
        single_eq_resid = std::unique_ptr<physics::lattices::Vector<hmc_float> > { new physics::lattices::Vector<hmc_float>(numberOfEquations, system) };
        single_eq_resid_host = std::unique_ptr<std::vector<hmc_float> > { new std::vector<hmc_float> };

    }
    USE_ASYNC_COPY = parametersInterface.getCgUseAsyncCopy();
    MINIMUM_ITERATIONS = parametersInterface.getCgMinimumIterationCount();
    if(USE_ASYNC_COPY)
        logger.warn() << "Asynchroneous copying in the CG-M is currently unimplemented!";
    if(MINIMUM_ITERATIONS)
        logger.warn() << "Minimum iterations set to " << MINIMUM_ITERATIONS << " -- should be used *only* for CGM benchmarking!";
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::setInitialConditions()
{
    zero.store(0.0);
    zeta_prev.store(std::vector<hmc_float>(numberOfEquations, 1.));   // zeta_prev[i] = 1
    zeta.store(std::vector<hmc_float>(numberOfEquations, 1.));        // zeta[i] = 1
    alpha_vec.store(std::vector<hmc_float>(numberOfEquations, 0.));   // alpha[i] = 0
    shift.store(sigma);
    for (int i = 0; i < numberOfEquations; i++) {
        x[i]->set_zero();    // x[i] = 0
        copyData(ps[i].get(), b);   // ps[i] = b
    }
    copyData(&r, b);                        // r = b
    copyData(&p, b);                        // p = b
    scalar_product_real_part(&tmp1, r, r);  // set tmp1 = (r, r) for the first iteration
    beta_scalar.store(1.0);                 // beta_scalar = 1, here I should set beta_scalar_prev
                                            // but in this way I can set beta_scalar_prev at the begin
                                            // of the loop over iterationNumber consistently with following iterations.
    alpha_scalar.store(0.0);                // alpha_scalar = 0. The same as beta_scalar above.
}

//<<<<<<< HEAD
template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::updateBetaScalar()
{
    //v=A.p and tmp1=(r,r) and tmp3=(p,v) ---> beta_scalar=(-1)*tmp1/tmp3
    copyData(&beta_scalar_prev, beta_scalar);   //before updating beta_scalar its value is saved
    A(&v, gf, p, &additionalParameters);
    log_squarenorm(createLogPrefix() + "v: ", v);
    scalar_product_real_part(&tmp3, p, v);
    divide(&beta_scalar, tmp1, tmp3);   //tmp1 is set from previous iteration
    subtract(&beta_scalar, zero, beta_scalar);
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::updateAuxiliaryFieldR()
{
    //r+=beta_scalar*A.p ---> r = r + beta_scalar*v
    saxpy(&r, beta_scalar, v, r);
    log_squarenorm(createLogPrefix() + "r: ", r);
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::updateAlphaScalar()
{
    //We store in tmp2 the quantity (r,r) that we use later. When we check the residuum, then it is already calculated.
    scalar_product_real_part(&tmp2, r, r);
    //tmp2=(rNew,rNew) and tmp1=(r,r) ---> alpha_scalar = tmp2/tmp1
    copyData(&alpha_scalar_prev, alpha_scalar);   //before updating alpha_scalar its value is saved
    divide(&alpha_scalar, tmp2, tmp1);
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::updateAuxiliaryFieldP()
{
    //p = r + alpha_scalar*p
    saxpy(&p, alpha_scalar, p, r);
    log_squarenorm(createLogPrefix() + "p: ", p);
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::updateVectorQuantities()
{
    physics::lattices::update_zeta_cgm(&zeta_foll, zeta, zeta_prev, beta_scalar_prev, beta_scalar, alpha_scalar_prev, shift, numberOfEquations);
    physics::lattices::update_beta_cgm(&beta_vec, beta_scalar, zeta_foll, zeta, numberOfEquations);
    physics::lattices::update_alpha_cgm(&alpha_vec, alpha_scalar, zeta_foll, beta_vec, zeta, beta_scalar, numberOfEquations);
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::updateSingleFieldOfSolution(unsigned int index)
{
    //x[k] = x[k] - beta[k]*ps[k] --->  remember that in beta_vec we store (- beta[k])
    saxpy(x[index].get(), beta_vec, index, *ps[index], *x[index]);
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::updateSingleFieldOfAuxiliaryFieldPs(unsigned int index)
{
    //ps[k] = zeta_iii[k]*r + alpha[k]*ps[k]
    saxpby(ps[index].get(), zeta_foll, index, r, alpha_vec, index, *ps[index]);
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::checkFiledsSquarenormsForPossibleNaN()
{
    if(logger.beDebug()){
        for(int k = 0; k < numberOfEquations; k++){
            if(std::isnan(squarenorm(*x[k]))){
                logger.fatal() << createLogPrefix() << "NAN occurred in x[" << k << "] squarenorm!";
                throw SolverStuck(iterationNumber, __FILE__, __LINE__);

            }
            if(std::isnan(squarenorm(*ps[k]))){
                logger.fatal() << createLogPrefix() << "NAN occurred in ps[" << k << "] squarenorm!";
                throw SolverStuck(iterationNumber, __FILE__, __LINE__);
            }
        }
        if(std::isnan(squarenorm(r))){
            logger.fatal() << createLogPrefix() << "NAN occurred in r squarenorm!";
            throw SolverStuck(iterationNumber, __FILE__, __LINE__);
        }
        if(std::isnan(squarenorm(p))){
            logger.fatal() << createLogPrefix() << "NAN occurred in p squarenorm!";
            throw SolverStuck(iterationNumber, __FILE__, __LINE__);
        }
        if(std::isnan(squarenorm(v))){
            logger.fatal() << createLogPrefix() << "NAN occurred in v squarenorm!";
            throw SolverStuck(iterationNumber, __FILE__, __LINE__);
        }
    }
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::updateQuantitiesForFollowingIteration()
{
    copyData(&zeta_prev, zeta);
    copyData(&zeta, zeta_foll);
    copyData(&tmp1, tmp2);
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::calculateResiduumSingleEquation(bool useMergedSpinorKernels, unsigned int index)
{

    if(useMergedSpinorKernels == true) {
        //Single equation Residuum = ||zeta_foll[k] * r||^2  ---> v = zeta_foll[k] * r for each k
        sax_vec_and_squarenorm(single_eq_resid.get(), zeta_foll, r);
        *single_eq_resid_host = single_eq_resid->get();
    } else {
        //Single equation Residuum = ||zeta_foll[k] * r||^2  ---> v = zeta_foll[k] * r for given k
        sax(&v, zeta_foll, index, r);
    }

}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::checkIfSingleEquationConverged(unsigned int index)
{
    if(iterationNumber % parametersInterface.getCgIterationBlockSize() == 0) {
        if((!parametersInterface.getUseMergeKernelsSpinor() && (squarenorm(v) < solverPrecision))
                || (parametersInterface.getUseMergeKernelsSpinor() && ((*single_eq_resid_host)[index] < solverPrecision))) {
            single_system_converged[index] = true;
            single_system_iter.push_back((uint) iterationNumber);
            logger.debug() << " ===> System number " << index << " converged after " << iterationNumber << " iterations! resid = " << tmp2.get();
        }
    }
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
bool physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::hasSingleSystemConverged(unsigned int index)
{
    return single_system_converged[index];
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
bool physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::hasSystemConverged()
{
    if(single_system_iter.size() != (uint) numberOfEquations)
        return false;
    residuumValue = tmp2.get();
    logger.debug() << createLogPrefix() << "resid: " << residuumValue;
    if(std::isnan(residuumValue)) {
        logger.fatal() << createLogPrefix() << "NAN occured!";
        throw SolverStuck(iterationNumber, __FILE__, __LINE__);
    }
    return (residuumValue < solverPrecision);
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::makePerformanceReport()
{
    logger.debug() << createLogPrefix() << "Solver converged in " << iterationNumber << " iterations! resid:\t" << residuumValue;
    if(logger.beInfo()) {
        // we are always synchroneous here, as we had to recieve the residuum from the device
        const uint64_t duration = timer.getTime();
        const uint64_t duration_noWarmup = timer_noWarmup.getTime();

        const cl_ulong matrix_flops = A.get_flops();
        const int sum_of_partial_iter = std::accumulate(single_system_iter.begin(), single_system_iter.end(), 0);
        logger.debug() << "matrix_flops: " << matrix_flops;
        cl_ulong flops_per_iter_no_inner_loop = matrix_flops
                + 2 * physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::scalar_product_real_part>(system) + 3
                + 2 * physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::saxpy>(system)
                + physics::lattices::get_flops_update_cgm("alpha", numberOfEquations, system)
                + physics::lattices::get_flops_update_cgm("beta", numberOfEquations, system)
                + physics::lattices::get_flops_update_cgm("zeta", numberOfEquations, system);
        cl_ulong flops_per_iter_only_inner_loop = physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::saxpy>(system)
                + physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::saxpby>(system)
                + physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::sax>(system)
                + physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, physics::lattices::squarenorm>(system);
        cl_ulong total_flops = iterationNumber * flops_per_iter_no_inner_loop + sum_of_partial_iter * flops_per_iter_only_inner_loop;
        cl_ulong noWarmup_flops = (iterationNumber - 1) * flops_per_iter_no_inner_loop + (sum_of_partial_iter - numberOfEquations) * flops_per_iter_only_inner_loop;

        logger.debug() << "total_flops: " << total_flops;
        logger.info() << createLogPrefix() << "CG-M completed in " << std::setprecision(6) << duration / 1000.f << " ms @ "
                << ((hmc_float) total_flops / duration / 1000.f) << " Gflops. Performed " << iterationNumber << " iterations. Performance after warmup: "
                << ((hmc_float) noWarmup_flops / duration_noWarmup / 1000.f) << " Gflops.";
    }
    debugLogSquarenormSetOfFields(createLogPrefix() + "x (final): ", x, numberOfEquations);
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::resetNoWarmupTimer()
{
    timer_noWarmup.reset();
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::debugMakeReportOfFieldsSquarenorm()
{
    if(logger.beDebug()){
        const int numberOfComponentsToBePrinted = numberOfEquations;
        if(numberOfComponentsToBePrinted > numberOfEquations)
            throw std::invalid_argument("In cg-m numberOfComponentsToBePrinted cannot be bigger than numberOfEquations!");
        if(iterationNumber == 0){
            log_squarenorm(createLogPrefix() + "b: ", b);
            log_squarenorm(createLogPrefix() + "r: ", r);
            log_squarenorm(createLogPrefix() + "p: ", p);
        }
        debugLogSquarenormSetOfFields(createLogPrefix() + "x", x, numberOfComponentsToBePrinted);
        debugLogSquarenormSetOfFields(createLogPrefix() + "ps", ps, numberOfComponentsToBePrinted);
    }
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::debugCalculateSquarenormOfResultField()
{
    if(logger.beDebug()) {
        for (uint i = 0; i < x.size(); i++)
            resultSquarenorm[i] = squarenorm(*x[i]);
    }
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
std::string physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::createLogPrefix()
{
    std::string separator_big = "\t";
    std::string separator_small = " ";
    std::string label = "SOLVER";
    std::string name = "CG-M";
    std::stringstream strnumber;
    strnumber.fill('0');
    strnumber.width(std::to_string(parametersInterface.getCgMax()).length());
    strnumber << std::right << iterationNumber;
    std::stringstream logPrefix;
    logPrefix << separator_big << label << separator_small << "[" << name << "]" << separator_small << "[" << strnumber.str() << "]:" << separator_big;
    return logPrefix.str();
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::debugLogSquarenormSetOfFields(const std::string& message,
                                                 const std::vector<std::shared_ptr<physics::lattices::Staggeredfield_eo> > setOfFields, const int reportNumber)
{
    if(logger.beDebug()) {
        for (int i = 0; i < reportNumber; i++) {
            std::ostringstream messageComplete(message);
            messageComplete <<  "[field_" << i << "]: ";
            physics::lattices::log_squarenorm(messageComplete.str(), *setOfFields[i]);
        }
    }
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
void physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::debugCompareSquarenormsOfResultFieldsBetweenBeginAndEndOfIteration()
{
    if(logger.beDebug()){
        hmc_float tmp;
        logger.debug() << "===============================================";
        for (uint i = 0; i < x.size(); i++) {
            tmp = squarenorm(*x[i]);
            logger.debug() << ((i < 10) ? " " : "") << "delta_sqnorm[field_" << i << "]: " << std::scientific << std::setprecision(16) << tmp - resultSquarenorm[i];
        }
        logger.debug() << "===============================================";
    }
}


template<typename FERMIONFIELD, typename FERMIONMATRIX>
unsigned int physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::getNumberOfIterationsDone()
{
    if(hasSystemBeSolved)
        return iterationNumber;
    else
        throw Print_Error_Message("SolverShifted not solved but asked for the number of iterations done!");
}

template<typename FERMIONFIELD, typename FERMIONMATRIX>
const std::vector<std::shared_ptr<FERMIONFIELD> > physics::algorithms::solvers::SolverShifted<FERMIONFIELD, FERMIONMATRIX>::solve()
{
    if(hasSystemBeSolved) return x;
    if(squarenorm(b) == 0) {
        logger.warn() << "CG-M solver called with zero field as r.h.s. -> trivial solution!";
        for (uint i = 0; i < sigma.size(); i++) x[i]->set_zero();
        hasSystemBeSolved = true;
        return x;
    }
    setInitialConditions();
    debugMakeReportOfFieldsSquarenorm();
    while(iterationNumber < parametersInterface.getCgMax() || iterationNumber < MINIMUM_ITERATIONS){
        updateBetaScalar();
        updateAuxiliaryFieldR();
        updateAlphaScalar();
        debugCalculateSquarenormOfResultField();
        updateAuxiliaryFieldP();
        updateVectorQuantities();
        calculateResiduumSingleEquation(parametersInterface.getUseMergeKernelsSpinor());
        for(uint indexEquation = 0; indexEquation < numberOfEquations; indexEquation++){
            if(hasSingleSystemConverged(indexEquation) == false){
                updateSingleFieldOfSolution(indexEquation);
                updateSingleFieldOfAuxiliaryFieldPs(indexEquation);
                calculateResiduumSingleEquation(parametersInterface.getUseMergeKernelsSpinor(), indexEquation);
                checkIfSingleEquationConverged(indexEquation);
            }
        }
        checkFiledsSquarenormsForPossibleNaN();
        updateQuantitiesForFollowingIteration();
        debugCompareSquarenormsOfResultFieldsBetweenBeginAndEndOfIteration();
        debugMakeReportOfFieldsSquarenorm();
        if(hasSystemConverged() && iterationNumber >= MINIMUM_ITERATIONS) {
            makePerformanceReport();
            hasSystemBeSolved = true;
            return x;
        }
        if(iterationNumber == 0) resetNoWarmupTimer();
        iterationNumber++;
    }
    logger.fatal() << createLogPrefix() << "Solver did not solve in " << parametersInterface.getCgMax() << " iterations. Last resid: " << residuumValue;
    throw SolverDidNotSolve(iterationNumber, __FILE__, __LINE__);
}















