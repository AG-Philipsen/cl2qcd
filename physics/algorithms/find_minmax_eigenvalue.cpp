/** @file
 * Implementation of the algorithm to find min and max eigenvalues of an operator
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

#include "find_minmax_eigenvalue.hpp"
#include "../lattices/staggeredfield_eo.hpp"
#include "../lattices/util.hpp"
#include "../lattices/scalar.hpp"
#include "../../host_functionality/logger.hpp"
#include "cmath"
#include "sstream"
#include "solvers/exceptions.hpp"

static std::string create_log_prefix_find_max(int number) noexcept;
static std::string create_log_prefix_find_min(int number) noexcept;
static hmc_float find_min_knowing_max(const hmc_float max, const physics::fermionmatrix::Fermionmatrix_stagg_eo& A, const physics::lattices::Gaugefield& gf,
                                      const hardware::System& system,  physics::InterfacesHandler& interfacesHandler,
                                      hmc_float prec, const physics::AdditionalParameters& additionalParameters);

hmc_float physics::algorithms::find_max_eigenvalue(const physics::fermionmatrix::Fermionmatrix_stagg_eo& A, const physics::lattices::Gaugefield& gf,
                                                   const hardware::System& system, physics::InterfacesHandler& interfacesHandler, hmc_float prec,
                                                   const physics::AdditionalParameters& additionalParameters)
{
    using namespace physics::lattices;
    using namespace physics::algorithms;

    if(!(A.isHermitian()))
        throw std::invalid_argument("Unable to deal with non-hermitian matrices in find_max_eigenvalue!");

    Scalar<hmc_complex> max(system);
    hmc_float resid;

    //This timer is to know how long this function takes
    klepsydra::Monotonic timer;

    //This field is the starting point and it must be random (we have to be sure
    //to have a non zero component along the eigenvectors referring to the biggest eigenvalue)
    Staggeredfield_eo v1(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
    pseudo_randomize<Staggeredfield_eo, su3vec>(&v1, 123);
    sax(&v1, { 1. / sqrt(squarenorm(v1)), 0. }, v1);   //v1 is now normalized
    //Auxiliary field
    Staggeredfield_eo v2(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
    //How often to check resid
    const physics::algorithms::MinMaxEigenvalueParametersInterface & parametersInterface = interfacesHandler.getMinMaxEigenvalueParametersInterface();
    const int RESID_CHECK_FREQUENCY = parametersInterface.getFindMinMaxIterationBlockSize();

    log_squarenorm(create_log_prefix_find_max(0) + "v1 (initial): ", v1);
    log_squarenorm(create_log_prefix_find_max(0) + "v2 (initial) [not-initialized]: ", v2);

    for (unsigned int i = 0; i < parametersInterface.getFindMinMaxMaxValue(); i++) {
        //Apply A onto v1
        A(&v2, gf, v1, &additionalParameters);
        if(i % 100 == 0)
            log_squarenorm(create_log_prefix_find_max(i) + "v2: ", v2);
        //Normalize v2
        sax(&v2, { 1. / sqrt(squarenorm(v2)), 0. }, v2);
        if(i % 100 == 0)
            log_squarenorm(create_log_prefix_find_max(i) + "v2: ", v2);
        //Check whether the algorithm converged
        if(i % RESID_CHECK_FREQUENCY == 0) {
            saxpy(&v1, { -1., 0. }, v1, v2);   //Here I can use v1 as container, the important is not
                                               //to modify v2 that will be copied in v1
            log_squarenorm(create_log_prefix_find_max(i) + "v1: ", v1);
            resid = sqrt(squarenorm(v1));

            logger.debug() << create_log_prefix_find_max(i) << "resid: " << std::setprecision(8) << resid;

            if(resid < prec) {
                A(&v1, gf, v2, &additionalParameters);
                scalar_product(&max, v2, v1);
                hmc_complex result = max.get();
                logger.debug() << "max.im = " << result.im;
                if((result.im) > 1.e-12) {
                    logger.fatal() << "Power Method found complex eigenvalue! result.im = " << result.im;
                    throw solvers::SolverStuck(i, __FILE__, __LINE__);
                }
                //Here we are sure the eigenvalue is correctly found, then we get the duration
                const uint64_t duration = timer.getTime();
                logger.debug() << "Find_max_eig completed in " << duration / 1000.f << " ms. Performed " << i << " iterations (resid = " << resid << ").";
                return result.re;
            }
        }
        copyData(&v1, v2);
    }

    logger.fatal() << "Power Method failed in finding max_eig in " << parametersInterface.getFindMinMaxMaxValue() << " iterations. Last resid: " << resid;
    throw solvers::SolverDidNotSolve(parametersInterface.getFindMinMaxMaxValue(), __FILE__, __LINE__);

}

hmc_float physics::algorithms::find_min_eigenvalue(const physics::fermionmatrix::Fermionmatrix_stagg_eo& A, const physics::lattices::Gaugefield& gf,
                                                   const hardware::System& system, physics::InterfacesHandler& interfacesHandler, hmc_float prec,
                                                   const physics::AdditionalParameters& additionalParameters)
{
    if(additionalParameters.getConservative())
        return A.getThresholdForMinimumEigenvalue(additionalParameters.getMass());

    if(!(A.isHermitian()))
        throw std::invalid_argument("Unable to deal with non-hermitian matrices in find_max_eigenvalue!");

    hmc_float max = find_max_eigenvalue(A, gf, system, interfacesHandler, prec, additionalParameters);

    return find_min_knowing_max(max, A, gf, system, interfacesHandler, prec, additionalParameters);

}

void physics::algorithms::find_maxmin_eigenvalue(hmc_float& max, hmc_float& min, const physics::fermionmatrix::Fermionmatrix_stagg_eo& A,
                                                 const physics::lattices::Gaugefield& gf, const hardware::System& system,
                                                 physics::InterfacesHandler& interfacesHandler, hmc_float prec,
                                                 const physics::AdditionalParameters& additionalParameters)
{
    //This timer is to know how long this function takes
    klepsydra::Monotonic timer;

    max = find_max_eigenvalue(A, gf, system, interfacesHandler, prec, additionalParameters);

    if(additionalParameters.getConservative()){
        min = A.getThresholdForMinimumEigenvalue(additionalParameters.getMass());
        max *= 1.05;
    }else{
        min = find_min_knowing_max(max, A, gf, system, interfacesHandler, prec, additionalParameters);
    }

    //Here we are sure the eigenvalue is correctly found, then we get the duration
    const uint64_t duration = timer.getTime();
    logger.debug() << "Find_maxmin_eig completed in " << duration / 1000.f << " ms.";

}

static hmc_float find_min_knowing_max(const hmc_float max, const physics::fermionmatrix::Fermionmatrix_stagg_eo& A, const physics::lattices::Gaugefield& gf,
                                      const hardware::System& system, physics::InterfacesHandler& interfacesHandler,
                                      hmc_float prec, const physics::AdditionalParameters& additionalParameters)
{
    using namespace physics::lattices;
    using namespace physics::algorithms;

    Scalar<hmc_complex> min(system);
    hmc_float resid;

    //This timer is to know how long this function takes
    klepsydra::Monotonic timer;

    //This field is the starting point and it must be random (we have to be sure
    //to have a non zero component along the eigenvectors referring to the smallest eigenvalue)
    Staggeredfield_eo v1(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
    pseudo_randomize<Staggeredfield_eo, su3vec>(&v1, 321);
    log_squarenorm(create_log_prefix_find_min(0) + "v1 (initial): ", v1);
    sax(&v1, { 1. / sqrt(squarenorm(v1)), 0. }, v1);   //v1 is now normalized
    //Auxiliary field
    Staggeredfield_eo v2(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
    //How often to check resid
    const physics::algorithms::MinMaxEigenvalueParametersInterface & parametersInterface = interfacesHandler.getMinMaxEigenvalueParametersInterface();
    const int RESID_CHECK_FREQUENCY = parametersInterface.getFindMinMaxIterationBlockSize();

    log_squarenorm(create_log_prefix_find_min(0) + "v1 (initial): ", v1);
    log_squarenorm(create_log_prefix_find_min(0) + "v2 (initial) [not-initialized]: ", v2);

    for (unsigned int i = 0; i < parametersInterface.getFindMinMaxMaxValue(); i++) {
        //Apply (max-A) onto v1
        A(&v2, gf, v1, &additionalParameters);
        saxpby(&v2, { max, 0. }, v1, { -1., 0. }, v2);   //Now in v2 there is (max-A)*v1
        if(i % 100 == 0)
            log_squarenorm(create_log_prefix_find_min(i) + "v2: ", v2);
        //Normalize v2
        sax(&v2, { 1. / sqrt(squarenorm(v2)), 0. }, v2);
        if(i % 100 == 0)
            log_squarenorm(create_log_prefix_find_min(i) + "v2: ", v2);
        //Check whether the algorithm converged
        if(i % RESID_CHECK_FREQUENCY == 0) {
            saxpy(&v1, { -1., 0. }, v1, v2);   //Here I can use v1 as container, the important is not
                                               //to modify v2 that will be copied in v1
            log_squarenorm(create_log_prefix_find_min(i) + "v1: ", v1);
            resid = sqrt(squarenorm(v1));

            logger.debug() << create_log_prefix_find_min(i) << "resid: " << std::setprecision(8) << resid;

            if(resid < prec) {
                //Apply (max-A) onto v2
                A(&v1, gf, v2, &additionalParameters);
                saxpby(&v1, { max, 0. }, v2, { -1., 0. }, v1);   //Now in v1 there is (max-A)*v2
                scalar_product(&min, v2, v1);
                hmc_complex result = min.get();
                logger.debug() << "min.im = " << result.im;
                if(abs(result.im) > 1.e-12) {
                    logger.fatal() << "Power Method found complex eigenvalue!";
                    throw solvers::SolverStuck(i, __FILE__, __LINE__);
                }
                //Here we are sure the eigenvalue is correctly found, then we get the duration
                const uint64_t duration = timer.getTime();
                logger.debug() << "Find_min_eig completed in " << duration / 1000.f << " ms. Performed " << i << " iterations (resid = " << resid << ").";
                return max - result.re;
            }
        }
        copyData(&v1, v2);
    }

    logger.fatal() << "Power Method failed in finding min_eig in " << parametersInterface.getFindMinMaxMaxValue() << " iterations. Last resid: " << resid;
    throw solvers::SolverDidNotSolve(parametersInterface.getFindMinMaxMaxValue(), __FILE__, __LINE__);

}

static std::string create_log_prefix_find(std::string name, int number) noexcept
{
    using namespace std;
    string separator_big = "\t";
    string separator_small = " ";
    string label = "FIND";

    stringstream strnumber;
    strnumber.fill('0');
    /// @todo this should be length(findminmax_max)
    strnumber.width(6);
    strnumber << right << number;
    stringstream outfilename;
    outfilename << separator_big << label << separator_small << "[" << name << "]" << separator_small << "[" << strnumber.str() << "]:" << separator_big;
    string outputfile = outfilename.str();
    return outputfile;
}

static std::string create_log_prefix_find_max(int number) noexcept
{
    return create_log_prefix_find("MAX", number);
}

static std::string create_log_prefix_find_min(int number) noexcept
{
    return create_log_prefix_find("MIN", number);
}

