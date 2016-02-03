/** @file
 * Declaration of the algorithm to find min and max eigenvalues of an operator
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

#ifndef _PHYSICS_ALGORITHMS_FINDMINMAX_
#define _PHYSICS_ALGORITHMS_FINDMINMAX_

#include "../lattices/gaugefield.hpp"
#include "../fermionmatrix/fermionmatrix_stagg.hpp"
#include "../interfacesHandler.hpp"

namespace physics {

    namespace algorithms {

        /**
         * This function returns the maximum eigenvalue of the operator A up to some specified precision.
         *  @param A The operator whose maximum eigenvalue has to be found
         *  @param gf The gaugefield on which A depends
         *  @param system The system it is operating on
         *  @param prec The precision up to which the maximum eigenvalue is found
         *
         *  @note The algorithm here implemented is the so called "Power Method" (with the
         *        Rayleigh quotient estimate to get the largest eigenvalue of A). Hence to be sure
         *        that the algorithm works, the eigenvalues must be real. If this is not the case,
         *        the algorithm could give meaningless results or not converge. Actually, even if
         *        the eigenvalues are real, if there is not only ONE dominant eigenvalue, then
         *        the algorithm won't converge (for example if there are two opposite biggest eigenvalues).
         *        However, to be honest, in the RHMC with staggered fermions, one needs to find the
         *        maximum and the minimum eigenvalues of the matrix MdagM that is hermitian and then with
         *        real positive eigenvalue. For this reason, this function will calculate the eigenvalue
         *        if and only if the operator A is hermitian.
         */
        hmc_float find_max_eigenvalue(const physics::fermionmatrix::Fermionmatrix_stagg_eo& A, const physics::lattices::Gaugefield& gf,
                                      const hardware::System& system, physics::InterfacesHandler& interfacesHandler, hmc_float prec,
                                      const physics::AdditionalParameters& additionalParameters);

        /**
         * This function returns the minimum eigenvalue of the operator A up to some specified precision.
         *  @param A The operator whose minimum eigenvalue has to be found
         *  @param gf The gaugefield on which A depends
         *  @param system The system it is operating on
         *  @param prec The precision up to which the minimum eigenvalue is found
         *  @param conservative If true, the method returns the return value of the method getThresholdForMinimumEigenvalue
         *                      of the operator A. If false, then the minimum eigenvalue is calculated.
         *
         *  @note The note of the function find_max_eigenvalue is still valid for this function.
         *        Here to estimate the minimum eigenvalue, the maximum one is at first found, and then
         *        we find the maximum eigenvalue (max) of the operator (max*ID - A) where ID is the identity
         *        matrix. This, subtracted to max, will give the minimum eigenvalue of A.
         *        Furthermore, in the RHMC with staggered fermions, one needs to find the
         *        maximum and the minimum eigenvalues of the matrix MdagM that, beyond being hermitian,
         *        has all eigenvalues bigger than or equal to the mass squared of the quarks. This means
         *        that one could be conservative and return the mass^2 as minimum eigenvalue.
         */
        hmc_float find_min_eigenvalue(const physics::fermionmatrix::Fermionmatrix_stagg_eo& A, const physics::lattices::Gaugefield& gf,
                                      const hardware::System& system, physics::InterfacesHandler& interfacesHandler, hmc_float prec,
                                      const physics::AdditionalParameters& additionalParameters);

        /**
         * This function finds both the minimum and the maximum eigenvalue of the operator A
         * up to some specified precision.
         *  @param max The maximum eigenvalue
         *  @param min The minimum eigenvalue
         *  @param A The operator whose minimum eigenvalue has to be found
         *  @param gf The gaugefield on which A depends
         *  @param system The system it is operating on
         *  @param prec The precision up to which the minimum eigenvalue is found
         *  @param conservative If true, the method sets the return value of the method getThresholdForMinimumEigenvalue
         *                      of the operator A as minimum eigenvalue. If false, then the minimum eigenvalue is calculated.
         *                      The maximum eigenvalue is calculated and then increased vy 5%.
         *
         *  @note The note of the functions find_max_eigenvalue and find_min_eigenvalue are
         *        still valid for this function.
         */
        void find_maxmin_eigenvalue(hmc_float& max, hmc_float& min, const physics::fermionmatrix::Fermionmatrix_stagg_eo& A,
                                    const physics::lattices::Gaugefield& gf, const hardware::System& system,  physics::InterfacesHandler& interfacesHandler,
                                    hmc_float prec, const physics::AdditionalParameters& additionalParameters);

    }

}

#endif /* _PHYSICS_ALGORITHMS_FINDMINMAX_ */
