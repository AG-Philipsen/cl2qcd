/** @file
 * Declaration of the solver algorithms
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

#ifndef _PHYSICS_ALGORITHMS_SOLVER_SHIFTED_
#define _PHYSICS_ALGORITHMS_SOLVER_SHIFTED_

#include "../fermionmatrix/fermionmatrix_stagg.hpp"
#include "../interfacesHandler.hpp"

namespace physics {
    namespace algorithms {
        namespace solvers {

            /**
             * Solve the linear system (A + sigma) * x = b for x and a whole set of values of sigma,
             * using the CG-M algorithm (multi-shifted inverter) for even-odd preconditioned b and x.
             * For further details, see B.Jegerlehner arXiv:hep-lat/9612014v1, even if here there are some missing information
             * and also some wrong indeces. Then compare with "Accurate conjugate gradient methods for families of
             * shifted systems" that is a work by Jasper van den Eshof and Gerard L. G. Sleijpen of December 2003:
             * here the notation is completely different, but one can build the complete CG-M algorithm. This was done
             * and summarized in the file CG-M.pdf (see Feature #482)
             *
             * \return The number of iterations performed
             *
             * \exception SolverStuck if the solver gets stuck. Contains information on performed iterations
             *
             * \exception SolverDidNotSolve if the solver did not solve (hit iteration limit).
             */
            int cg_m(const std::vector<std::shared_ptr<physics::lattices::Staggeredfield_eo> > x, const physics::fermionmatrix::Fermionmatrix_stagg_eo& A,
                     const physics::lattices::Gaugefield& gf, const std::vector<hmc_float> sigma,
                     const physics::lattices::Staggeredfield_eo& b, const hardware::System& system,
                     physics::InterfacesHandler& interfacesHandler, hmc_float prec, const physics::AdditionalParameters& additionalParameters);


            template<typename FERMIONFIELD, typename FERMIONMATRIX>
            class SolverShifted {
                public:
                    SolverShifted(const std::vector<std::shared_ptr<FERMIONFIELD> >, const FERMIONMATRIX&,
                                  const physics::lattices::Gaugefield&, const std::vector<hmc_float>, const FERMIONFIELD&,
                                  const hardware::System&, physics::InterfacesHandler&, hmc_float, const physics::AdditionalParameters&);
                    ~SolverShifted(){};

                    SolverShifted() = delete;
                    SolverShifted(const SolverShifted&) = delete;
                    SolverShifted& operator=(const SolverShifted&) = delete;
                    unsigned int getNumberOfIterationsDone();
                    const std::vector<std::shared_ptr<FERMIONFIELD> > solve();

                private:
                    void setInitialConditions();
                    void updateBetaScalar();
                    void updateAuxiliaryFieldR();
                    void updateAlphaScalar();
                    void updateAuxiliaryFieldP();
                    void updateVectorQuantities();
                    void updateSingleFieldOfSolution(unsigned int);
                    void updateSingleFieldOfAuxiliaryFieldPs(unsigned int);
                    void checkFiledsSquarenormsForPossibleNaN();
                    void updateQuantitiesForFollowingIteration();
                    void calculateResiduumSingleEquation(bool, unsigned int = 0);
                    void checkIfSingleEquationConverged(unsigned int);
                    bool hasSingleSystemConverged(unsigned int);
                    bool hasSystemConverged();
                    void makePerformanceReport();
                    void resetNoWarmupTimer();
                    std::string createLogPrefix();
                    void debugLogSquarenormSetOfFields(const std::string&, const std::vector<std::shared_ptr<physics::lattices::Staggeredfield_eo> >, const int);
                    void debugCompareSquarenormsOfResultFieldsBetweenBeginAndEndOfIteration();
                    void debugMakeReportOfFieldsSquarenorm();
                    void debugCalculateSquarenormOfResultField();


                    //Parameter read from outside
                    const std::vector<std::shared_ptr<FERMIONFIELD> > x;
                    const FERMIONMATRIX& A;
                    const physics::lattices::Gaugefield& gf;
                    const std::vector<hmc_float> sigma;
                    const FERMIONFIELD& b;
                    const hardware::System& system;
                    hmc_float solverPrecision;
                    const physics::AdditionalParameters& additionalParameters;

                    //Variable of the solver
                    const physics::algorithms::SolversParametersInterface & parametersInterface;
                    bool hasSystemBeSolved;
                    unsigned int numberOfEquations;
                    unsigned int iterationNumber;
                    hmc_float residuumValue;
                    klepsydra::Monotonic timer; /// @TODO: start timers synchronized with device(s)
                    klepsydra::Monotonic timer_noWarmup;
                    bool USE_ASYNC_COPY;     /// @todo make configurable from outside
                    unsigned int MINIMUM_ITERATIONS;  /// @todo make configurable from outside


                    //Auxiliary staggered fields
                    const FERMIONFIELD r;
                    const FERMIONFIELD p;
                    std::vector<std::shared_ptr<FERMIONFIELD> > ps;

                    //Auxiliary scalar vectors
                    const physics::lattices::Vector<hmc_float> zeta_prev;    //This is zeta at the step iter-1
                    const physics::lattices::Vector<hmc_float> zeta;         //This is zeta at the step iter
                    const physics::lattices::Vector<hmc_float> zeta_foll;    //This is zeta at the step iter+1
                    physics::lattices::Vector<hmc_float> alpha_vec;
                    physics::lattices::Vector<hmc_float> beta_vec;
                    physics::lattices::Vector<hmc_float> shift;              //This is to store constants sigma
                    std::vector<bool> single_system_converged;               //This is to stop calculation on single system
                    std::vector<uint> single_system_iter;                    //This is to calculate performance properly
                    std::vector<hmc_float> resultSquarenorm;

                    //Auxiliary scalars
                    const physics::lattices::Scalar<hmc_float> alpha_scalar_prev;       //This is alpha_scalar at the step iter-1
                    const physics::lattices::Scalar<hmc_float> alpha_scalar;            //This is alpha_scalar at the step iter
                    const physics::lattices::Scalar<hmc_float> beta_scalar_prev;        //This is beta_scalar at the step iter-1
                    const physics::lattices::Scalar<hmc_float> beta_scalar;             //This is beta_scalar at the step iter
                    const physics::lattices::Scalar<hmc_float> zero;

                    //Auxiliary containers for temporary saving
                    const FERMIONFIELD v;                               //This is to store A.p
                    const physics::lattices::Scalar<hmc_float> tmp1;    //This is to store (r,r) before updating r
                    const physics::lattices::Scalar<hmc_float> tmp2;    //This is to store (r,r) after updating r
                    const physics::lattices::Scalar<hmc_float> tmp3;    //This is to store (p,v) as Scalar

                    //Only if merged kernels are used (this is to store the single eq. residuum all at once)
                    std::unique_ptr<physics::lattices::Vector<hmc_float> > single_eq_resid;
                    std::unique_ptr<std::vector<hmc_float> > single_eq_resid_host;

            };

        }
    }
}
#endif /* _PHYSICS_ALGORITHMS_SOLVER_SHIFTED_ */
