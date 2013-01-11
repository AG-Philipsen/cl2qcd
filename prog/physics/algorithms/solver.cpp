/** @file
 * Implementation of the solver algorithms
 *
 * (c) 2012-2013 Christopher Pinke <pinke@th.uni-frankfurt.de>
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "solver.hpp"

#include "../../logger.hpp"
#include "../../operations_complex.h"
#include "../lattices/util.hpp"

static std::string create_solver_stuck_message(int iterations);
/**
 * A "save" version of the bicgstab algorithm.
 *
 * This is chosen if "bicgstab_save" is selected in the inputfile
 */
static int bicgstab_save(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, hmc_float prec);

/**
 * BICGstab implementation with a different structure than "save" one, similar to tmlqcd. This should be the default bicgstab.
 * In particular this version does not perform the check if the "real" residuum is sufficiently small!
 */
static int bicgstab_fast(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, hmc_float prec);

physics::algorithms::solvers::SolverStuck::SolverStuck(int iterations, std::string filename, int linenumber) : SolverException(create_solver_stuck_message(iterations), iterations, filename, linenumber) { };

static std::string create_solver_stuck_message(int iterations)
{
	std::ostringstream tmp;
	tmp << "Solver got stuck after " << iterations << " iterations";
	return tmp.str();
}

int physics::algorithms::solvers::bicgstab(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, hmc_float prec)
{
	auto params = system.get_inputparameters();

	//"save" version, with comments. this is called if "bicgstab_save" is choosen.
	if (params.get_solver() == meta::Inputparameters::bicgstab_save) {
		return bicgstab_save(x, A, gf, b, system, prec);
	} else { /*if (get_parameters().get_solver() == meta::Inputparameters::bicgstab)*/
		// NOTE: I commented out the if, since one runs into trouble if one uses CG in the HMC.
		// Then, no bicgstab type is chosen, however, one still uses it "hardcoded".
		// Then one gets a fatal, which is not really meaningful. In this way, it is like in the eo case.
		return bicgstab_fast(x, A, gf, b, system, prec);
	}
}

static int bicgstab_save(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, hmc_float prec)
{
	// TODO this function often contains -1 in the comment but 1 in the code...
	using physics::lattices::Spinorfield;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	auto params = system.get_inputparameters();
	hmc_float resid;
	hmc_complex alpha, omega, rho;

	Spinorfield v(system);
	Spinorfield p(system);
	Spinorfield rn(system);
	Spinorfield rhat(system);
	Spinorfield s(system);
	Spinorfield t(system);
	Spinorfield aux(system);

	int iter;
	for(iter = 0; iter < params.get_cgmax(); iter++) {
		if(iter % params.get_iter_refresh() == 0) {
			v.zero();
			p.zero();

			//initial r_n
			f(&rn, gf, *x);
			saxpy(&rn, {1., 0}, rn, b);
			//rhat = r_n
			copyData(&rhat, rn);
			//set some constants to 1
			alpha = {1., 0.};
			omega = {1., 0.};
			rho = {1., 0.};
		}
		//rho_next = (rhat, rn)
		hmc_complex rho_next = scalar_product(rhat, rn);

		//check if algorithm is stuck
		//if rho is too small the algorithm will get stuck and will never converge!!
		if(abs(rho_next.re) < 1e-25 && abs(rho_next.im) < 1e-25 ) {
			//print the last residuum
			logger.fatal() << "\t\t\tsolver stuck at resid:\t" << resid;
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

		//v = A*p
		f(&v, gf, p);
		//tmp1 = (rhat, v)
		tmp1 = scalar_product(rhat, v);
		//alpha = rho/tmp1 = (..)/(rhat, v)
		alpha = complexdivide(rho, tmp1);
		//s = - alpha * v - r_n
		saxpy(&s, alpha, v, rn);
		//t = A s
		f(&t, gf, s);
		//tmp1 = (t, s)
		tmp1 = scalar_product(t, s);
		//!!CP: this can also be global_squarenorm, but one needs a complex number here
		//tmp2 = (t,t)
		tmp2 = scalar_product(t, t);
		//omega = tmp1/tmp2 = (t,s)/(t,t)
		omega = complexdivide(tmp1, tmp2);
		//r_n = - omega*t - s
		saxpy(&rn, omega, t, s);
		//inout = alpha*p + omega * s + inout
		saxsbypz(x, alpha, p, omega, s, *x);
		//resid = (rn,rn)
		resid = squarenorm(rn);

		logger.debug() << "resid: " << resid;
		//test if resid is NAN
		if(resid != resid) {
			logger.fatal() << "\tNAN occured in bicgstab!";
			// TODO throw a NAN exception
			throw SolverStuck(iter, __FILE__, __LINE__);
		}
		if(resid < prec) {
			//aux = A inout
			f(&aux, gf, *x);
			//aux = -aux + source
			saxpy(&aux, {1., 0.}, aux, *x);
			//trueresid = (aux, aux)
			hmc_float trueresid = squarenorm(aux);
			logger.debug() << "\tsolver converged! true resid:\t" << trueresid;
			if(trueresid < prec)
				return iter;
		}
	}
	throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}

static int bicgstab_fast(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, hmc_float prec, hmc_float kappa, hmc_float)
{
	using physics::lattices::Spinorfield;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	auto params = system.get_inputparameters();
	hmc_float resid;
	hmc_complex rho;

	Spinorfield p(system);
	Spinorfield rn(system);
	Spinorfield rhat(system);
	Spinorfield v(system);
	Spinorfield s(system);
	Spinorfield t(system);

	int iter;
	for(iter = 0; iter < params.get_cgmax(); iter++) {
		if(iter % params.get_iter_refresh() == 0) {
			//initial r_n, saved in p
			f(&rn, gf, *x);
			saxpy(&p, {1.0, 0}, rn, b);
			//rhat = p
			copyData(&rhat, p);
			//r_n = p
			copyData(&rn, p);
			//rho = (rhat, rn)
			rho = scalar_product(rhat, rn);
		}
		//resid = (rn,rn)
		resid = squarenorm(rn);
		//test if resid is NAN
		if(resid != resid) {
			logger.fatal() << "\tNAN occured in bicgstab!";
			throw SolverStuck(iter, __FILE__, __LINE__);
		}
		if(resid < prec) {
			return iter;
		}
		//v = A*p
		f(&v, gf, p);
		//tmp1 = (rhat, v)
		hmc_complex tmp1 = scalar_product(rhat, v);
		//alpha = rho/tmp1 = (rhat, rn)/(rhat, v)
		hmc_complex alpha = complexdivide(rho, tmp1);
		//s = - alpha * v - r_n
		saxpy(&s, alpha, v, rn);
		//t = A s
		f(&t, gf, s);
		//tmp1 = (t, s)
		tmp1 = scalar_product(t, s);
		//!!CP: this can also be global_squarenorm, but one needs a complex number here
		//tmp2 = (t,t)
		hmc_complex tmp2 = scalar_product(t, t);
		//omega = tmp1/tmp2 = (t,s)/(t,t)
		hmc_complex omega = complexdivide(tmp1, tmp2);
		//inout = alpha*p + omega * s + inout
		saxsbypz(x, alpha, p, omega, s, *x);
		//r_n = - omega*t - s
		saxpy(&rn, omega, t, s);
		//rho_next = (rhat, rn)
		hmc_complex rho_next = scalar_product(rhat, rn);
		//check if algorithm is stuck
		//if rho is too small the algorithm will get stuck and will never converge!!
		if(abs(rho_next.re) < 1e-25 && abs(rho_next.im) < 1e-25 ) {
			//print the last residuum
			logger.fatal() << "\t\t\tsolver stuck at resid:\t" << resid;
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
		tmp2 = complexsubtract( {0., 0.}, tmp1);
		//p = beta*p + tmp2*v + r_n = beta*p - beta*omega*v + r_n
		saxsbypz(&p, beta, p, tmp2, v, rn);
		//rho_next = rho
		rho = rho_next;
	}
	throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}
