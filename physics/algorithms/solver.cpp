/** @file
 * Implementation of the solver algorithms
 *
 * Copyright (c) 2012-2013 Christopher Pinke <pinke@th.uni-frankfurt.de>
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

#include "solver.hpp"

#include "../../host_functionality/logger.hpp"
#include "../../common_header_files/operations_complex.h"
#include "../../meta/type_ops.hpp"
#include "../lattices/util.hpp"
#include "../lattices/scalar_complex.hpp"
#include <cmath>
#include <sstream>

static std::string create_solver_stuck_message(int iterations);
/**
 * A "save" version of the bicgstab algorithm.
 * It is explicitely checked if the true residuum also fullfills the break condition.
 * This is chosen if "bicgstab_save" is selected.
 */
static int bicgstab_save(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, hmc_float prec);
static int bicgstab_save(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, hmc_float prec);

/**
 * BICGstab implementation with a different structure than "save" one, similar to tmlqcd. This should be the default bicgstab.
 * In particular this version does not perform the check if the "real" residuum is sufficiently small!
 */
static int bicgstab_fast(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, hmc_float prec);
static int bicgstab_fast(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, hmc_float prec);

namespace {

int cg_singledev(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, const hmc_float prec);
int cg_multidev(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, const hmc_float prec);

}

physics::algorithms::solvers::SolverStuck::SolverStuck(int iterations, std::string filename, int linenumber) : SolverException(create_solver_stuck_message(iterations), iterations, filename, linenumber) { }

static std::string create_solver_stuck_message(int iterations)
{
	std::ostringstream tmp;
	tmp << "Solver got stuck after " << iterations << " iterations";
	return tmp.str();
}

static std::string create_log_prefix_solver(std::string name, int number) noexcept;
static std::string create_log_prefix_cg(int number) noexcept;
static std::string create_log_prefix_bicgstab(int number) noexcept;

int physics::algorithms::solvers::bicgstab(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, hmc_float prec)
{
	const auto & params = system.get_inputparameters();

	if (params.get_solver() == meta::Inputparameters::bicgstab_save) {
		return bicgstab_save(x, A, gf, b, system, prec);
	} else {
		return bicgstab_fast(x, A, gf, b, system, prec);
	}
}

static int bicgstab_save(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, const hmc_float prec)
{
	// @todo this function often contains -1 in the comment but 1 in the code...
	using physics::lattices::Spinorfield;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	const auto & params = system.get_inputparameters();

	/// @todo start timer synchronized with device(s)
	klepsydra::Monotonic timer;

	hmc_float resid;
	hmc_complex alpha, omega, rho{std::nan(""), std::nan("")};

	const Spinorfield v(system);
	const Spinorfield p(system);
	const Spinorfield rn(system);
	const Spinorfield rhat(system);
	const Spinorfield s(system);
	const Spinorfield t(system);
	const Spinorfield aux(system);

	int iter=0;
	log_squarenorm(create_log_prefix_bicgstab(iter) + "b: ", b);
	log_squarenorm(create_log_prefix_bicgstab(iter) + "x (initial): ", *x);

	for(iter = 0; iter < params.get_cgmax(); iter++) {
		if(iter % params.get_iter_refresh() == 0) {
			v.zero();
			p.zero();

			//initial r_n
			f(&rn, gf, *x);
			log_squarenorm(create_log_prefix_bicgstab(iter) + "rn: ", rn);
			saxpy(&rn, {1., 0}, rn, b);
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
		f(&v, gf, p);
		log_squarenorm(create_log_prefix_bicgstab(iter) + "v: ", v);

		//tmp1 = (rhat, v)
		tmp1 = scalar_product(rhat, v);
		//alpha = rho/tmp1 = (..)/(rhat, v)
		alpha = complexdivide(rho, tmp1);
		//s = - alpha * v - r_n
		saxpy(&s, alpha, v, rn);
		log_squarenorm(create_log_prefix_bicgstab(iter) + "s: ", s);

		//t = A s
		f(&t, gf, s);
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
		logger.debug() << create_log_prefix_bicgstab(iter) <<  "resid: " << resid;

		//test if resid is NAN
		if(resid != resid) {
		  logger.fatal() << create_log_prefix_bicgstab(iter) << "\tNAN occured!";
		  throw SolverStuck(iter, __FILE__, __LINE__);
		}

		if(resid < prec) {
			//aux = A inout
			f(&aux, gf, *x);
			log_squarenorm(create_log_prefix_bicgstab(iter) + "aux: ", aux);
		
			//aux = -aux + source
			saxpy(&aux, {1., 0.}, aux, b);
			log_squarenorm(create_log_prefix_bicgstab(iter) + "aux: ", aux);

			//trueresid = (aux, aux)
			hmc_float trueresid = squarenorm(aux);
			if(trueresid < prec){
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

	logger.fatal() << create_log_prefix_bicgstab(iter) << "Solver did not solve in " << params.get_cgmax() << " iterations. Last resid: " << resid;
	throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}

static int bicgstab_fast(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, const hmc_float prec)
{
	using physics::lattices::Spinorfield;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	const auto & params = system.get_inputparameters();

	/// @todo start timer synchronized with device(s)
	klepsydra::Monotonic timer;

	hmc_float resid;
	hmc_complex rho;

	const Spinorfield p(system);
	const Spinorfield rn(system);
	const Spinorfield rhat(system);
	const Spinorfield v(system);
	const Spinorfield s(system);
	const Spinorfield t(system);

	int iter=0;
	log_squarenorm(create_log_prefix_bicgstab(iter) + "b: ", b);
	log_squarenorm(create_log_prefix_bicgstab(iter) + "x (initial): ", *x);

	for(iter = 0; iter < params.get_cgmax(); iter++) {
		if(iter % params.get_iter_refresh() == 0) {
			//initial r_n, saved in p
			f(&rn, gf, *x);
			log_squarenorm(create_log_prefix_bicgstab(iter) + "rn: ", rn);
			saxpy(&p, {1.0, 0}, rn, b);
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
		f(&v, gf, p);
		log_squarenorm(create_log_prefix_bicgstab(iter) + "v: ", v);

		//tmp1 = (rhat, v)
		hmc_complex tmp1 = scalar_product(rhat, v);
		//alpha = rho/tmp1 = (rhat, rn)/(rhat, v)
		hmc_complex alpha = complexdivide(rho, tmp1);

		//s = - alpha * v - r_n
		saxpy(&s, alpha, v, rn);
		log_squarenorm(create_log_prefix_bicgstab(iter) + "s: ", s);

		//t = A s
		f(&t, gf, s);
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
		if(std::abs(rho_next.re) < 1e-25 && std::abs(rho_next.im) < 1e-25 ) {
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
		tmp2 = complexsubtract( {0., 0.}, tmp1);

		//p = beta*p + tmp2*v + r_n = beta*p - beta*omega*v + r_n
		saxsbypz(&p, beta, p, tmp2, v, rn);
		log_squarenorm(create_log_prefix_bicgstab(iter) + "p: ", p);

		//rho_next = rho
		rho = rho_next;
	}

	logger.fatal() << create_log_prefix_bicgstab(iter) << "Solver did not solve in " << params.get_cgmax() << " iterations. Last resid: " << resid;
	throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}

int physics::algorithms::solvers::cg(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, const hmc_float prec)
{
	using physics::lattices::Spinorfield;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	const auto & params = system.get_inputparameters();

	/// @todo start timer synchronized with device(s)
	klepsydra::Monotonic timer;

	const Spinorfield rn(system);
	const Spinorfield p(system);
	const Spinorfield v(system);

	hmc_complex rho_next;
	hmc_float resid;
	int iter = 0;

	log_squarenorm(create_log_prefix_cg(iter) + "b: ", b);
	log_squarenorm(create_log_prefix_cg(iter) + "x (initial): ", *x);

	// NOTE: here, most of the complex numbers may also be just hmc_floats. However, for this one would need some add. functions...
	for(iter = 0; iter < params.get_cgmax(); iter ++) {
		hmc_complex omega;
		if(iter % params.get_iter_refresh() == 0) {
			//rn = A*inout
			f(&rn, gf, *x);
			log_squarenorm(create_log_prefix_cg(iter) + "rn: ", rn);
			//rn = source - A*inout
			saxpy(&rn, {1., 0.}, rn, b);
			log_squarenorm(create_log_prefix_cg(iter) + "rn: ", rn);
			//p = rn
			copyData(&p, rn);
			log_squarenorm(create_log_prefix_cg(iter) + "p: ", p);
			//omega = (rn,rn)
			omega = scalar_product(rn, rn);
		} else {
			//update omega
			omega = rho_next;
		}
		//v = A pn
		f(&v, gf, p);
		log_squarenorm(create_log_prefix_cg(iter) + "v: ", v);

		//alpha = (rn, rn)/(pn, Apn) --> alpha = omega/rho
		hmc_complex rho = scalar_product(p, v);
		hmc_complex alpha = complexdivide(omega, rho);
		hmc_complex tmp1 = complexsubtract( {0., 0.}, alpha);

		//xn+1 = xn + alpha*p = xn - tmp1*p = xn - (-tmp1)*p
		saxpy(x, tmp1, p, *x);
		log_squarenorm(create_log_prefix_cg(iter) + "x: ", *x);

		//rn+1 = rn - alpha*v -> rhat
		saxpy(&rn, alpha, v, rn);
		log_squarenorm(create_log_prefix_cg(iter) + "rn: ", rn);

		//calc residuum
		//NOTE: for beta one needs a complex number at the moment, therefore, this is done with "rho_next" instead of "resid"
		rho_next = scalar_product(rn, rn);
		resid = rho_next.re;

		logger.debug() << create_log_prefix_cg(iter) << "resid: " << resid;

		//test if resid is NAN
		if(resid != resid) {
		  logger.fatal() << create_log_prefix_cg(iter) << "NAN occured!";
		  throw SolverStuck(iter, __FILE__, __LINE__);
		}

		if(resid < prec){
		  logger.debug() << create_log_prefix_cg(iter) << "Solver converged in " << iter << " iterations! resid:\t" << resid;

		  // report on performance
		  if(logger.beInfo()) {
		    // we are always synchroneous here, as we had to recieve the residium from the device
		    const uint64_t duration = timer.getTime();
		    
		    // calculate flops
		    /**
		     * @todo this is not implemented since for non-eo this should not be of interest
		     */
		    // report performance
		    logger.info() << create_log_prefix_cg(iter) << "Solver completed in " << duration / 1000 << " ms. Performed " << iter << " iterations";
		  }
		  
		  //report on solution
		  log_squarenorm(create_log_prefix_cg(iter) + "x (final): ", *x);
		  return iter;
		}

		//beta = (rn+1, rn+1)/(rn, rn) --> alpha = rho_next/omega
		hmc_complex beta = complexdivide(rho_next, omega);

		//pn+1 = rn+1 + beta*pn
		hmc_complex tmp2 = complexsubtract( {0., 0.}, beta);
		saxpy(&p, tmp2, p, rn);
		log_squarenorm(create_log_prefix_cg(iter) + "p: ", p);
	}

	logger.fatal() << create_log_prefix_cg(iter) << "Solver did not solve in " << params.get_cgmax() << " iterations. Last resid: " << resid;
	throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}

int physics::algorithms::solvers::bicgstab(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, hmc_float prec)
{
	const auto & params = system.get_inputparameters();

	if (params.get_solver() == meta::Inputparameters::bicgstab_save) {
	  return bicgstab_save(x, A, gf, b, system, prec);
	} else { 
	  return bicgstab_fast(x, A, gf, b, system, prec);
	}
}

static int bicgstab_save(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, const hmc_float prec)
{
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;
	using namespace physics::lattices;

	// TODO start timer synchronized with device(s)
	klepsydra::Monotonic timer;

	auto & params = system.get_inputparameters();

	const Spinorfield_eo s(system);
	const Spinorfield_eo t(system);
	const Spinorfield_eo v(system);
	const Spinorfield_eo p(system);
	const Spinorfield_eo rn(system);
	const Spinorfield_eo rhat(system);
	const Spinorfield_eo aux(system);

	unsigned retests = 0;
	int cgmax = params.get_cgmax();

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
	int iter =  0;

	// report source and initial solution
	log_squarenorm(create_log_prefix_bicgstab(iter) + "b (initial): ", b);
	log_squarenorm(create_log_prefix_bicgstab(iter) + "x (initial): ", *x);

	// comments correspond to the bicgstab_fast version
	for(iter = 0; iter < cgmax; iter++) {
		if(iter % params.get_iter_refresh() == 0) {
			v.zero();
			p.zero();

			f(&rn, gf, *x);
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
		if(std::abs(test.re) < 1e-25 && std::abs(test.im) < 1e-25 ) {
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

		f(&v, gf, p);
		log_squarenorm(create_log_prefix_bicgstab(iter) + "v: ", v);

		scalar_product(&tmp1, rhat, v);
		divide(&alpha, rho, tmp1);

		saxpy(&s, alpha, v, rn);
		log_squarenorm(create_log_prefix_bicgstab(iter) + "s: ", s);

		f(&t, gf, s);
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

			f(&aux, gf, *x);
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
			    const unsigned refreshs = iter / params.get_iter_refresh() + 1;
			    const size_t mf_flops = f.get_flops();
			    
			    cl_ulong total_flops = 4 * get_flops<Spinorfield_eo, scalar_product>(system) + 4 * get_flops<hmc_complex, complexdivide>()
			      + 3 * get_flops<hmc_complex, complexmult>() + 2 * get_flops<Spinorfield_eo, saxsbypz>(system)
			      + 2 * mf_flops + 2 * get_flops<Spinorfield_eo, saxpy>(system)
			      + get_flops<Spinorfield_eo, squarenorm>(system);
			    total_flops *= iter;
			    
			    total_flops += refreshs * (mf_flops + get_flops<Spinorfield_eo, saxpy>(system));
			    total_flops += retests * (mf_flops + get_flops<Spinorfield_eo, saxpy>(system) + get_flops<Spinorfield_eo, squarenorm>(system));
			    
			    // report performanc
			    logger.info() << create_log_prefix_bicgstab(iter) << "BiCGstab_save completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f) << " Gflops. Performed " << iter << " iterations";
			  }
			  
			  log_squarenorm(create_log_prefix_bicgstab(iter) + "x (final): ", *x);
			  return iter;
			}
		}
	}

	logger.fatal() << create_log_prefix_bicgstab(iter) << "Solver did not solve in " << params.get_cgmax() << " iterations. Last resid: " << resid;
	throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}

static int bicgstab_fast(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, hmc_float prec)
{
	using namespace physics::lattices;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	// TODO start timer synchronized with device(s)
	klepsydra::Monotonic timer;

	auto & params = system.get_inputparameters();

	const Spinorfield_eo p(system);
	const Spinorfield_eo rn(system);
	const Spinorfield_eo rhat(system);
	const Spinorfield_eo s(system);
	const Spinorfield_eo t(system);
	const Spinorfield_eo v(system);

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
	int iter = 0;

	// report source and initial solution
	log_squarenorm(create_log_prefix_bicgstab(iter) + "b (initial): ", b);
	log_squarenorm(create_log_prefix_bicgstab(iter) + "x (initial): ", *x);

	for(iter = 0; iter < params.get_cgmax(); iter++) {
		if(iter % params.get_iter_refresh() == 0) {
			//initial r_n, saved in p
			f(&rn, gf, *x);
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
				const unsigned refreshs = iter / params.get_iter_refresh() + 1;
				const cl_ulong mf_flops = f.get_flops();

				cl_ulong total_flops = get_flops<Spinorfield_eo, squarenorm>(system) + 2 * mf_flops
				                       + 4 * get_flops<Spinorfield_eo, scalar_product>(system) + 4 * get_flops<hmc_complex, complexdivide>()
				                       + 2 * get_flops<Spinorfield_eo, saxpy>(system) + 2 * get_flops<Spinorfield_eo, saxsbypz>(system)
				                       + 3 * get_flops<hmc_complex, complexmult>();
				total_flops *= iter;

				total_flops += refreshs * (mf_flops + get_flops<Spinorfield_eo, saxpy>(system) + get_flops<Spinorfield_eo, scalar_product>(system));

				// report performanc
				logger.info() << create_log_prefix_bicgstab(iter) << "Solver completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f) << " Gflops. Performed " << iter << " iterations";
			}

			log_squarenorm(create_log_prefix_bicgstab(iter) + "x (final): ", *x);
			return iter;
		}
		//v = A*p
		f(&v, gf, p);
		log_squarenorm(create_log_prefix_bicgstab(iter) + "v: ", v);

		//tmp1 = (rhat, v)
		scalar_product(&tmp1, rhat, v);
		//alpha = rho/tmp1 = (rhat, rn)/(rhat, v)
		divide(&alpha, rho, tmp1);
		//s = - alpha * v - r_n
		saxpy(&s, alpha, v, rn);
		log_squarenorm(create_log_prefix_bicgstab(iter) + "s: ", s);

		//t = A s
		f(&t, gf, s);
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
		if(std::abs(test.re) < 1e-25 && std::abs(test.im) < 1e-25 ) {
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

	logger.fatal() << create_log_prefix_bicgstab(iter) << "Solver did not solve in " << params.get_cgmax() << " iterations. Last resid: " << resid;
	throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}

int physics::algorithms::solvers::cg(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, const hmc_float prec)
{
	if(system.get_devices().size() > 1) {
		return cg_multidev(x, f, gf, b, system, prec);
	} else {
		return cg_singledev(x, f, gf, b, system, prec);
	}
}


namespace {

void testIfResiduumIsNan(hmc_float resid, int iter)
{
	if(resid != resid) {
		logger.fatal() << create_log_prefix_cg(iter) << "NAN occured!";
		throw physics::algorithms::solvers::SolverStuck(iter, __FILE__, __LINE__);
	}
}

void reportPerformance_cg(int iter, const uint64_t duration, const uint64_t duration_noWarmup, cl_ulong total_flops, cl_ulong noWarmup_flops )
{
	logger.info() << create_log_prefix_cg(iter) << "CG completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f) << " Gflops. Performed " << iter << " iterations. Performance after warmup: " << (noWarmup_flops / duration_noWarmup / 1000.f) << " Gflops.";
}
	
int cg_singledev(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, const hmc_float prec)
{
	using namespace physics::lattices;

	const auto & params = system.get_inputparameters();

	/// @todo start timer synchronized with device(s)
	klepsydra::Monotonic timer;
	klepsydra::Monotonic timer_noWarmup;
	/// @todo make configurable from outside
	const int RESID_CHECK_FREQUENCY = params.get_cg_iteration_block_size();
	const bool USE_ASYNC_COPY = params.get_cg_use_async_copy();
	const int MINIMUM_ITERATIONS = params.get_cg_minimum_iteration_count();
	if(USE_ASYNC_COPY) {
		logger.warn() << "Asynchroneous copying in the CG is currently unimplemented!";
	}
	if(MINIMUM_ITERATIONS) {
		logger.warn() << "Minimum iterations set to " << MINIMUM_ITERATIONS << " -- should be used *only* for inverter benchmarking!";
	}

	const Spinorfield_eo p(system);
	const Spinorfield_eo rn(system);
	const Spinorfield_eo v(system);

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
	int iter = 0;

	// report source and initial solution
	log_squarenorm(create_log_prefix_cg(iter) + "b (initial): ", b);
	log_squarenorm(create_log_prefix_cg(iter) + "x (initial): ", *x);

	//NOTE: here, most of the complex numbers may also be just hmc_floats. However, for this one would need some add. functions...
	for(iter = 0; iter < params.get_cgmax() || iter < MINIMUM_ITERATIONS; iter ++) {
		if(iter % params.get_iter_refresh() == 0) {
			f(&rn, gf, *x); //rn = A*inout
			log_squarenorm(create_log_prefix_cg(iter) + "rn: ", rn);
			
			saxpy(&rn, one, rn, b); //rn = source - A*inout
			log_squarenorm(create_log_prefix_cg(iter) + "rn: ", rn);
			
			copyData(&p, rn); //p = rn
			log_squarenorm(create_log_prefix_cg(iter) + "p: ", p);
			
			scalar_product(&omega, rn, rn); //omega = (rn,rn)
		} else {
			copyData(&omega, rho_next);
		}
		f(&v, gf, p); //v = A pn
		log_squarenorm(create_log_prefix_cg(iter) + "v: ", v);

		scalar_product(&rho, p, v);
		divide(&alpha, omega, rho);
		multiply(&tmp1, minus_one, alpha); //alpha = (rn, rn)/(pn, Apn) --> alpha = omega/rho

		saxpy(x, tmp1, p, *x); //xn+1 = xn + alpha*p = xn - tmp1*p = xn - (-tmp1)*p
		log_squarenorm(create_log_prefix_cg(iter) + "x: ", *x);

		//rn+1 = rn - alpha*v -> rhat
		//NOTE: for beta one needs a complex number at the moment, therefore, this is done with "rho_next" instead of "resid"
		if(params.get_use_merge_kernels_spinor()) {
			physics::lattices::saxpy_AND_squarenorm(&rn, alpha, v, rn, rho_next);
			log_squarenorm(create_log_prefix_cg(iter) + "rn: ", rn);
		} else {
			saxpy(&rn, alpha, v, rn);
			scalar_product(&rho_next, rn, rn);
			log_squarenorm(create_log_prefix_cg(iter) + "rn: ", rn);
		}
		
		if(iter % RESID_CHECK_FREQUENCY == 0) {
			resid = rho_next.get().re;
			//if(USE_ASYNC_COPY) {
			//  if(iter) {
			//    resid_event.wait();
			//    resid = resid_rho.re;
			//  } else {
			//    // first iteration
			//    resid = prec;
			//  }
			//  resid_event = clmem_rho_next.dump_async(&resid_rho);
			//} else {
			//  clmem_rho_next.dump(&resid_rho);
			//  resid = resid_rho.re;
			//  //this is the orig. call
			//  //set_float_to_global_squarenorm_device(&clmem_rn, clmem_resid);
			//  //get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));
			//}

			logger.debug() << create_log_prefix_cg(iter) << "resid: " << resid;
			testIfResiduumIsNan(resid, iter);
			
			if(resid < prec && iter >= MINIMUM_ITERATIONS) {
				//if(USE_ASYNC_COPY) {
				//  // make sure everything using our event is completed
				//  resid_event.wait();
				//}
			  logger.debug() << create_log_prefix_cg(iter) << "Solver converged in " << iter << " iterations! resid:\t" << resid;

				// report on performance
				if(logger.beInfo()) {
					// we are always synchroneous here, as we had to recieve the residium from the device
					const uint64_t duration = timer.getTime();
					const uint64_t duration_noWarmup = timer_noWarmup.getTime();

					// calculate flops
					const unsigned refreshs = iter / params.get_iter_refresh() + 1;
					const cl_ulong mf_flops = f.get_flops();
					logger.trace() << "mf_flops: " << mf_flops;

					cl_ulong flops_per_iter = mf_flops + 2 * get_flops<Spinorfield_eo, scalar_product>(system)
					                          + 2 * ::get_flops<hmc_complex, complexdivide>() + 2 * ::get_flops<hmc_complex, complexmult>()
					                          + 3 * get_flops<Spinorfield_eo, saxpy>(system);
					cl_ulong flops_per_refresh = mf_flops + get_flops<Spinorfield_eo, saxpy>(system) + get_flops<Spinorfield_eo, scalar_product>(system);
					cl_ulong total_flops = iter * flops_per_iter + refreshs * flops_per_refresh;
					cl_ulong noWarmup_flops = (iter - 1) * flops_per_iter + (refreshs - 1) * flops_per_refresh;

					reportPerformance_cg(iter, duration, duration_noWarmup, total_flops, noWarmup_flops);
				}
				log_squarenorm(create_log_prefix_cg(iter) + "x (final): ", *x);
				return iter;
			}

			if(iter == 0) {
				timer_noWarmup.reset();
			}
		}

		divide(&beta, rho_next, omega); //beta = (rn+1, rn+1)/(rn, rn) --> alpha = rho_next/omega

		multiply(&tmp2, minus_one, beta);
		saxpy(&p, tmp2, p, rn); //pn+1 = rn+1 + beta*pn
		log_squarenorm(create_log_prefix_cg(iter) + "p: ", p);
	}

	logger.fatal() << create_log_prefix_cg(iter) << "Solver did not solve in " << params.get_cgmax() << " iterations. Last resid: " << resid;
	throw physics::algorithms::solvers::SolverDidNotSolve(iter, __FILE__, __LINE__);
}

int cg_multidev(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, const hmc_float prec)
{
	using namespace physics::lattices;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	const auto & params = system.get_inputparameters();
	const int MINIMUM_ITERATIONS = params.get_cg_minimum_iteration_count();
	if(MINIMUM_ITERATIONS) {
		logger.warn() << "Minimum iterations set to " << MINIMUM_ITERATIONS << " -- should be used *only* for inverter benchmarking!";
	}

	/// @todo start timer synchronized with device(s)
	klepsydra::Monotonic timer;
	klepsydra::Monotonic timer_noWarmup;

	const Spinorfield_eo p(system);
	const Spinorfield_eo rn(system);
	const Spinorfield_eo v(system);

	const Scalar<hmc_float> tmp_float(system);
	const Scalar<hmc_complex> tmp_complex(system);
	hmc_complex alpha;
	hmc_complex beta;
	hmc_complex omega;
	hmc_complex rho;
	hmc_complex rho_next{std::nan(""), std::nan("")};
	hmc_complex tmp1;
	hmc_complex tmp2;

	hmc_float resid;
	int iter = 0;

	// report source and initial solution
	log_squarenorm(create_log_prefix_cg(iter) + "b (initial): ", b);
	log_squarenorm(create_log_prefix_cg(iter) + "x (initial): ", *x);

	//NOTE: here, most of the complex numbers may also be just hmc_floats. However, for this one would need some add. functions...
	for(iter = 0; iter < params.get_cgmax(); iter ++) {
		if(iter % params.get_iter_refresh() == 0) {
			//rn = A*inout
			f(&rn, gf, *x);
			log_squarenorm(create_log_prefix_cg(iter) + "rn: ", rn);
			//rn = source - A*inout
			saxpy(&rn, hmc_complex_one, rn, b);
			log_squarenorm(create_log_prefix_cg(iter) + "rn: ", rn);
			//p = rn
			copyData(&p, rn);
			log_squarenorm(create_log_prefix_cg(iter) + "p: ", p);
			//omega = (rn,rn)
			omega = hmc_complex{squarenorm(rn, &tmp_float), 0};
		} else {
			//update omega
			omega = rho_next;
		}
		//v = A pn
		f(&v, gf, p);
		log_squarenorm(create_log_prefix_cg(iter) + "v: ", v);


		//alpha = (rn, rn)/(pn, Apn) --> alpha = omega/rho
		rho = scalar_product(p, v, &tmp_complex);
		alpha = complexdivide(omega, rho);
		tmp1 = complexmult(hmc_complex_minusone, alpha);

		//xn+1 = xn + alpha*p = xn - tmp1*p = xn - (-tmp1)*p
		saxpy(x, tmp1, p, *x);
		log_squarenorm(create_log_prefix_cg(iter) + "x: ", *x);

		//switch between original version and kernel merged one
		if(params.get_use_merge_kernels_spinor()) {
			//merge two calls:
			//rn+1 = rn - alpha*v -> rhat
			//and
			//rho_next = |rhat|^2
			//rho_next is a complex number, set its imag to zero
			//spinor_code->saxpy_AND_squarenorm_eo_device(&clmem_v_eo, &clmem_rn_eo, &clmem_alpha, &clmem_rn_eo, &clmem_rho_next);
			throw Print_Error_Message("Kernel merging currently not implemented for multidev", __FILE__, __LINE__);
			log_squarenorm(create_log_prefix_cg(iter) + "rn: ", rn);
		} else {
			//rn+1 = rn - alpha*v -> rhat
			saxpy(&rn, alpha, v, rn);
			log_squarenorm(create_log_prefix_cg(iter) + "rn: ", rn);

			//calc residuum
			//NOTE: for beta one needs a complex number at the moment, therefore, this is done with "rho_next" instead of "resid"
			resid = squarenorm(rn, &tmp_float);
			rho_next = hmc_complex{resid, 0.};
		}
		logger.debug() << create_log_prefix_cg(iter) << "resid: " << resid;
		//test if resid is NAN
		if(resid != resid) {
			logger.fatal() << create_log_prefix_cg(iter) << "NAN occured!";
			throw SolverStuck(iter, __FILE__, __LINE__);
		}
		if(resid < prec && iter >= MINIMUM_ITERATIONS) {
			logger.debug() << create_log_prefix_cg(iter) << "Solver converged in " << iter << " iterations! resid:\t" << resid;

			// report on performance
			if(logger.beInfo()) {
				// we are always synchroneous here, as we had to recieve the residium from the device
				const uint64_t duration = timer.getTime();
				const uint64_t duration_noWarmup = timer_noWarmup.getTime();

				// calculate flops
				const unsigned refreshs = iter / params.get_iter_refresh() + 1;
				const cl_ulong mf_flops = f.get_flops();
				logger.trace() << "mf_flops: " << mf_flops;

				cl_ulong flops_per_iter = mf_flops + 2 * get_flops<Spinorfield_eo, scalar_product>(system)
				                          + 2 * ::get_flops<hmc_complex, complexdivide>() + 2 * ::get_flops<hmc_complex, complexmult>()
				                          + 3 * get_flops<Spinorfield_eo, saxpy>(system);
				cl_ulong flops_per_refresh = mf_flops + get_flops<Spinorfield_eo, saxpy>(system) + get_flops<Spinorfield_eo, scalar_product>(system);
				cl_ulong total_flops = iter * flops_per_iter + refreshs * flops_per_refresh;
				cl_ulong noWarmup_flops = (iter - 1) * flops_per_iter + (refreshs - 1) * flops_per_refresh;

				// report performance
				logger.info() << create_log_prefix_cg(iter) << "CG completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f) << " Gflops. Performed " << iter << " iterations. Performance after warmup: " << (noWarmup_flops / duration_noWarmup / 1000.f) << " Gflops.";
			}
			// report on solution
			log_squarenorm(create_log_prefix_cg(iter) + "x (final): ", *x);
			return iter;
		}

		//beta = (rn+1, rn+1)/(rn, rn) --> alpha = rho_next/omega
		beta = complexdivide(rho_next, omega);

		//pn+1 = rn+1 + beta*pn
		tmp2 = complexmult(hmc_complex_minusone, beta);
		saxpy(&p, tmp2, p, rn);
		log_squarenorm(create_log_prefix_cg(iter) + "p: ", p);

		if(iter == 0) {
			timer_noWarmup.reset();
		}

	}

	logger.fatal() << create_log_prefix_cg(iter) << "Solver did not solve in " << params.get_cgmax() << " iterations. Last resid: " << resid;
	throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}

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

static std::string create_log_prefix_cg(int number) noexcept
{
  return create_log_prefix_solver("CG", number);
}

static std::string create_log_prefix_bicgstab(int number) noexcept
{
  return create_log_prefix_solver("BICGSTAB", number);
}
