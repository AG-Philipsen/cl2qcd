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
#include <cmath>
#include <sstream>

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

/**
 * A "save" version of the bicgstab algorithm.
 *
 * This is chosen if "bicgstab_save" is selected in the inputfile
 */
static int bicgstab_save(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, hmc_float prec);

/**
 * BICGstab implementation with a different structure than "save" one, similar to tmlqcd. This should be the default bicgstab.
 * In particular this version does not perform the check if the "real" residuum is sufficiently small!
 */
static int bicgstab_fast(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, hmc_float prec);

physics::algorithms::solvers::SolverStuck::SolverStuck(int iterations, std::string filename, int linenumber) : SolverException(create_solver_stuck_message(iterations), iterations, filename, linenumber) { };

static void trace_squarenorm(const std::string msg, const physics::lattices::Spinorfield_eo& x)
{
	if(logger.beTrace()) {
		hmc_float tmp = squarenorm(x);
		logger.trace() << msg << std::scientific << std::setprecision(10) << tmp;
	}
}

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

static int bicgstab_save(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, const hmc_float prec)
{
	// TODO this function often contains -1 in the comment but 1 in the code...
	using physics::lattices::Spinorfield;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	auto params = system.get_inputparameters();
	hmc_float resid;
	hmc_complex alpha, omega, rho;

	const Spinorfield v(system);
	const Spinorfield p(system);
	const Spinorfield rn(system);
	const Spinorfield rhat(system);
	const Spinorfield s(system);
	const Spinorfield t(system);
	const Spinorfield aux(system);

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
		if(std::abs(rho_next.re) < 1e-25 && std::abs(rho_next.im) < 1e-25 ) {
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
			saxpy(&aux, {1., 0.}, aux, b);
			//trueresid = (aux, aux)
			hmc_float trueresid = squarenorm(aux);
			logger.debug() << "\tsolver converged! true resid:\t" << trueresid;
			if(trueresid < prec)
				return iter;
		}
	}
	throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}

static int bicgstab_fast(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, const hmc_float prec)
{
	using physics::lattices::Spinorfield;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	auto params = system.get_inputparameters();
	hmc_float resid;
	hmc_complex rho;

	const Spinorfield p(system);
	const Spinorfield rn(system);
	const Spinorfield rhat(system);
	const Spinorfield v(system);
	const Spinorfield s(system);
	const Spinorfield t(system);

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
		if(std::abs(rho_next.re) < 1e-25 && std::abs(rho_next.im) < 1e-25 ) {
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

int physics::algorithms::solvers::cg(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, const hmc_float prec)
{
	using physics::lattices::Spinorfield;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	auto params = system.get_inputparameters();

	const Spinorfield rn(system);
	const Spinorfield p(system);
	const Spinorfield v(system);

	hmc_complex rho_next;

	//CP: here I do not use clmem_rnhat anymore and saved one scalar_product (omega)
	//NOTE: here, most of the complex numbers may also be just hmc_floats. However, for this one would need some add. functions...
	int iter;
	for(iter = 0; iter < params.get_cgmax(); iter ++) {
		hmc_complex omega;
		if(iter % params.get_iter_refresh() == 0) {
			//rn = A*inout
			f(&rn, gf, *x);
			//rn = source - A*inout
			saxpy(&rn, {1., 0.}, b, rn);
			//p = rn
			copyData(&p, rn);
			//omega = (rn,rn)
			omega = scalar_product(rn, rn);
		} else {
			//update omega
			omega = rho_next;
		}
		//v = A pn
		f(&v, gf, p);
		//alpha = (rn, rn)/(pn, Apn) --> alpha = omega/rho
		hmc_complex rho = scalar_product(p, v);
		hmc_complex alpha = complexdivide(omega, rho);
		hmc_complex tmp1 = complexsubtract( {0., 0.}, alpha);

		//xn+1 = xn + alpha*p = xn - tmp1*p = xn - (-tmp1)*p
		saxpy(x, tmp1, p, *x);
		//rn+1 = rn - alpha*v -> rhat
		saxpy(&rn, alpha, v, rn);

		//calc residuum
		//NOTE: for beta one needs a complex number at the moment, therefore, this is done with "rho_next" instead of "resid"
		rho_next = scalar_product(rn, rn);
		hmc_float resid = rho_next.re;
		//this is the orig. call
		//set_float_to_global_squarenorm_device(clmem_rn, clmem_resid);
		//get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));

		logger.debug() << "resid: " << resid;
		//test if resid is NAN
		if(resid != resid) {
			logger.fatal() << "\tNAN occured in cg!";
			throw SolverStuck(iter, __FILE__, __LINE__);
		}
		if(resid < prec)
			return iter;

		//beta = (rn+1, rn+1)/(rn, rn) --> alpha = rho_next/omega
		hmc_complex beta = complexdivide(rho_next, omega);

		//pn+1 = rn+1 + beta*pn
		hmc_complex tmp2 = complexsubtract( {0., 0.}, beta);
		saxpy(&p, tmp2, p, rn);
	}
	throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}

int physics::algorithms::solvers::bicgstab(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, hmc_float prec)
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

static int bicgstab_save(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, const hmc_float prec)
{
	using physics::lattices::Spinorfield_eo;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	// TODO start timer synchronized with device(s)
	klepsydra::Monotonic timer;

	auto params = system.get_inputparameters();

	const Spinorfield_eo s(system);
	const Spinorfield_eo t(system);
	const Spinorfield_eo v(system);
	const Spinorfield_eo p(system);
	const Spinorfield_eo rn(system);
	const Spinorfield_eo rhat(system);
	const Spinorfield_eo aux(system);

	//CP: these have to be on the host
	unsigned retests = 0;

	int cgmax = params.get_cgmax();

	hmc_complex alpha, omega, rho, rho_next;
	hmc_float resid;

	int iter;
	for(iter = 0; iter < cgmax; iter++) {
		if(iter % params.get_iter_refresh() == 0) {
			v.zero();
			p.zero();

			f(&rn, gf, *x);
			saxpy(&rn, {1., 0.}, rn, b);

			copyData(&rhat, rn);

			alpha = {1., 0.};
			omega = {1., 0.};
			rho = {1., 0.};
		}
		rho_next = scalar_product(rhat, rn);
		//check if algorithm is stuck
		if(std::abs(rho_next.re) < 1e-25 && std::abs(rho_next.im) < 1e-25 ) {
			//print the last residuum
			logger.fatal() << "\t\t\tsolver stuck at resid:\t" << resid;
			throw SolverStuck(iter, __FILE__, __LINE__);
		}
		hmc_complex tmp1 = complexdivide(rho_next, rho);
		rho = rho_next;
		hmc_complex tmp2 = complexdivide(alpha, omega);
		hmc_complex beta = complexmult(tmp1, tmp2);

		tmp1 = complexmult(beta, omega);
		tmp2 = complexsubtract( {0., 0.}, tmp1);
		saxsbypz(&p, beta, p, tmp2, v, rn);

		f(&v, gf, p);

		tmp1 = scalar_product(rhat, v);
		alpha = complexdivide(rho, tmp1);

		saxpy(&s, alpha, v, rn);

		f(&t, gf, s);

		tmp1 = scalar_product(t, s);
		//!!CP: can this also be global_squarenorm??
		tmp2 = scalar_product(t, t);
		omega = complexdivide(tmp1, tmp2);

		saxpy(&rn, omega, t, s);

		saxsbypz(x, alpha, p, omega, s, *x);

		resid = squarenorm(rn);

		logger.debug() << "resid: " << resid;
		//test if resid is NAN
		if(resid != resid) {
			logger.fatal() << "\tNAN occured in bicgstab_eo!";
			throw SolverStuck(iter, __FILE__, __LINE__);
		}
		if(resid < prec) {
			++retests;

			f(&aux, gf, *x);
			saxpy(&aux, {1., 0.}, aux, b);

			hmc_float trueresid = squarenorm(aux);
			if(trueresid < prec) {
				// report on performance
				if(logger.beInfo()) {
					// we are always synchroneous here, as we had to recieve the residium from the device
					uint64_t duration = timer.getTime();

					// calculate flops
					unsigned refreshs = iter / params.get_iter_refresh() + 1;
					cl_ulong mf_flops = f.get_flops();

					// TODO fix flop calculation
					//cl_ulong total_flops = 4 * get_flop_size("scalar_product_eoprec") + 4 * get_flop_size("ratio") + 3 * get_flop_size("product") + 2 * spinor_code->get_flop_size("saxsbypz_eoprec") + 2 * mf_flops + 2 * spinor_code->get_flop_size("saxpy_eoprec") + get_flop_size("global_squarenorm_eoprec");
					//total_flops *= iter;

					//total_flops += refreshs * (mf_flops + spinor_code->get_flop_size("saxpy_eoprec"));

					//total_flops += retests * (mf_flops + spinor_code->get_flop_size("saxpy_eoprec") + get_flop_size("global_squarenorm_eoprec"));
					cl_ulong total_flops = 0;

					// report performanc
					logger.info() << "BiCGstab_save completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f) << " Gflops. Performed " << iter << " iterations";
				}

				// we are done here
				return iter;
			}
		}
	}
	throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}

static int bicgstab_fast(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, hmc_float prec)
{
	using physics::lattices::Spinorfield_eo;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	// TODO start timer synchronized with device(s)
	klepsydra::Monotonic timer;

	auto params = system.get_inputparameters();

	const Spinorfield_eo p(system);
	const Spinorfield_eo rn(system);
	const Spinorfield_eo rhat(system);
	const Spinorfield_eo s(system);
	const Spinorfield_eo t(system);
	const Spinorfield_eo v(system);

	hmc_complex rho;

	int iter;
	for(iter = 0; iter < params.get_cgmax(); iter++) {
		if(iter % params.get_iter_refresh() == 0) {
			//initial r_n, saved in p
			f(&rn, gf, *x);
			saxpy(&p, {1., 0}, rn, b);
			//rhat = p
			copyData(&rhat, p);
			//r_n = p
			copyData(&rn, p);
			//rho = (rhat, rn)
			rho = scalar_product(rhat, rn);
		}
		//resid = (rn,rn)
		hmc_float resid = squarenorm(rn);

		logger.debug() << "resid: " << resid;
		//test if resid is NAN
		if(resid != resid) {
			logger.fatal() << "\tNAN occured in bicgstab_eo!";
			throw SolverStuck(iter, __FILE__, __LINE__);
		}
		if(resid < prec) {
			// report on performance
			if(logger.beInfo()) {
				// we are always synchroneous here, as we had to recieve the residium from the device
				uint64_t duration = timer.getTime();

				// calculate flops
				unsigned refreshs = iter / params.get_iter_refresh() + 1;
				cl_ulong mf_flops = f.get_flops();

				// TODO fix flop calculation
				//cl_ulong total_flops = get_flop_size("global_squarenorm_eoprec") + 2 * mf_flops + 4 * get_flop_size("scalar_product_eoprec") + 4 * get_flop_size("ratio") + 2 * spinor_code->get_flop_size("saxpy_eoprec") + 2 * spinor_code->get_flop_size("saxsbypz_eoprec") + 3 * get_flop_size("product");
				//total_flops *= iter;

				//total_flops += refreshs * (mf_flops + spinor_code->get_flop_size("saxpy_eoprec") + get_flop_size("scalar_product_eoprec"));
				cl_ulong total_flops = 0;

				// report performanc
				logger.info() << "BiCGstab completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f) << " Gflops. Performed " << iter << " iterations";
			}

			// we are done here
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
		if(std::abs(rho_next.re) < 1e-25 && std::abs(rho_next.im) < 1e-25 ) {
			//print the last residuum
			logger.fatal() << "\t\t\tsolver stuck at resid:\t" << resid;
			SolverStuck(iter, __FILE__, __LINE__);
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

int physics::algorithms::solvers::cg(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, const hmc_float prec)
{
	using physics::lattices::Spinorfield_eo;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	// TODO start timer synchronized with device(s)
	klepsydra::Monotonic timer;

	auto params = system.get_inputparameters();

	/// @todo make configurable from outside
	const int RESID_CHECK_FREQUENCY = params.get_cg_iteration_block_size();
	if(RESID_CHECK_FREQUENCY != 1) {
		logger.warn() << "The residuum is currently checked at every iteration, even though you selected a check only every " << RESID_CHECK_FREQUENCY << " iterations.";
	}
	const bool USE_ASYNC_COPY = params.get_cg_use_async_copy();
	if(USE_ASYNC_COPY) {
		logger.warn() << "Asynchroneous copying in the CG is currently unimplemented!";
	}

	hmc_complex rho_next;

	const Spinorfield_eo p(system);
	const Spinorfield_eo rn(system);
	const Spinorfield_eo v(system);

	trace_squarenorm("CG: b: ", b);
	trace_squarenorm("CG: x: ", *x);

	//this corresponds to the above function
	//NOTE: here, most of the complex numbers may also be just hmc_floats. However, for this one would need some add. functions...
	int iter;
	for(iter = 0; iter < params.get_cgmax(); iter ++) {
		hmc_complex omega;
		if(iter % params.get_iter_refresh() == 0) {
			//rn = A*inout
			f(&rn, gf, *x);
			trace_squarenorm("CG: rn: ", rn);
			//rn = source - A*inout
			saxpy(&rn, {1., 0.}, rn, b);
			trace_squarenorm("CG: rn: ", rn);
			//p = rn
			copyData(&p, rn);
			trace_squarenorm("CG: p: ", p);
			//omega = (rn,rn)
			omega = scalar_product(rn, rn);
		} else {
			//update omega
			omega = rho_next;
		}
		//v = A pn
		f(&v, gf, p);
		trace_squarenorm("CG: v: ", v);

		//alpha = (rn, rn)/(pn, Apn) --> alpha = omega/rho
		hmc_complex rho = scalar_product(p, v);
		hmc_complex alpha = complexdivide(omega, rho);
		hmc_complex tmp1 = complexsubtract( {0., 0.}, alpha);

		//xn+1 = xn + alpha*p = xn - tmp1*p = xn - (-tmp1)*p
		saxpy(x, tmp1, p, *x);
		trace_squarenorm("CG: x: ", *x);
		//switch between original version and kernel merged one
		if(params.get_use_merge_kernels_spinor()) {
			//merge two calls:
			//rn+1 = rn - alpha*v -> rhat
			//and
			//rho_next = |rhat|^2
			//rho_next is a complex number, set its imag to zero
			//spinor_code->saxpy_AND_squarenorm_eo_device(&clmem_v_eo, &clmem_rn_eo, &clmem_alpha, &clmem_rn_eo, &clmem_rho_next);
			throw Print_Error_Message("Kernel merging currently not implemented", __FILE__, __LINE__);
		} else {
			//rn+1 = rn - alpha*v -> rhat
			saxpy(&rn, alpha, v, rn);
			trace_squarenorm("CG: rn: ", rn);

			//calc residuum
			//NOTE: for beta one needs a complex number at the moment, therefore, this is done with "rho_next" instead of "resid"
			rho_next = scalar_product(rn, rn);
		}
		//if(iter % RESID_CHECK_FREQUENCY == 0) {
		hmc_float resid = rho_next.re;
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

		logger.debug() << "resid: " << resid;
		//test if resid is NAN
		if(resid != resid) {
			logger.fatal() << "\tNAN occured in cg_eo!";
			throw SolverStuck(iter, __FILE__, __LINE__);
		}
		if(resid < prec) {
			//if(USE_ASYNC_COPY) {
			//  // make sure everything using our event is completed
			//  resid_event.wait();
			//}
			// report on performance
			if(logger.beInfo()) {
				// we are always synchroneous here, as we had to recieve the residium from the device
				uint64_t duration = timer.getTime();

				// calculate flops
				unsigned refreshs = iter / params.get_iter_refresh() + 1;
				cl_ulong mf_flops = f.get_flops();

				// TODO reenable flop calculation
				//cl_ulong total_flops = mf_flops + 3 * get_flop_size("scalar_product_eoprec") + 2 * get_flop_size("ratio") + 2 * get_flop_size("product") + 3 * spinor_code->get_flop_size("saxpy_eoprec");
				//total_flops *= iter;

				//total_flops += refreshs * (mf_flops + spinor_code->get_flop_size("saxpy_eoprec") + get_flop_size("scalar_product_eoprec"));
				cl_ulong total_flops = 0;

				// report performanc
				logger.info() << "CG completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f) << " Gflops. Performed " << iter << " iterations";
			}

			return iter;
		}
		//}

		//beta = (rn+1, rn+1)/(rn, rn) --> alpha = rho_next/omega
		hmc_complex beta = complexdivide(rho_next, omega);

		//pn+1 = rn+1 + beta*pn
		hmc_complex tmp2 = complexsubtract( {0., 0.}, beta);
		saxpy(&p, tmp2, rn, p);
		trace_squarenorm("CG: p: ", p);
	}
	throw SolverDidNotSolve(iter, __FILE__, __LINE__);
}
