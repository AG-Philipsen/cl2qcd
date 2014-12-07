/** @file
 * Implementation of the cg algorithm
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

#include "cg.hpp"

#include "../../lattices/scalar_complex.hpp"
#include "../../lattices/util.hpp"

static std::string create_log_prefix_cg(int number) noexcept;
//@todo: move to own file
static std::string create_log_prefix_solver(std::string name, int number) noexcept;

namespace {

int cg_singledev(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, const hmc_float prec);
int cg_multidev(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& f, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, const hmc_float prec);

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
