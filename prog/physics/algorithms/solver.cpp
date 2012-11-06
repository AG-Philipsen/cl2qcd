/** @file
 * Implementation of the physics::algorithms::solver class
 *
 * (c) 2012 Christopher Pinke <pinke@th.uni-frankfurt.de>
 */

#include "solver.hpp"

#include "../../hardware/system.hpp"
#include "../../logger.hpp"

#include "solver.hpp"

int physics::algorithms::solver::Solver::get_iter_max() const noexcept
{
	return iter_max;
}

int physics::algorithms::solver::Solver::get_iter_refresh() const noexcept
{
	return iter_refresh;
}

int physics::algorithms::solver::Solver::get_iter()
{
	return iter;
}

hmc_float physics::algorithms::solver::Solver::get_prec() const noexcept
{
	return acc_prec;
}

bool physics::algorithms::solver::Cg_eo::solve() const
{
	//NOTE: here, most of the complex numbers may also be just hmc_floats. However, for this one would need some add. functions...
	klepsydra::Monotonic timer;
	if(logger.beInfo()) {
		cl_event start_event;
		spinor_code->get_device()->enqueue_marker(&start_event);
		clSetEventCallback(start_event, CL_COMPLETE, resetTimerOnComplete, &timer);
		clReleaseEvent(start_event);
	}
	int iter = 0;
	for(iter = 0; iter < iter_max; iter ++) {
		if(iter % iter_refresh == 0) {
			//rn = A*inout
			f(x, rn, gf);
			//rn = source - A*inout
			spinor_code->saxpy_eoprec_device(rn, b, one, rn);
			//p = rn
			hardware::buffers::copyData(p, rn);
			//omega = (rn,rn)
			spinor_code->set_complex_to_scalar_product_eoprec_device(rn, rn, omega);
		} else {
			//update omega
			hardware::buffers::copyData(omega, rho_next);
		}
		//v = A pn
		f(p, v, get_gf());

		//alpha = (rn, rn)/(pn, Apn) --> alpha = omega/rho
		spinor_code->set_complex_to_scalar_product_eoprec_device(p, v, rho);
		spinor_code->set_complex_to_ratio_device(omega, rho, alpha);
		spinor_code->set_complex_to_product_device(alpha, minusone, tmp1);

		//xn+1 = xn + alpha*p = xn - tmp1*p = xn - (-tmp1)*p
		spinor_code->saxpy_eoprec_device(p, x, tmp1, x);
		//switch between original version and kernel merged one
		if(merge_kernels == false) {
			//rn+1 = rn - alpha*v -> rhat
			spinor_code->saxpy_eoprec_device(v, rn, alpha, rn);

			//calc residuum
			//NOTE: for beta one needs a complex number at the moment, therefore, this is done with "rho_next" instead of "resid"
			spinor_code->set_complex_to_scalar_product_eoprec_device(rn, rn, rho_next);
		} else {
			//merge two calls:
			//rn+1 = rn - alpha*v -> rhat
			//and
			//rho_next = |rhat|^2
			//rho_next is a complex number, set its imag to zero
			spinor_code->saxpy_AND_squarenorm_eo_device(v, rn, alpha, rn, rho_next);
		}
		hmc_complex tmp;
		rho_next->dump(&tmp);
		hmc_float resid = tmp.re;
		//this is the orig. call
		//set_float_to_global_squarenorm_device(&clmem_rn, clmem_resid);
		//get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));

		logger.debug() << "resid: " << resid;
		//test if resid is NAN
		if(resid != resid) {
			logger.fatal() << "\tNAN occured in cg_eo!";
			return -iter;
		}
		if(resid < acc_prec) {
			// report on performance
			if(logger.beInfo()) {
				// we are always synchroneous here, as we had to recieve the residium from the device
				uint64_t duration = timer.getTime();

				// calculate flops
				unsigned refreshs = iter / iter_refresh + 1;
				cl_ulong mf_flops = f.get_Flops();

				cl_ulong total_flops = mf_flops + 3 * spinor_code->get_flop_size("scalar_product_eoprec") + 2 * spinor_code->get_flop_size("ratio") + 2 * spinor_code->get_flop_size("product") + 3 * spinor_code->get_flop_size("saxpy_eoprec");
				total_flops *= iter;

				total_flops += refreshs * (mf_flops + spinor_code->get_flop_size("saxpy_eoprec") + spinor_code->get_flop_size("scalar_product_eoprec"));

				// report performanc
				logger.info() << "CG completed in " << duration / 1000 << " ms @ " << (total_flops / duration / 1000.f) << " Gflops. Performed " << iter << " iterations";
			}

			return iter;
		}

		//beta = (rn+1, rn+1)/(rn, rn) --> alpha = rho_next/omega
		spinor_code->set_complex_to_ratio_device(rho_next, omega, beta);

		//pn+1 = rn+1 + beta*pn
		spinor_code->set_complex_to_product_device(beta, minusone, tmp2);
		spinor_code->saxpy_eoprec_device(p, rn, tmp2, p);
	}
	return -1;
}
