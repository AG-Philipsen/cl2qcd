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

#include "solver.hpp"
#include "solver_shifted.hpp"

#include "../../host_functionality/logger.hpp"
//#include "../../operations_complex.h"
#include "../../meta/type_ops.hpp"
#include "../../meta/util.hpp"
#include "../lattices/scalar_complex.hpp"
#include "../lattices/algebra_real.hpp"
#include "../lattices/staggeredfield_eo.hpp"
//#include <cmath>
#include <sstream>
#include <vector>
#include <numeric>

static std::string create_log_prefix_cgm(int number) noexcept;
static void log_squarenorm_aux(const std::string& msg, const std::vector<physics::lattices::Staggeredfield_eo *> x, const int n) noexcept;
static void compare_sqnorm(const std::vector<hmc_float> a, const std::vector<physics::lattices::Staggeredfield_eo *> x) noexcept;

int physics::algorithms::solvers::cg_m(const std::vector<physics::lattices::Staggeredfield_eo *> x, const std::vector<hmc_float> sigma, const physics::fermionmatrix::Fermionmatrix_stagg_eo& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Staggeredfield_eo& b, const hardware::System& system, hmc_float prec)
{
	using namespace physics::lattices;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;
	
	auto params = system.get_inputparameters();
  
	/// @todo start timer synchronized with device(s)
	klepsydra::Monotonic timer;
	klepsydra::Monotonic timer_noWarmup;
	/// @todo make configurable from outside
	const bool USE_ASYNC_COPY = params.get_cg_use_async_copy();
	const int MINIMUM_ITERATIONS = params.get_cg_minimum_iteration_count();
	if(USE_ASYNC_COPY) {
		logger.warn() << "Asynchroneous copying in the CG-M is currently unimplemented!";
	}
	if(MINIMUM_ITERATIONS) {
		logger.warn() << "Minimum iterations set to " << MINIMUM_ITERATIONS << " -- should be used *only* for CGM benchmarking!";
	}
	
	if(squarenorm(b)==0){
	  for(uint i=0; i<sigma.size(); i++)
		x[i]->set_zero();   
	  return 0;
	}
	
	std::vector<hmc_float> xsq(x.size(), 0);
	
	if(sigma.size() != x.size())
		throw std::invalid_argument("Wrong size of multi-shifted inverter parameters!");
	
	//The number of values of constants sigma as well as the number of fields x
	const int Neqs = sigma.size();
	
	//Auxiliary staggered fields
	const Staggeredfield_eo r(system);
	const Staggeredfield_eo p(system);
	std::vector<Staggeredfield_eo*> ps;
	
	//Auxiliary scalar vectors
	std::vector<const Scalar<hmc_float>*> alpha;
	std::vector<const Scalar<hmc_float>*> beta;
	std::vector<const Scalar<hmc_float>*> zeta_i;   //This is zeta at the step iter-1
	std::vector<const Scalar<hmc_float>*> zeta_ii;  //This is zeta at the step iter
	std::vector<const Scalar<hmc_float>*> zeta_iii; //This is zeta at the step iter+1
	std::vector<const Scalar<hmc_float>*> shift;    //This is to store constants sigma
	std::vector<bool> single_system_converged;        //This is to stop calculation on single system
	std::vector<uint> single_system_iter;             //This is to calculate performance properly
	//Auxiliary scalars
	const Scalar<hmc_float> alpha_scalar_prev(system);   //This is alpha_scalar at the step iter-1
	const Scalar<hmc_float> alpha_scalar(system);        //This is alpha_scalar at the step iter
	const Scalar<hmc_float> beta_scalar_prev(system);    //This is beta_scalar at the step iter-1
	const Scalar<hmc_float> beta_scalar(system);         //This is beta_scalar at the step iter
	
	//Auxiliary containers for temporary saving
	const Staggeredfield_eo v(system); //this is to store A.p
	const Scalar<hmc_float> tmp1(system);                     //this is to store (r,r) before updating r
	const Scalar<hmc_float> tmp2(system);                     //this is to store (r,r) after updating r
	const Scalar<hmc_float> tmp3(system);                     //this is to store (p,v) as Scalar
	const Scalar<hmc_float> num(system);                      //this is to store constants numerators
	const Scalar<hmc_float> den(system);                      //this is to store constants denumerators
	bool converged=false;                                       //this is to start to check the residuum

	//Auxiliary constants as Scalar
	const Scalar<hmc_float> zero(system);
	zero.store(0.0);
	const Scalar<hmc_float> one(system);
	one.store(1.0);
	
	
	hmc_float resid;
	int iter = 0;
	
	//Initialization auxilary and output quantities
	for(int i=0; i<Neqs; i++){
		x[i]->set_zero();                                     // x[i] = 0 
		ps.push_back(new Staggeredfield_eo(system));    
		copyData(ps[i], b);                                    // ps[i] = b
		beta.push_back(new Scalar<hmc_float>(system));
		alpha.push_back(new Scalar<hmc_float>(system));
		alpha[i]->store(0.0);                    // alpha[i] = 0
		zeta_i.push_back(new Scalar<hmc_float>(system));
		zeta_ii.push_back(new Scalar<hmc_float>(system));
		zeta_iii.push_back(new Scalar<hmc_float>(system));
		zeta_i[i]->store(1.0);                    // zeta_i[i] = 1
		zeta_ii[i]->store(1.0);                   // zeta_ii[i] = 1
		shift.push_back(new Scalar<hmc_float>(system));
		shift[i]->store(sigma[i]);
		single_system_converged.push_back(false);             //no system converged
	}
	copyData(&r, b);                          // r = b
	copyData(&p, b);                          // p = b
	scalar_product_real_part(&tmp1, r, r);           // set tmp1 = (r, r) for the first iteration
	beta_scalar.store(1.0);      // beta_scalar = 1, here I should set beta_scalar_prev
						  // but in this way I can set beta_scalar_prev at the begin
						  // of the loop over iter recursively.
	alpha_scalar.store(0.0);    // alpha_scalar = 0. The same as beta_scalar above.
	
	//At first, in the log string I will report all the data about the system.
	//To avoid to print to shell tens of lines per time, since Neqs can be 
	//also 20 or something like that, reduce report_num.
	const int report_num = Neqs;
	if(report_num>Neqs)
		throw std::invalid_argument("In cg-m report_num cannot be bigger than Neqs!");
	log_squarenorm(create_log_prefix_cgm(iter) + "b (initial): ", b);
	log_squarenorm(create_log_prefix_cgm(iter) + "r (initial): ", r);
	log_squarenorm(create_log_prefix_cgm(iter) + "p (initial): ", p);
	log_squarenorm_aux(create_log_prefix_cgm(iter) + "x (initial)", x, report_num);
	log_squarenorm_aux(create_log_prefix_cgm(iter) + "ps (initial)", ps, report_num);
	
	for(iter = 0; iter < params.get_cgmax() || iter < MINIMUM_ITERATIONS; iter ++) {
		//Update beta_scalar: v=A.p and tmp1=(r,r) and tmp3=(p,v) ---> beta_scalar=(-1)*tmp1/tmp3
		copyData(&beta_scalar_prev,beta_scalar);  //before updating beta_scalar its value is saved
		A(&v, gf, p);
		log_squarenorm(create_log_prefix_cgm(iter) + "v: ", v);
		scalar_product_real_part(&tmp3, p, v); 
		divide(&beta_scalar, tmp1, tmp3);  //tmp1 is set from previous iteration
		subtract(&beta_scalar, zero, beta_scalar);
		//Update field r: r+=beta_scalar*A.p ---> r = r + beta_scalar*v
		saxpy(&r, beta_scalar, v, r);
		log_squarenorm(create_log_prefix_cgm(iter) + "r: ", r);
		//We store in tmp2 the quantity (r,r) that we use later. When we check
		//the residuum, then it is already calculated.
		scalar_product_real_part(&tmp2, r, r);
		if(logger.beDebug()){
			//Calculate squarenorm of the output field
			for(uint i=0; i<x.size(); i++)
				xsq[i] = squarenorm(*x[i]);
		}
		//Update alpha_scalar: alpha_scalar = tmp2/tmp1
		copyData(&alpha_scalar_prev,alpha_scalar); //before updating alpha_scalar its value is saved
		divide(&alpha_scalar, tmp2, tmp1);
		//Update field p: p = r + alpha_scalar*p
		saxpy(&p, alpha_scalar, p, r);
		//saxpy(p, alpha_scalar, *p, *r);
		log_squarenorm(create_log_prefix_cgm(iter) + "p: ", p);


		Vector<hmc_float> zeta_prev(Neqs, system);
		Vector<hmc_float> zeta(Neqs, system);
		Vector<hmc_float> zeta_foll(Neqs, system);
		Vector<hmc_float> masses(Neqs, system);
		std::vector<hmc_float> aux1, aux2, aux3;
		for(int k=0; k<Neqs; k++){
		  aux1.push_back(zeta_i[k]->get());
		  aux2.push_back(zeta_ii[k]->get());
		  aux3.push_back(zeta_iii[k]->get());
		}
		zeta_prev.store(aux1);
		zeta.store(aux2);
		zeta_foll.store(aux3);
		masses.store(sigma);
		
		update_zeta_cgm(&zeta_foll, zeta, zeta_prev, beta_scalar_prev, beta_scalar, alpha_scalar_prev, masses, Neqs);
		
		for(int k=0; k<Neqs; k++)
		  zeta_iii[k]->store((zeta_foll.get())[k]);
		
		//Loop over the system equations, namely over the set of sigma values
		for(int k=0; k<Neqs; k++){
			if(single_system_converged[k]==false){
				//Update zeta_iii[k]: num = zeta_i[k]*zeta_ii[k]*beta_scalar_prev
				//                    den = beta_scalar*alpha_scalar_prev*(zeta_i[k]-zeta_ii[k])+
				//                        + zeta_i[k]*beta_scalar_prev*(1-sigma[k]*beta_scalar)
				// ---> zeta_iii[k] = num/den
				//I calculate before the denominator to can use num as auxiliary variable
// 				subtract(&num, *zeta_i[k], *zeta_ii[k]);
// 				multiply(&num, num, alpha_scalar_prev);
// 				multiply(&num, num, beta_scalar);
// 				multiply(&den, *shift[k], beta_scalar);
// 				subtract(&den, one, den);
// 				multiply(&den, den, beta_scalar_prev);
// 				multiply(&den, den, *zeta_i[k]);
// 				add(&den, num, den);
// 				//Calculation of the numerator
// 				multiply(&num, *zeta_i[k], *zeta_ii[k]);
// 				multiply(&num, num, beta_scalar_prev);
// 				//Update of zeta_iii[k]
// 				divide(zeta_iii[k], num, den);
// 				logger.warn() << "zeta_iii[" << k << "] = " << zeta_iii[k]
				//Update beta[k]: beta[k] = beta_scalar*zeta_iii[k]/zeta_ii[k]
				multiply(beta[k], beta_scalar, *zeta_iii[k]);
				divide(beta[k], *beta[k], *zeta_ii[k]);
				//Update x[k]: x[k] = x[k] - beta[k]*ps[k]
				// ---> use num to store (- beta[k]) for saxpy
				subtract(&num, zero, *beta[k]);
				saxpy(x[k], num, *ps[k], *x[k]);
				//Update alpha[k]: num = alpha_scalar*zeta_iii[k]*beta[k]
				//                 den = zeta_ii[k]*beta_scalar
				// ---> alpha[k] = num/den
				multiply(&num, alpha_scalar, *zeta_iii[k]);
				multiply(&num, num, *beta[k]);
				multiply(&den, *zeta_ii[k], beta_scalar);
				divide(alpha[k], num, den);
				//Update ps[k]: ps[k] = zeta_iii[k]*r + alpha[k]*ps[k]
				saxpby(ps[k], *zeta_iii[k], r, *alpha[k], *ps[k]);
				//Check fields squarenorm for possible nan
				if(logger.beDebug()){
					if((squarenorm(*x[k]) != squarenorm(*x[k])) ||
					  (squarenorm(*ps[k]) != squarenorm(*ps[k]))){
						logger.fatal() << create_log_prefix_cgm(iter) << "NAN occured!";
						throw SolverStuck(iter, __FILE__, __LINE__);
					}
				}
				//Check if single system converged: ||zeta_iii[k] * r||^2 < prec
				// ---> v = zeta_iii[k] * r
				sax(&v, *zeta_iii[k], r);
				//Check if single system converged: (zeta_ii[k] * zeta_ii[k] * tmp1) < prec
				//multiply(&tmp3, *zeta_ii[k], *zeta_ii[k]);
				//multiply(&tmp3, tmp3, tmp1);
				if(squarenorm(v) < prec){
				//if(tmp3.get().re < prec){
					single_system_converged[k] = true;
					single_system_iter.push_back((uint)iter);
					logger.debug() << " ===> System number " << k << " converged after " << iter << " iterations! resid = " << tmp2.get();
				}
				//Adjust zeta for the following iteration
				copyData(zeta_i[k], zeta_ii[k]);
				copyData(zeta_ii[k], zeta_iii[k]);
			}
		}

		if(single_system_iter.size()==(uint)Neqs)
			converged = true;
		if(logger.beDebug()){
			compare_sqnorm(xsq, x);
		}
		log_squarenorm_aux(create_log_prefix_cgm(iter) + "x", x, report_num);
		log_squarenorm_aux(create_log_prefix_cgm(iter) + "ps", ps, report_num);
		//Check whether the algorithm converged
		if(converged) {
			//Calculate resid: 
			resid = tmp2.get();
			logger.debug() << create_log_prefix_cgm(iter) << "resid: " << resid;
			//test if resid is NAN
			if(resid != resid) {
				logger.fatal() << create_log_prefix_cgm(iter) << "NAN occured!";
				throw SolverStuck(iter, __FILE__, __LINE__);
			}
			if(resid < prec && iter >= MINIMUM_ITERATIONS) {
				logger.debug() << create_log_prefix_cgm(iter) << "Solver converged in " << iter << " iterations! resid:\t" << resid;
				// report on performance
				if(logger.beInfo()) {
					// we are always synchroneous here, as we had to recieve the residuum from the device
					const uint64_t duration = timer.getTime();
					const uint64_t duration_noWarmup = timer_noWarmup.getTime();
					// calculate flops
					const cl_ulong matrix_flops = A.get_flops();
					const int sum_of_partial_iter = std::accumulate(single_system_iter.begin(), single_system_iter.end(),0);
					logger.debug() << "matrix_flops: " << matrix_flops;
					
					cl_ulong flops_per_iter_no_inner_loop=(matrix_flops +
						2 * get_flops<Staggeredfield_eo, hmc_float, scalar_product_real_part>(system) +
						2 * 1 +
						    1 +
						2 * get_flops<Staggeredfield_eo, hmc_float, saxpy>(system));
					cl_ulong flops_per_iter_only_inner_loop=(
							 3 * 1 +
							11 * 1 +
							     1 +
							 3 * 1 +
							       get_flops<Staggeredfield_eo, hmc_float, saxpy>(system) +
							       get_flops<Staggeredfield_eo, hmc_float, saxpby>(system) +
							       get_flops<Staggeredfield_eo, hmc_float, sax>(system) +
							   5 * get_flops<Staggeredfield_eo, squarenorm>(system));
					
					cl_ulong total_flops = iter * flops_per_iter_no_inner_loop +
						sum_of_partial_iter * flops_per_iter_only_inner_loop;
					cl_ulong noWarmup_flops = (iter-1) * flops_per_iter_no_inner_loop +
									      flops_per_iter_only_inner_loop;
					
					logger.debug() << "total_flops: " << total_flops;
					// report performance
					logger.info() << create_log_prefix_cgm(iter) << "CG-M completed in " << std::setprecision(6) << duration / 1000.f << " ms @ " << ((hmc_float)total_flops / duration / 1000.f) << " Gflops. Performed " << iter << " iterations. Performance after warmup: " << ((hmc_float)noWarmup_flops / duration_noWarmup / 1000.f) << " Gflops.";
				}
				// report on solution
				log_squarenorm_aux(create_log_prefix_cgm(iter) + "x (final): ", x, report_num);
				
				//for(Staggeredfield_eo* i : x)
				//  logger.warn() << "sqnorm(x[i]) = " << squarenorm(*i);
				
				//Before returning I have to clean all the memory!!!
				meta::free_container(ps);
				meta::free_container(alpha);
				meta::free_container(beta);
				meta::free_container(zeta_i);
				meta::free_container(zeta_ii);
				meta::free_container(zeta_iii);
				meta::free_container(shift);
				
				return iter;
			}
		}
		//Set tmp1 from tmp2 for following iteration
		copyData(&tmp1, tmp2);
		if(iter == 0) {
			timer_noWarmup.reset();
		}
	//logger.info() << "cgm not properly executed in " << timer.getTime() / 1000.f << " ms";
	//return 0;
	}
	
	logger.fatal() << create_log_prefix_cgm(iter) << "Solver did not solve in " << params.get_cgmax() << " iterations. Last resid: " << tmp2.get();
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

static std::string create_log_prefix_cgm(int number) noexcept
{
  return create_log_prefix_solver("CG-M", number);
}

static void log_squarenorm_aux(const std::string& msg, const std::vector<physics::lattices::Staggeredfield_eo *> x, const int n) noexcept
{
	if(logger.beDebug()) {
		hmc_float tmp;
		for(int i=0; i<n; i++){
			tmp = squarenorm(*x[i]);
			logger.debug() << msg << "[field_" << i << "]: " << std::scientific << std::setprecision(16) << tmp;
		}
	}
}

static void compare_sqnorm(const std::vector<hmc_float> a, const std::vector<physics::lattices::Staggeredfield_eo *> x) noexcept
{
	hmc_float tmp;
	logger.debug() << "===============================================";
	for(uint i=0; i<x.size(); i++){
		tmp = squarenorm(*x[i]);
		logger.debug() << ((i<10) ? " " : "") << "delta_sqnorm[field_" << i << "]: " << std::scientific << std::setprecision(16) << tmp-a[i];		  
	}
	logger.debug() << "===============================================";
}
