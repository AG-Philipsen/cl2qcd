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

#include "../../logger.hpp"
//#include "../../operations_complex.h"
#include "../../meta/type_ops.hpp"
//#include "../lattices/util.hpp"
#include "../lattices/scalar_complex.hpp"
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
	if(squarenorm(b)==0){
	  for(uint i=0; i<sigma.size(); i++)
		x[i]->set_zero();   
	  return 0;
	}
	
	std::vector<hmc_float> xsq(x.size(), 0);
	
	if(sigma.size() != x.size())
		throw std::invalid_argument("Wrong size of multi-shifted inverter parameters!");
	
	using namespace physics::lattices;
	using physics::algorithms::solvers::SolverStuck;
	using physics::algorithms::solvers::SolverDidNotSolve;

	auto params = system.get_inputparameters();
	
	/// @todo start timer synchronized with device(s)
	klepsydra::Monotonic timer;
	/// @todo make configurable from outside
	const bool USE_ASYNC_COPY = params.get_cg_use_async_copy();
	if(USE_ASYNC_COPY) {
		logger.warn() << "Asynchroneous copying in the CG-M is currently unimplemented!";
	}

	//The number of values of constants sigma as well as the number of fields x
	const int Neqs = sigma.size();
	
	//Auxiliary staggered fields
	const Staggeredfield_eo* r = new Staggeredfield_eo(system);
	const Staggeredfield_eo* p = new Staggeredfield_eo(system);
	std::vector<Staggeredfield_eo*> ps;
	
	//Auxiliary scalar vectors
	std::vector<const Scalar<hmc_complex>*> alpha;
	std::vector<const Scalar<hmc_complex>*> beta;
	std::vector<const Scalar<hmc_complex>*> zeta_i;   //This is zeta at the step iter-1
	std::vector<const Scalar<hmc_complex>*> zeta_ii;  //This is zeta at the step iter
	std::vector<const Scalar<hmc_complex>*> zeta_iii; //This is zeta at the step iter+1
	std::vector<const Scalar<hmc_complex>*> shift;    //This is to store constants sigma
	std::vector<bool> single_system_converged;        //This is to stop calculation on single system
	std::vector<uint> single_system_iter;             //This is to calculate performance properly
	//Auxiliary scalars
	const Scalar<hmc_complex> alpha_scalar_prev(system);   //This is alpha_scalar at the step iter-1
	const Scalar<hmc_complex> alpha_scalar(system);        //This is alpha_scalar at the step iter
	const Scalar<hmc_complex> beta_scalar_prev(system);    //This is beta_scalar at the step iter-1
	const Scalar<hmc_complex> beta_scalar(system);         //This is beta_scalar at the step iter
	
	//Auxiliary containers for temporary saving
	const Staggeredfield_eo* v = new Staggeredfield_eo(system); //this is to store A.p
	const Scalar<hmc_complex> tmp1(system);                     //this is to store (r,r) before updating r
	const Scalar<hmc_complex> tmp2(system);                     //this is to store (r,r) after updating r
	const Scalar<hmc_complex> tmp3(system);                     //this is to store (p,v) as Scalar
	const Scalar<hmc_complex> num(system);                      //this is to store constants numerators
	const Scalar<hmc_complex> den(system);                      //this is to store constants denumerators
	bool converged=false;                                       //this is to start to check the residuum

	//Auxiliary constants as Scalar
	const Scalar<hmc_complex> zero(system);
	zero.store(hmc_complex_zero);
	const Scalar<hmc_complex> one(system);
	one.store(hmc_complex_one);
	
	
	hmc_float resid;
	int iter = 0;
	
	//Initialization auxilary and output quantities
	for(int i=0; i<Neqs; i++){
		x[i]->set_zero();                                     // x[i] = 0 
		ps.push_back(new Staggeredfield_eo(system));    
		copyData(ps[i], b);                                    // ps[i] = b
		beta.push_back(new Scalar<hmc_complex>(system));
		alpha.push_back(new Scalar<hmc_complex>(system));
		alpha[i]->store(hmc_complex_zero);                    // alpha[i] = 0
		zeta_i.push_back(new Scalar<hmc_complex>(system));
		zeta_ii.push_back(new Scalar<hmc_complex>(system));
		zeta_iii.push_back(new Scalar<hmc_complex>(system));
		zeta_i[i]->store(hmc_complex_one);                    // zeta_i[i] = 1
		zeta_ii[i]->store(hmc_complex_one);                   // zeta_ii[i] = 1
		shift.push_back(new Scalar<hmc_complex>(system));
		shift[i]->store({sigma[i],0.});
		single_system_converged.push_back(false);             //no system converged
	}
	copyData(r, b);                          // r = b
	copyData(p, b);                          // p = b
	scalar_product(&tmp1, *r, *r);           // set tmp1 = (r, r) for the first iteration
	beta_scalar.store(hmc_complex_one);      // beta_scalar = 1, here I should set beta_scalar_prev
						  // but in this way I can set beta_scalar_prev at the begin
						  // of the loop over iter recursively.
	alpha_scalar.store(hmc_complex_zero);    // alpha_scalar = 0. The same as beta_scalar above.
	
	//At first, in the log string I will report all the data about the system.
	//To avoid to print to shell tens of lines per time, since Neqs can be 
	//also 20 or something like that, reduce report_num.
	const int report_num = Neqs;
	if(report_num>Neqs)
		throw std::invalid_argument("In cg-m report_num cannot be bigger than Neqs!");
	log_squarenorm(create_log_prefix_cgm(iter) + "b (initial): ", b);
	log_squarenorm(create_log_prefix_cgm(iter) + "r (initial): ", *r);
	log_squarenorm(create_log_prefix_cgm(iter) + "p (initial): ", *p);
	log_squarenorm_aux(create_log_prefix_cgm(iter) + "x (initial)", x, report_num);
	log_squarenorm_aux(create_log_prefix_cgm(iter) + "ps (initial)", ps, report_num);
	
	//NOTE: Here, most of the complex numbers may also be just hmc_floats.
	//      However, for this one would need some additional functions in hardware::code...
	for(iter = 0; iter < params.get_cgmax(); iter ++) {
		//Update beta_scalar: v=A.p and tmp1=(r,r) and tmp3=(p,v) ---> beta_scalar=(-1)*tmp1/tmp3
		copyData(&beta_scalar_prev,beta_scalar);  //before updating beta_scalar its value is saved
		A(v, gf, *p);
		log_squarenorm(create_log_prefix_cgm(iter) + "v: ", *v);
		scalar_product(&tmp3, *p, *v); 
		divide(&beta_scalar, tmp1, tmp3);  //tmp1 is set from previous iteration
		subtract(&beta_scalar, zero, beta_scalar);
		//Update field r: r+=beta_scalar*A.p ---> r = r + beta_scalar*v
		saxpy(r, beta_scalar, *v, *r);
		log_squarenorm(create_log_prefix_cgm(iter) + "r: ", *r);
		//We store in tmp2 the quantity (r,r) that we use later. When we check
		//the residuum, then it is already calculated.
		scalar_product(&tmp2, *r, *r);
		if(logger.beDebug()){
			//Calculate squarenorm of the output field
			for(uint i=0; i<x.size(); i++)
				xsq[i] = squarenorm(*x[i]);
		}
		//Update alpha_scalar: alpha_scalar = tmp2/tmp1
		copyData(&alpha_scalar_prev,alpha_scalar); //before updating alpha_scalar its value is saved
		divide(&alpha_scalar, tmp2, tmp1);
		//Update field p: p = r + alpha_scalar*p
		saxpy(p, alpha_scalar, *p, *r);
		//saxpy(p, alpha_scalar, *p, *r);
		log_squarenorm(create_log_prefix_cgm(iter) + "p: ", *p);
		//Loop over the system equations, namely over the set of sigma values
		for(int k=0; k<Neqs; k++){
			if(single_system_converged[k]==false){
				//Update zeta_iii[k]: num = zeta_i[k]*zeta_ii[k]*beta_scalar_prev
				//                    den = beta_scalar*alpha_scalar*(zeta_i[k]-zeta_ii[k])+
				//                        + zeta_i[k]*beta_scalar_prev*(1-sigma[k]*beta_scalar)
				// ---> zeta_iii[k] = num/den
				//I calculate before the denominator to can use num as auxiliary variable
				subtract(&num, *zeta_i[k], *zeta_ii[k]);
				multiply(&num, num, alpha_scalar_prev);
				multiply(&num, num, beta_scalar);
				multiply(&den, *shift[k], beta_scalar);
				subtract(&den, one, den);
				multiply(&den, den, beta_scalar_prev);
				multiply(&den, den, *zeta_i[k]);
				add(&den, num, den);
				//Calculation of the numerator
				multiply(&num, *zeta_i[k], *zeta_ii[k]);
				multiply(&num, num, beta_scalar_prev);
				//Update of zeta_iii[k]
				divide(zeta_iii[k], num, den);
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
				saxpby(ps[k], *zeta_iii[k], *r, *alpha[k], *ps[k]);
				//Check fields squarenorm form possible nan
				if((squarenorm(*x[k]) != squarenorm(*x[k])) ||
				   (squarenorm(*ps[k]) != squarenorm(*ps[k]))){
					logger.fatal() << create_log_prefix_cgm(iter) << "NAN occured!";
					throw SolverStuck(iter, __FILE__, __LINE__);
				}
				//Check if single system converged: ||zeta_iii[k] * r||^2 < prec
				// ---> v = zeta_iii[k] * r
				sax(v, *zeta_iii[k], *r);
				//Check if single system converged: (zeta_ii[k] * zeta_ii[k] * tmp1) < prec
				//multiply(&tmp3, *zeta_ii[k], *zeta_ii[k]);
				//multiply(&tmp3, tmp3, tmp1);
				if(squarenorm(*v) < prec){
				//if(tmp3.get().re < prec){
					single_system_converged[k] = true;
					single_system_iter.push_back((uint)iter);
					logger.debug() << " ===> System number " << k << " converged after " << iter << " iterations! resid = " << tmp2.get().re;
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
			resid = tmp2.get().re;
			logger.debug() << create_log_prefix_cgm(iter) << "resid: " << resid;
			//test if resid is NAN
			if(resid != resid) {
				logger.fatal() << create_log_prefix_cgm(iter) << "NAN occured!";
				throw SolverStuck(iter, __FILE__, __LINE__);
			}
			if(resid < prec) {
				logger.debug() << create_log_prefix_cgm(iter) << "Solver converged in " << iter << " iterations! resid:\t" << resid;
				// report on performance
				if(logger.beInfo()) {
					// we are always synchroneous here, as we had to recieve the residuum from the device
					const uint64_t duration = timer.getTime();
					// calculate flops
					const cl_ulong matrix_flops = A.get_flops();
					const int sum_of_partial_iter = std::accumulate(single_system_iter.begin(), single_system_iter.end(),0);
					logger.debug() << "matrix_flops: " << matrix_flops;
					cl_ulong total_flops = iter * (matrix_flops +
						2 * get_flops<Staggeredfield_eo, scalar_product>(system) +
						2 * ::get_flops<hmc_complex, complexdivide>() +
						    ::get_flops<hmc_complex, complexsubtract>() +
						2 * get_flops<Staggeredfield_eo, saxpy>(system)) +
						sum_of_partial_iter * (
							 3 * ::get_flops<hmc_complex, complexdivide>() +
							11 * ::get_flops<hmc_complex, complexmult>() +
							     ::get_flops<hmc_complex, complexadd>() +
							 3 * ::get_flops<hmc_complex, complexsubtract>() +
							       get_flops<Staggeredfield_eo, saxpy>(system) +
							       get_flops<Staggeredfield_eo, saxpby>(system) +
							       get_flops<Staggeredfield_eo, sax>(system) +
							   5 * get_flops<Staggeredfield_eo, squarenorm>(system));
					logger.debug() << "total_flops: " << total_flops;
					// report performance
					logger.info() << create_log_prefix_cgm(iter) << "CG-M completed in " << std::setprecision(6) << duration / 1000.f << " ms @ " << ((hmc_float)total_flops / duration / 1000.f) << " Gflops. Performed " << iter << " iterations.";
				}
				// report on solution
				log_squarenorm_aux(create_log_prefix_cgm(iter) + "x (final): ", x, report_num);
				
				//for(Staggeredfield_eo* i : x)
				//  logger.warn() << "sqnorm(x[i]) = " << squarenorm(*i);
				return iter;
			}
		}
		//Set tmp1 from tmp2 for following iteration
		copyData(&tmp1, tmp2);
	}
	
	logger.fatal() << create_log_prefix_cgm(iter) << "Solver did not solve in " << params.get_cgmax() << " iterations. Last resid: " << tmp2.get().re;
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
  
  
