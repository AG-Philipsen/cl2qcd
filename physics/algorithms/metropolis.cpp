/** @file
 * Implementation of the metropolis algorithm
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2012-2014 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
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

#include "metropolis.hpp"

#include "solvers/solvers.hpp"
#include "solver_shifted.hpp"
#include "../lattices/util.hpp"
#include "../../meta/util.hpp"
#include <cmath>
#include "../observables/gaugeObservables.h"

static void print_info_debug(const meta::Inputparameters& params, std::string metropolis_part, hmc_float value, bool info=true);

hmc_float physics::algorithms::calc_s_fermion(const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
	using physics::lattices::Spinorfield;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;

	const auto & params = system.get_inputparameters();

	const Spinorfield phi_inv(system);

	logger.debug() << "\t\t\tstart solver";

	/** 
	 * @todo at the moment, we can only put in a cold spinorfield
	 * or a point-source spinorfield as trial-solution
	 */
	const Spinorfield solution(system);
	solution.cold();
	int iterations = 0;

	if(params.get_solver() == common::cg) {
		const QplusQminus fm(kappa, mubar, system);
		iterations  = cg(&solution, fm, gf, phi, system, params.get_solver_prec());
		const Qminus qminus(kappa, mubar, system);
		qminus(&phi_inv, gf, solution);

	} else  {
		const Qplus fm(kappa, mubar, system);
		iterations = bicgstab(&solution, fm, gf, phi, system, params.get_solver_prec());
		copyData(&phi_inv, solution);
	}
	logger.trace() << "Calulated S_fermion, solver took " << iterations << " iterations.";
	return squarenorm(phi_inv);
}

hmc_float physics::algorithms::calc_s_fermion(const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
	using physics::lattices::Spinorfield_eo;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;

	const auto & params = system.get_inputparameters();

	const Spinorfield_eo phi_inv(system);

	logger.debug() << "\t\t\tstart solver";

	/** @todo at the moment, we can only put in a cold spinorfield
	 * or a point-source spinorfield as trial-solution
	 */
	const Spinorfield_eo solution(system);
	int iterations = 0;
	logger.debug() << "\t\t\tstart solver";

	//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
	if(params.get_solver() == common::cg) {
		solution.cold();

		const QplusQminus_eo fm(kappa, mubar, system);
		iterations = cg(&solution, fm, gf, phi, system, params.get_solver_prec());
		const Qminus_eo qminus(kappa, mubar, system);
		qminus(&phi_inv, gf, solution);
	} else {
		solution.zero();
		solution.gamma5();
		const Qplus_eo fm(kappa, mubar, system);
		iterations = bicgstab(&solution, fm, gf, phi, system, params.get_solver_prec());
		copyData(&phi_inv, solution);
	}
	logger.trace() << "Calulated S_fermion, solver took " << iterations << " iterations.";
	return squarenorm(phi_inv);
}

hmc_float physics::algorithms::calc_s_fermion_mp(const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system)
{
	//this function essentially performs the same steps as in the non mass-prec case, however, one has to apply one more matrix multiplication
	//  therefore, comments are deleted here...
	//  Furthermore, in the bicgstab-case, the second inversions are not needed

	using physics::lattices::Spinorfield;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;

	logger.trace() << "\tHMC [DH]:\tcalc final fermion energy...";

	const Spinorfield phi_inv(system);

	const auto & params = system.get_inputparameters();
	const hmc_float kappa = params.get_kappa();
	const hmc_float mubar = meta::get_mubar(params);

	const Spinorfield tmp(system);
	const Qplus qplus_mp(params.get_kappa_mp(), meta::get_mubar_mp(params), system);
	qplus_mp(&tmp, gf, phi);

	const Spinorfield solution(system);
	solution.cold();

	logger.debug() << "\t\t\tstart solver";
	int iterations = 0;

	if(params.get_solver() == common::cg) {
		const QplusQminus fm(kappa, mubar, system);
		iterations = cg(&solution, fm, gf, tmp, system, params.get_solver_prec());
		const Qminus qminus(kappa, mubar, system);
		qminus(&phi_inv, gf, solution);
	} else  {
		const Qplus fm(kappa, mubar, system);
		iterations = bicgstab(&solution, fm, gf, tmp, system, params.get_solver_prec());
		copyData(&phi_inv, solution);
	}
	logger.trace() << "Calulated S_fermion, solver took " << iterations << " iterations.";
	return squarenorm(phi_inv);
}

hmc_float physics::algorithms::calc_s_fermion_mp(const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system)
{
	//this function essentially performs the same steps as in the non mass-prec case, however, one has to apply one more matrix multiplication
	//  therefore, comments are deleted here...
	//  Furthermore, in the bicgstab-case, the second inversions are not needed
	using physics::lattices::Spinorfield_eo;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;

	logger.trace() << "\tHMC [DH]:\tcalc final fermion energy...";

	const auto & params = system.get_inputparameters();
	const hmc_float kappa = params.get_kappa();
	const hmc_float mubar = meta::get_mubar(params);

	const Spinorfield_eo phi_inv(system);

	logger.debug() << "\t\t\tstart solver";
	int iterations = 0;

	//sf_tmp = Qplus(light_mass) phi_mp
	const Spinorfield_eo tmp(system);
	const Qplus_eo qplus_mp(params.get_kappa_mp(), meta::get_mubar_mp(params), system);
	qplus_mp(&tmp, gf, phi);

	const Spinorfield_eo solution(system);

	logger.debug() << "\t\t\tstart solver";

	if(params.get_solver() == common::cg) {
		solution.cold();
		const QplusQminus_eo fm(kappa, mubar, system);
		iterations = cg(&solution, fm, gf, tmp, system, params.get_solver_prec());
		const Qminus_eo qminus(kappa, mubar, system);
		qminus(&phi_inv, gf, solution);
	} else {
		solution.zero();
		solution.gamma5();

		const Qplus_eo fm(kappa, mubar, system);
		iterations = bicgstab(&solution, fm, gf, tmp, system, params.get_solver_prec());

		copyData(&phi_inv, solution);
	}
	logger.trace() << "Calulated S_fermion, solver took " << iterations << " iterations.";
	return squarenorm(phi_inv);
}

hmc_float physics::algorithms::calc_s_fermion_mp(const physics::lattices::Gaugefield&, const physics::lattices::Rooted_Staggeredfield_eo&, const hardware::System&)
{
	throw std::runtime_error("Not implemented!");
}

/**
 * This function returns the value of the fermionic part of the action for the RHMC, i.e.
 * \f$ \phi^*\,(M^\dag\,M)^{-\frac{N_f}{4}}\,\phi \f$.
 * 
 * Here, we use the coefficients of the rational expansion of \f$ x^{-\frac{N_f}{4}} \f$ 
 * included in the Rooted_Staggeredfield_eo object to calculate \f$ (M^\dag\,M)^{-\frac{N_f}{4}}\,\phi \f$
 * (using the multi-shifted inverter). Then a scalar product give the returning value.
 * 
 */
hmc_float physics::algorithms::calc_s_fermion(const physics::lattices::Gaugefield& gf, const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system, const hmc_float mass, const hmc_float)
{
	using physics::lattices::Rooted_Staggeredfield_eo;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;
	
	logger.trace() << "\tRHMC [DH]:\tcalc final fermion energy...";
	
	const auto & params = system.get_inputparameters();
	const physics::fermionmatrix::MdagM_eo fm(system, mass);
	int iterations = 0;

	//Temporary fields for shifted inverter
	logger.debug() << "\t\tstart solver...";
	std::vector<physics::lattices::Staggeredfield_eo *> X;
	for(int i=0; i<phi.Get_order(); i++)
		X.push_back(new physics::lattices::Staggeredfield_eo(system));
	//Here the inversion must be performed with high precision, because it'll be used for Metropolis test
	iterations = physics::algorithms::solvers::cg_m(X, phi.Get_b(), fm, gf, phi, system, params.get_solver_prec());
	logger.debug() << "\t\t...end solver in " << iterations << " iterations";
	
	physics::lattices::Staggeredfield_eo tmp(system); //this is to reconstruct (MdagM)^{-\frac{N_f}{4}}\,\phi
	sax(&tmp, {phi.Get_a0(), 0.}, phi);
	for(int i=0; i<phi.Get_order(); i++){
		saxpy(&tmp, {(phi.Get_a())[i], 0.}, *X[i], tmp);
	}
	
	meta::free_container(X);
	logger.trace() << "Calulated S_fermion, solver took " << iterations << " iterations.";
	return scalar_product(phi, tmp).re;
  
}

template <class SPINORFIELD> static hmc_observables metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf, const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p, const physics::lattices::Gaugemomenta& new_p, const SPINORFIELD& phi, const hmc_float spinor_energy_init, const SPINORFIELD * const phi_mp, const hmc_float spinor_energy_mp_init, const hardware::System& system)
{
	using namespace physics::algorithms;

	const auto & params = system.get_inputparameters();

	//Calc Hamiltonian
	print_info_debug(params, "[DH]:\tCalculate Hamiltonian", sqrt(-1.), false);
	hmc_float deltaH = 0.;
	hmc_float s_old = 0.;
	hmc_float s_new = 0.;

	//Gauge-Part
	hmc_float plaq = physics::observables::measurePlaquetteWithoutNormalization(&gf);
	hmc_float plaq_new = physics::observables::measurePlaquetteWithoutNormalization(&new_u);

	hmc_float rect_new = 0.;
	hmc_float rect = 0.;

	if(meta::get_use_rectangles(params) == true) {
	  rect = physics::observables::measureRectangles(&gf);
	  rect_new = physics::observables::measureRectangles(&new_u);
		hmc_float c0 = meta::get_c0(params);
		hmc_float c1 = meta::get_c1(params);
		deltaH = - beta * ( c0 * (plaq - plaq_new)  + c1 * ( rect - rect_new )  );
		s_old = - beta * ( c0 * ( plaq )  + c1 * ( rect )  );
		s_new = - beta * ( c0 * (plaq_new)  + c1 * ( rect_new )  );
	} else {
		/** NOTE: the minus here is introduced to fit tmlqcd!!! */
		deltaH = -(plaq - plaq_new) * beta ;
		s_old = -(plaq ) * beta ;
		s_new = -(plaq_new) * beta ;
	}

	print_info_debug(params, "[DH]:\tS[GF]_0:\t", s_old, false);
	print_info_debug(params, "[DH]:\tS[GF]_1:\t", s_new, false);
	print_info_debug(params, "[DH]:\tdS[GF]: \t", deltaH);
	//check on NANs
	if (s_old != s_old || s_new != s_new || deltaH != deltaH) {
		throw Print_Error_Message("NAN occured in Metropolis! Aborting!", __FILE__, __LINE__);
	}

	//Gaugemomentum-Part
	hmc_float p2 = squarenorm(p);
	hmc_float new_p2 = squarenorm(new_p);
	if(params.get_fermact() != common::action::rooted_stagg){
		//the energy is half the squarenorm
		deltaH += 0.5 * (p2 - new_p2);
	}else{
		//The gaugemomenta part of the action in the staggered code is:
		// \sum_{n,\mu} 0.5 * Tr[P_\mu(n)P_\mu(n)] = 0.25 \sum_{n,\mu,A} ||P^A_\mu(n)||^2
		    p2 *= 0.5;
		new_p2 *= 0.5; //I multiply here by 0.5 so that in the following cout we have the correct numbers
		deltaH += 0.5 * (p2 - new_p2);
	}
	
	print_info_debug(params, "[DH]:\tS[GM]_0:\t", 0.5 * p2, false);
	print_info_debug(params, "[DH]:\tS[GM]_1:\t", 0.5 * new_p2, false);
	print_info_debug(params, "[DH]:\tdS[GM]: \t", 0.5 * (p2 - new_p2));
	//check on NANs
	if (p2 != p2 || new_p2 != new_p2 || deltaH != deltaH) {
		throw Print_Error_Message("NAN occured in Metropolis! Aborting!", __FILE__, __LINE__);
	}

	//Fermion-Part:
	if(! params.get_use_gauge_only() ) {
		if( params.get_use_mp() ) {
			if(params.get_fermact() == common::action::rooted_stagg) {
				throw Invalid_Parameters("Mass preconditioning not implemented for staggered fermions!", "NOT rooted_stagg", params.get_fermact());
			}
			//in this case one has contributions from det(m_light/m_heavy) and det(m_heavy)
			// det(m_heavy)
			hmc_float s_fermion_final;
			//initial energy has been computed in the beginning...
			s_fermion_final = calc_s_fermion(new_u, phi, system, params.get_kappa_mp(),  meta::get_mubar_mp(params));
			deltaH += spinor_energy_init - s_fermion_final;

			print_info_debug(params, "[DH]:\tS[DET]_0:\t", spinor_energy_init, false);
			print_info_debug(params, "[DH]:\tS[DET]_1:\t", s_fermion_final, false);
			print_info_debug(params, "[DH]:\tdS[DET]:\t" , spinor_energy_init - s_fermion_final);
			//check on NANs
			if (spinor_energy_init != spinor_energy_init || s_fermion_final != s_fermion_final || deltaH != deltaH) {
				throw Print_Error_Message("NAN occured in Metropolis! Aborting!", __FILE__, __LINE__);
			}

			// det(m_light/m_heavy)
			//initial energy has been computed in the beginning...
			hmc_float s_fermion_mp_final = calc_s_fermion_mp(new_u, *phi_mp, system);
			deltaH += spinor_energy_mp_init - s_fermion_mp_final;

			print_info_debug(params, "[DH]:\tS[DETRAT]_0:\t", spinor_energy_mp_init, false);
			print_info_debug(params, "[DH]:\tS[DETRAT]_1:\t", s_fermion_mp_final, false);
			print_info_debug(params, "[DH]:\tdS[DETRAT]:\t", spinor_energy_mp_init - s_fermion_mp_final);
			//check on NANs
			if (spinor_energy_mp_init != spinor_energy_mp_init || s_fermion_mp_final != s_fermion_mp_final || deltaH != deltaH) {
				throw Print_Error_Message("NAN occured in Metropolis! Aborting!", __FILE__, __LINE__);
			}
		} else {
			hmc_float s_fermion_final = calc_s_fermion(new_u, phi, system);
			deltaH += spinor_energy_init - s_fermion_final;

			print_info_debug(params, "[DH]:\tS[DET]_0:\t", spinor_energy_init, false);
			print_info_debug(params, "[DH]:\tS[DET]_1:\t", s_fermion_final, false);
			print_info_debug(params, "[DH]:\tdS[DET]: \t", spinor_energy_init - s_fermion_final);
			//check on NANs
			if (spinor_energy_init != spinor_energy_init || s_fermion_final != s_fermion_final || deltaH != deltaH) {
				throw Print_Error_Message("NAN occured in Metropolis! Aborting!", __FILE__, __LINE__);
			}
		}
	}
	//Metropolis-Part
	hmc_float compare_prob;
	if(deltaH < 0) {
		compare_prob = std::exp(deltaH);
	} else {
		compare_prob = 1.0;
	}
	print_info_debug(params, "[DH]:\tdH:\t\t", deltaH);
	print_info_debug(params, "[MET]:\tAcc-Prop:\t", compare_prob);

	//calc gaugeobservables
	//todo: calc only of final configuration
	auto plaqs = physics::observables::measureAllPlaquettes(&gf);
	plaq = plaqs.plaquette; 
	hmc_float splaq = plaqs.spatialPlaquette;
	hmc_float tplaq = plaqs.temporalPlaquette;

	plaqs = physics::observables::measureAllPlaquettes(&new_u);
	plaq_new = plaqs.plaquette; 
	hmc_float splaq_new = plaqs.spatialPlaquette;
	hmc_float tplaq_new = plaqs.temporalPlaquette;

	hmc_complex poly = physics::observables::measurePolyakovloop(&gf);
	hmc_complex poly_new = physics::observables::measurePolyakovloop(&new_u);

	hmc_observables tmp;
	if(rnd <= compare_prob) {
		tmp.accept = 1;
		tmp.plaq = plaq_new;
		tmp.tplaq = tplaq_new;
		tmp.splaq = splaq_new;
		tmp.poly = poly_new;
		tmp.deltaH = deltaH;
		tmp.prob = compare_prob;
		if(meta::get_use_rectangles(params) ) tmp.rectangles = rect_new / get_rect_norm(params);
	} else {
		tmp.accept = 0;
		tmp.plaq = plaq;
		tmp.tplaq = tplaq;
		tmp.splaq = splaq;
		tmp.poly = poly;
		tmp.deltaH = deltaH;
		tmp.prob = compare_prob;
		if(meta::get_use_rectangles(params) ) tmp.rectangles = rect / get_rect_norm(params);
	}

	return tmp;
}

hmc_observables physics::algorithms::metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf, const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p, const physics::lattices::Gaugemomenta& new_p, const physics::lattices::Spinorfield& phi, const hmc_float spinor_energy_init, const physics::lattices::Spinorfield * const phi_mp, const hmc_float spinor_energy_mp_init, const hardware::System& system)
{
	return ::metropolis(rnd, beta, gf, new_u, p, new_p, phi, spinor_energy_init, phi_mp, spinor_energy_mp_init, system);
}

hmc_observables physics::algorithms::metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf, const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p, const physics::lattices::Gaugemomenta& new_p, const physics::lattices::Spinorfield_eo& phi, const hmc_float spinor_energy_init, const physics::lattices::Spinorfield_eo * const phi_mp, const hmc_float spinor_energy_mp_init, const hardware::System& system)
{
	return ::metropolis(rnd, beta, gf, new_u, p, new_p, phi, spinor_energy_init, phi_mp, spinor_energy_mp_init, system);
}

hmc_observables physics::algorithms::metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf, const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p, const physics::lattices::Gaugemomenta& new_p, const physics::lattices::Rooted_Staggeredfield_eo& phi, const hmc_float spinor_energy_init, const physics::lattices::Rooted_Staggeredfield_eo* const phi_mp, const hmc_float spinor_energy_mp_init, const hardware::System& system)
{
	return ::metropolis(rnd, beta, gf, new_u, p, new_p, phi, spinor_energy_init, phi_mp, spinor_energy_mp_init, system);
}

static void print_info_debug(const meta::Inputparameters& params, std::string metropolis_part, hmc_float value, bool info){

	if(info == false){
	  if(logger.beDebug()){
	    if(params.get_fermact() != common::action::rooted_stagg){
	      if(value == value)
		logger.debug() << "\tHMC " << metropolis_part << std::setprecision(10) << value;
	      else
		logger.debug() << "\tHMC " << metropolis_part;
	    }else{
	      if(value == value)
		logger.debug() << "\tRHMC " << metropolis_part << std::setprecision(10) << value;
	      else
		logger.debug() << "\tRHMC " << metropolis_part;
	    }
	  }
	}else{
	  if(params.get_fermact() != common::action::rooted_stagg){
	    if(value == value)
	      logger.info() << "\tHMC " << metropolis_part << std::setprecision(10) << value;
	    else
	      logger.info() << "\tHMC " << metropolis_part;
	  }else{
	    if(value == value)
	      logger.info() << "\tRHMC " << metropolis_part << std::setprecision(10) << value;
	    else
	      logger.info() << "\tRHMC " << metropolis_part;
	  }
	}

}
