/** @file
 * Implementation of the fermion_force functions
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#include "fermion_force.hpp"

#include "../fermionmatrix/fermionmatrix.hpp"
#include "solver.hpp"
#include "../lattices/util.hpp"
#include "molecular_dynamics.hpp"
#include "../../meta/util.hpp"
#include "../../hardware/device.hpp"
#include "../../hardware/code/molecular_dynamics.hpp"

//this function takes to args kappa and mubar because one has to use it with different masses when mass-prec is used and when not

void physics::algorithms::calc_fermion_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
	using physics::lattices::Spinorfield_eo;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;

	auto params = system.get_inputparameters();

	logger.debug() << "\t\tcalc fermion_force...";
	//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
	Spinorfield_eo solution(system);
	Spinorfield_eo phi_inv(system);
	if(params.get_solver() == meta::Inputparameters::cg) {
		/**
		 * The first inversion calculates
		 * X_even = phi = (Qplusminus_eo)^-1 psi
		 * out of
		 * Qplusminus_eo phi_even = psi
		 */
		logger.debug() << "\t\t\tstart solver";

		/** @todo at the moment, we can only put in a cold spinorfield
		 * or a point-source spinorfield as trial-solution
		 */
		/**
		 * Trial solution for the spinorfield
		 */
		solution.cold();

		const QplusQminus_eo fm(kappa, mubar, system);
		const int iterations = cg(&solution, fm, gf, phi, system, params.get_force_prec());

		/**
		 * Y_even is now just
		 *  Y_even = (Qminus_eo) X_even = (Qminus_eo) (Qplusminus_eo)^-1 psi =
		 *    = (Qplus_eo)^-1 ps"\tinv. field before inversion i
		 */
		const Qminus_eo qminus(kappa, mubar, system);
		qminus(&phi_inv, gf, solution);
	} else {
		///@todo if wanted, solvertimer has to be used here..
		//logger.debug() << "\t\tcalc fermion force ingredients using bicgstab with eo.";
		/**
		 * The first inversion calculates
		 * Y_even = phi = (Qplus_eo)^-1 psi
		 * out of
		 * Qplus_eo phi = psi
		 * This is also the energy of the final field!
		*/
		//logger.debug() << "\t\tcalc Y_even...";
		//logger.debug() << "\t\t\tstart solver";

		/** @todo at the moment, we can only put in a cold spinorfield
		 * or a point-source spinorfield as trial-solution
		 */
		/**
		 * Trial solution for the spinorfield
		 */
		solution.zero();
		solution.gamma5();

		const Qplus_eo qplus(kappa, mubar, system);
		int iterations = bicgstab(&solution, qplus, gf, phi, system, params.get_force_prec());

		copyData(&phi_inv, solution);

		/**
		 * Now, one has to calculate
		 * X = (Qminus_eo)^-1 Y = (Qminus_eo)^-1 (Qplus_eo)^-1 psi = (QplusQminus_eo)^-1 psi ??
		 * out of
		 * Qminus_eo clmem_inout_eo = clmem_phi_inv_eo
		 * Therefore, when calculating the final energy of the spinorfield,
		 *  one can just take clmem_phi_inv_eo (see also above)!!
		 */

		//logger.debug() << "\t\tcalc X_even...";
		//copy former solution to clmem_source
		Spinorfield_eo source_even(system);
		copyData(&source_even, solution);

		//this sets clmem_inout cold as trial-solution
		solution.cold();

		const Qminus_eo qminus(kappa, mubar, system);
		iterations = bicgstab(&solution, qminus, gf, source_even, system, params.get_force_prec());
	}
	/**
	 * At this point, one has calculated X_odd and Y_odd.
	 * If one has a fermionmatrix
	 *  M = R + D
	 * these are:
	 *  X_odd = -R(-mu)_inv D X_even
	 *  Y_odd = -R(mu)_inv D Y_even
	 */

	///@NOTE the following calculations could also go in a new function for convenience
	//calculate X_odd
	//therefore, clmem_tmp_eo_1 is used as intermediate state. The result is saved in clmem_inout, since
	//  this is used as a default in the force-function.
	Spinorfield_eo tmp_1(system);
	if(params.get_fermact() == meta::action::wilson) {
		dslash(&tmp_1, gf, solution, ODD, kappa);
		sax(&tmp_1, { -1., 0.}, tmp_1);
	} else if(params.get_fermact() == meta::action::twistedmass) {
		Spinorfield_eo tmp_2(system);
		dslash(&tmp_1, gf, solution, ODD, kappa);
		M_tm_inverse_sitediagonal_minus(&tmp_2, tmp_1, mubar);
		sax(&tmp_1, { -1., 0.}, tmp_2);
	} else {
		throw Print_Error_Message("The selected fermion force has not been implemented.", __FILE__, __LINE__);
	}
	//logger.debug() << "\t\tcalc eo fermion_force F(Y_even, X_odd)...";
	//Calc F(Y_even, X_odd) = F(clmem_phi_inv_eo, clmem_tmp_eo_1)
	fermion_force(force, phi_inv, tmp_1, EVEN, gf, kappa);

	//calculate Y_odd
	//therefore, clmem_tmp_eo_1 is used as intermediate state. The result is saved in clmem_phi_inv, since
	//  this is used as a default in the force-function.
	if(params.get_fermact() == meta::action::wilson) {
		dslash(&tmp_1, gf, phi_inv, ODD, kappa);
		sax(&tmp_1, { -1., 0.}, tmp_1);
	} else if(params.get_fermact() == meta::action::twistedmass) {
		Spinorfield_eo tmp_2(system);
		dslash(&tmp_1, gf, phi_inv, ODD, kappa);
		M_tm_inverse_sitediagonal(&tmp_2, tmp_1, mubar);
		sax(&tmp_1, { -1., 0.}, tmp_2);
	} else {
		throw Print_Error_Message("The selected fermion force has not been implemented.", __FILE__, __LINE__);
	}
	//logger.debug() << "\t\tcalc eoprec fermion_force F(Y_odd, X_even)...";
	//Calc F(Y_odd, X_even) = F(clmem_tmp_eo_1, clmem_inout_eo)
	fermion_force(force, tmp_1, solution, ODD, gf, kappa);
}


void physics::algorithms::calc_fermion_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
	using physics::lattices::Spinorfield;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;

	auto params = system.get_inputparameters();

	Spinorfield solution(system);
	Spinorfield phi_inv(system);

	logger.debug() << "\t\tcalc fermion_force...";
	//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
	if(params.get_solver() == meta::Inputparameters::cg) {
		/**
		 * The first inversion calculates
		 * X = phi = (Qplusminus)^-1 psi
		 * out of
		 * Qplusminus phi = psi
		 */
		logger.debug() << "\t\t\tstart solver";

		/** @todo at the moment, we can only put in a cold spinorfield
		 * or a point-source spinorfield as trial-solution
		 */
		/**
		 * Trial solution for the spinorfield
		 */
		solution.cold();

		//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
		const QplusQminus fm(kappa, mubar, system);
		const int iterations = cg(&solution, fm, gf, phi, system, params.get_force_prec());

		/**
		 * Y is now just
		 *  Y = (Qminus) X = (Qminus) (Qplusminus)^-1 psi =
		 *    = (Qplus)^-1 psi
		 */
		const Qminus qminus(kappa, mubar, system);
		qminus(&phi_inv, gf, solution);

	} else  {
		logger.debug() << "\t\tcalc fermion force ingredients using bicgstab without eo";

		/**
		 * The first inversion calculates
		 * Y = phi = (Qplus)^-1 psi
		 * out of
		 * Qplus phi = psi
		 * This is also the energy of the final field!
		 */
		logger.debug() << "\t\t\tstart solver";

		/** @todo at the moment, we can only put in a cold spinorfield
		 * or a point-source spinorfield as trial-solution
		 */
		/**
		 * Trial solution for the spinorfield
		 */
		solution.cold();

		//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
		const Qplus q_plus(kappa, mubar, system);
		int iterations = bicgstab(&solution, q_plus, gf, phi, system, params.get_force_prec());

		//store this result in clmem_phi_inv
		copyData(&phi_inv, solution);

		/**
		 * Now, one has to calculate
		 * X = (Qminus)^-1 Y = (Qminus)^-1 (Qplus)^-1 psi = (QplusQminus)^-1 psi ??
		 * out of
		 * Qminus clmem_inout = clmem_phi_inv
		 * Therefore, when calculating the final energy of the spinorfield,
		 *  one can just take clmem_phi_inv (see also above)!!
		 */

		//copy former solution to clmem_source
		Spinorfield source(system);
		copyData(&source, solution);

		logger.debug() << "\t\t\tstart solver";

		solution.cold();
		const Qminus q_minus(kappa, mubar, system);
		iterations = bicgstab(&solution, q_minus, gf, source, system, params.get_force_prec());
	}
	log_squarenorm("\tY ", phi_inv);
	log_squarenorm("\tX ", solution);

	logger.debug() << "\t\tcalc fermion_force...";
	fermion_force(force, phi_inv, solution, gf, kappa);
}


void physics::algorithms::calc_fermion_force_detratio(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi_mp, const hardware::System& system)
{
	using physics::lattices::Spinorfield;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;

	logger.debug() << "\t\tcalc fermion_force_detratio...";

	auto params = system.get_inputparameters();

	/**
	 * For detratio = det(kappa, mubar) / det(kappa2, mubar2) = det(Q_1^+Q_1^-) / det(Q_2^+Q_2^-)
	 * the force has almost the same ingredients as in the above case:
	 *   F(detratio) = - ( - phi^+ deriv(Q_2) X + Y^+ deriv(Q_1) X  ) + h.c.;
	 * where deriv(Q_i) is the same fct. as above with different parameters,
	 * X = (Q_1^+ Q_1^-)^-1 Q_2^+ phi
	 * Y = (Q_1^+)^-1 Q_2^+ phi
	 * (In the case of Q_2 = 1 = const., one recovers the expression for the "normal" force
	 * The main differences are:
	 *   - invert Q_2^+ phi, not phi for X and Y
	 *   - one additional force term with different mass-parameters and without Y
	 */
	const hmc_float kappa = params.get_kappa();
	const hmc_float mubar = meta::get_mubar(params);
	const hmc_float kappa2 = params.get_kappa_mp();
	const hmc_float mubar2 = meta::get_mubar_mp(params);

	const Spinorfield phi_inv(system);
	const Spinorfield solution(system);

	//CP: Init tmp spinorfield
	const Spinorfield tmp(system);
	//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
	//the source is now Q_2^+ phi = sf_tmp
	const Qplus q_plus_mp(kappa2, mubar2, system);
	q_plus_mp(&tmp, gf, phi_mp);

	if(params.get_solver() == meta::Inputparameters::cg) {
		/**
		 * The first inversion calculates
		 * X = phi = (Qplusminus)^-1 sf_tmp
		 * out of
		 * Qplusminus phi = sf_tmp
		 */
		logger.debug() << "\t\t\tstart solver";

		/** @todo at the moment, we can only put in a cold spinorfield
		 * or a point-source spinorfield as trial-solution
		 */
		/**
		 * Trial solution for the spinorfield
		 */
		solution.cold();

		//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
		const QplusQminus fm(kappa, mubar, system);
		const int iterations = cg(&solution, fm, gf, tmp, system, params.get_force_prec());

		/**
		 * Y is now just
		 *  Y = (Qminus) X = (Qminus) (Qplusminus)^-1 sf_tmp =
		 *    = (Qplus)^-1 sf_tmp
		 */
		const Qminus q_minus(kappa, mubar, system);
		q_minus(&phi_inv, gf, solution);

	} else {
		logger.debug() << "\t\tcalc fermion force ingredients using bicgstab without eo";

		/**
		 * The first inversion calculates
		 * Y = phi = (Qplus)^-1 sf_tmp
		 * out of
		 * Qplus phi = sf_tmp
		 * This is also the energy of the final field!
		 */
		logger.debug() << "\t\t\tstart solver";

		/** @todo at the moment, we can only put in a cold spinorfield
		 * or a point-source spinorfield as trial-solution
		 */
		/**
		 * Trial solution for the spinorfield
		 */
		solution.cold();

		//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
		Qplus q_plus(kappa, mubar, system);
		int iterations = bicgstab(&solution, q_plus, gf, tmp, system, params.get_force_prec());

		copyData(&phi_inv, solution);

		/**
		 * Now, one has to calculate
		 * X = (Qminus)^-1 Y = (Qminus)^-1 (Qplus)^-1 psi = (QplusQminus)^-1 psi ??
		 * out of
		 * Qminus clmem_inout = clmem_phi_inv
		 * Therefore, when calculating the final energy of the spinorfield,
		 *  one can just take clmem_phi_inv (see also above)!!
		 */

		//copy former solution to clmem_source
		const Spinorfield source(system);
		copyData(&source, solution);
		logger.debug() << "\t\t\tstart solver";

		//this sets clmem_inout cold as trial-solution
		solution.cold();

		Qminus q_minus(kappa, mubar, system);
		iterations = bicgstab(&solution, q_minus, gf, source, system, params.get_force_prec());
	}

	log_squarenorm("\tY ", phi_inv);
	log_squarenorm("\tX ", solution);

	logger.debug() << "\t\tcalc fermion_force...";
	fermion_force(force, phi_inv, solution, gf, kappa);

	/**
	 *Now, one has the additional term - phi^+ deriv(Q_2) X
	 *This works exactly the same as above, except replacing Y -> -phi and different mass parameters
	 */

	//Y is not needed anymore, therefore use clmem_phi_inv_eo to store -phi
	sax(&phi_inv, { -1., 0.}, phi_mp);

	fermion_force(force, phi_inv, solution, gf, kappa2);
}

void physics::algorithms::calc_fermion_force_detratio(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system)
{
	using physics::lattices::Spinorfield_eo;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;

	logger.debug() << "\t\tcalc fermion_force_detratio...";

	auto params = system.get_inputparameters();

	/**
	 * For detratio = det(kappa, mubar) / det(kappa2, mubar2) = det(Q_1^+Q_1^-) / det(Q_2^+Q_2^-)
	 * the force has almost the same ingredients as in the above case:
	 *   F(detratio) = - ( - phi^+ deriv(Q_2) X + Y^+ deriv(Q_1) X  ) + h.c.;
	 * where deriv(Q_i) is the same fct. as above with different parameters,
	 * X = (Q_1^+ Q_1^-)^-1 Q_2^+ phi
	 * Y = (Q_1^+)^-1 Q_2^+ phi
	 * (In the case of Q_2 = 1 = const., one recovers the expression for the "normal" force
	 * The main differences are:
	 *   - invert Q_2^+ phi, not phi for X and Y
	 *   - one additional force term with different mass-parameters and without Y
	 */
	const hmc_float kappa = params.get_kappa();
	const hmc_float mubar = get_mubar(params);
	const hmc_float kappa2 = params.get_kappa_mp();
	const hmc_float mubar2 = meta::get_mubar_mp(params);

	const Spinorfield_eo solution(system);
	const Spinorfield_eo phi_inv(system);

	const physics::lattices::Scalar<hmc_complex> mone(system);
	mone.store( { -1., 0.});

	const Spinorfield_eo tmp(system);
	//the source is now Q_2^+ phi = sf_eo_tmp
	const Qplus_eo q_plus_mp(kappa2, mubar2, system);
	q_plus_mp(&tmp, gf, phi_mp);
	if(params.get_solver() == meta::Inputparameters::cg) {
		/**
		 * The first inversion calculates
		 * X_even = phi = (Qplusminus_eo)^-1 sf_eo_tmp = (Qplusminus_eo)^-1 Q_2^+ phi
		 * out of
		 * Qplusminus_eo phi = sf_eo_tmp
		 */
		logger.debug() << "\t\t\tstart solver";

		/** @todo at the moment, we can only put in a cold spinorfield
		 * or a point-source spinorfield as trial-solution
		 */
		/**
		 * Trial solution for the spinorfield
		 */
		solution.cold();

		const QplusQminus_eo fm(kappa, mubar, system);
		const int iterations = cg(&solution, fm, gf, tmp, system, params.get_force_prec());

		/**
		 * Y_even is now just
		 *  Y_even = (Qminus_eo) X_even = (Qminus_eo) (Qplusminus_eo)^-1 sf_eo_tmp =
		 *    = (Qplus_eo)^-1 Q_2^+ psi
		 */
		const Qminus_eo q_minus(kappa, mubar, system);
		q_minus(&phi_inv, gf, solution);
	} else {
		///@todo if wanted, solvertimer has to be used here..
		//logger.debug() << "\t\tcalc fermion force ingredients using bicgstab with eo.";
		/**
		 * The first inversion calculates
		 * Y_even = phi = (Qplus_eo)^-1 psi
		 * out of
		 * Qplus_eo phi = psi
		 */
		logger.debug() << "\t\tcalc Y_even...";
		logger.debug() << "\t\t\tstart solver";

		/** @todo at the moment, we can only put in a cold spinorfield
		 * or a point-source spinorfield as trial-solution
		 */
		/**
		 * Trial solution for the spinorfield
		 */
		solution.zero();
		solution.gamma5();

		const Qplus_eo q_plus(kappa, mubar, system);
		int iterations = bicgstab(&solution, q_plus, gf, tmp, system, params.get_force_prec());

		//store this result in clmem_phi_inv
		copyData(&phi_inv, solution);

		/**
		 * Now, one has to calculate
		 * X = (Qminus_eo)^-1 Y = (Qminus_eo)^-1 (Qplus_eo)^-1 sf_eo_tmp = (QplusQminus_eo)^-1 sf_eo_tmp = (QplusQminus_eo)^-1 Q_2^+ psi
		 * out of
		 * Qminus_eo clmem_inout_eo = clmem_phi_inv_eo
		 */

		logger.debug() << "\t\tcalc X_even...";
		//copy former solution to clmem_source
		const Spinorfield_eo source_even(system);
		copyData(&source_even, solution);

		//this sets clmem_inout cold as trial-solution
		solution.cold();

		Qminus_eo q_minus(kappa, mubar, system);
		iterations = bicgstab(&solution, q_minus, gf, source_even, system, params.get_force_prec());
	}
	/**
	 * At this point, one has to calculate X_odd and Y_odd.
	 * If one has a fermionmatrix
	 *  M = R + D
	 * these are:
	 *  X_odd = -R(-mu)_inv D X_even
	 *  Y_odd = -R(mu)_inv D Y_even
	 */

	const Spinorfield_eo tmp1(system);
	const Spinorfield_eo tmp2(system);

	///@NOTE the following calculations could also go in a new function for convenience
	//calculate X_odd
	//therefore, sf_eo_tmp is used as intermediate state.
	if(params.get_fermact() == meta::action::wilson) {
		dslash(&tmp1, gf, solution, ODD);
		sax(&tmp, mone, tmp1);
	} else if(params.get_fermact() == meta::action::twistedmass) {
		dslash(&tmp1, gf, solution, ODD);
		M_tm_inverse_sitediagonal_minus(&tmp2, tmp1);
		sax(&tmp, mone, tmp2);
	} else {
		throw Print_Error_Message("Selected fermion action is not implemented", __FILE__, __LINE__);
	}

	//logger.debug() << "\t\tcalc eo fermion_force F(Y_even, X_odd)...";
	//Calc F(Y_even, X_odd) = F(clmem_phi_inv_eo, sf_eo_tmp)
	fermion_force(force, phi_inv, tmp, EVEN, gf, kappa);

	//calculate Y_odd
	//therefore, clmem_tmp_eo_1 is used as intermediate state.
	if(params.get_fermact() == meta::action::wilson) {
		dslash(&tmp1, gf, phi_inv, ODD);
		sax(&tmp1, mone, tmp1);
	} else if(params.get_fermact() == meta::action::twistedmass) {
		dslash(&tmp1, gf, phi_inv, ODD);
		M_tm_inverse_sitediagonal(&tmp2, tmp1);
		sax(&tmp1, mone, tmp2);
	} else {
		throw Print_Error_Message("Selected fermion action is not implemented", __FILE__, __LINE__);
	}

	//logger.debug() << "\t\tcalc eoprec fermion_force F(Y_odd, X_even)...";
	//Calc F(Y_odd, X_even) = F(clmem_tmp_eo_1, clmem_inout_eo)
	fermion_force(force, tmp1, solution, ODD, gf, kappa);

	/**
	 *Now, one has the additional term - phi^+ deriv(Q_2) X
	 *This works exactly the same as above, except replacing Y -> -phi and different mass parameters
	 */

	//Y is not needed anymore, therefore use clmem_phi_inv_eo to store -phi
	sax(&phi_inv, mone, phi_mp);

	//(re-) calculate X_odd (with other masses)
	//therefore, sf_eo_tmp is used as intermediate state.
	if(params.get_fermact() == meta::action::wilson) {
		dslash(&tmp1, gf, solution, ODD, kappa2);
		sax(&tmp, mone, tmp1);
	} else if(params.get_fermact() == meta::action::twistedmass) {
		dslash(&tmp1, gf, solution, ODD, kappa2);
		M_tm_inverse_sitediagonal_minus(&tmp2, tmp1, mubar2);
		sax(&tmp, mone, tmp2);
	} else {
		throw Print_Error_Message("Selected fermion action is not implemented", __FILE__, __LINE__);
	}

	//logger.debug() << "\t\tcalc eo fermion_force F(Y_even, X_odd)...";
	//Calc F(Y_even, X_odd) = F(clmem_phi_inv_eo, sf_eo_tmp)
	fermion_force(force, phi_inv, tmp, EVEN, gf, kappa2);

	//calculate phi_odd
	//this works in the same way as with Y above, since -phi_even is saved in the same buffer as Y_even
	//therefore, clmem_tmp_eo_1 is used as intermediate state.
	if(params.get_fermact() == meta::action::wilson) {
		dslash(&tmp1, gf, phi_inv, ODD, kappa2);
		sax(&tmp1, mone, tmp1);
	} else if(params.get_fermact() == meta::action::twistedmass) {
		dslash(&tmp1, gf, phi_inv, ODD, kappa2);
		M_tm_inverse_sitediagonal(&tmp2, tmp1, mubar2);
		sax(&tmp1, mone, tmp2);
	} else {
		throw Print_Error_Message("Selected fermion action is not implemented", __FILE__, __LINE__);
	}

	//logger.debug() << "\t\tcalc eoprec fermion_force F(Y_odd, X_even)...";
	//Calc F(Y_odd, X_even) = F(clmem_tmp_eo_1, clmem_inout_eo)
	fermion_force(force, tmp1, solution, ODD, gf, kappa2);
}

template<class SPINORFIELD> static void calc_fermion_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const SPINORFIELD& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	using physics::lattices::Gaugefield;
	using namespace physics::algorithms;

	auto params = system.get_inputparameters();
	//in case of stout-smearing we need every intermediate field for the force calculation
	//NOTE: if smearing is not used, this is just 0
	const int rho_iter = params.get_rho_iter();
	//array to save the intermediate fields
	//NOTE: One needs only rho_iter -1 here since the last iteration is saved in gf...
	//NOTE: If the original gf is also needed in the force calculation, one has to add it here
	//  or use the intermediate cl_mem obj gf_unsmeared. This is initialized in the smear_gaugefield function
	calc_fermion_force(force, gf, phi, system, kappa, mubar);
	if(params.get_use_smearing() == true) {
		throw Print_Error_Message("Smeared Gaugefield force is not implemented.", __FILE__, __LINE__);
		//  mol_dyn_code->stout_smeared_fermion_force_device(smeared_gfs);
		//  gf_code->unsmear_gaugefield(hmc_code->get_new_u());
	}
}

void physics::algorithms::calc_fermion_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	::calc_fermion_forces(force, gf, phi, system, kappa, mubar);
}

void physics::algorithms::calc_fermion_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	::calc_fermion_forces(force, gf, phi, system, kappa, mubar);
}

void physics::algorithms::fermion_force(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Spinorfield& Y, const physics::lattices::Spinorfield& X, const physics::lattices::Gaugefield& gf, const hmc_float kappa)
{
	auto gm_bufs = gm->get_buffers();
	auto Y_bufs = Y.get_buffers();
	auto X_bufs = X.get_buffers();
	auto gf_bufs = gf.get_buffers();
	size_t num_bufs = gm_bufs.size();
	if(num_bufs != Y_bufs.size() || num_bufs != X_bufs.size() || num_bufs != gf_bufs.size()) {
		throw Print_Error_Message(std::string(__func__) + " is only implemented for a single device.", __FILE__, __LINE__);
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto gm_buf = gm_bufs[i];
		auto Y_buf = Y_bufs[i];
		auto X_buf = X_bufs[i];
		auto gf_buf = gf_bufs[i];
		auto code = gm_buf->get_device()->get_molecular_dynamics_code();
		code->fermion_force_device(Y_buf, X_buf, gf_buf, gm_buf, kappa);
	}
	gm->update_halo();
}

void physics::algorithms::fermion_force(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Spinorfield_eo& Y, const physics::lattices::Spinorfield_eo& X, const int evenodd, const physics::lattices::Gaugefield& gf, const hmc_float kappa)
{
	Y.require_halo();
	X.require_halo();

	auto gm_bufs = gm->get_buffers();
	auto Y_bufs = Y.get_buffers();
	auto X_bufs = X.get_buffers();
	auto gf_bufs = gf.get_buffers();
	size_t num_bufs = gm_bufs.size();
	if(num_bufs != Y_bufs.size() || num_bufs != X_bufs.size() || num_bufs != gf_bufs.size()) {
		throw Print_Error_Message(std::string(__func__) + " is only implemented for a single device.", __FILE__, __LINE__);
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto gm_buf = gm_bufs[i];
		auto Y_buf = Y_bufs[i];
		auto X_buf = X_bufs[i];
		auto gf_buf = gf_bufs[i];
		auto code = gm_buf->get_device()->get_molecular_dynamics_code();
		code->fermion_force_eo_device(Y_buf, X_buf, gf_buf, gm_buf, evenodd, kappa);
	}

	gm->update_halo();
}


template<class SPINORFIELD> static void calc_detratio_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const SPINORFIELD& phi_mp, const hardware::System& system)
{
	using physics::lattices::Gaugefield;
	using namespace physics::algorithms;

	auto params = system.get_inputparameters();
	//in case of stout-smearing we need every intermediate field for the force calculation
	//NOTE: if smearing is not used, this is just 0
	const int rho_iter = params.get_rho_iter();
	//array to save the intermediate fields
	//NOTE: One needs only rho_iter -1 here since the last iteration is saved in gf...
	//NOTE: If the original gf is also needed in the force calculation, one has to add it here
	//  or use the intermediate cl_mem obj gf_unsmeared. This is initialized in the smear_gaugefield function
	calc_fermion_force_detratio(force, gf, phi_mp, system);
	if(params.get_use_smearing() == true) {
		throw Print_Error_Message("Smeared Gaugefield force is not implemented.", __FILE__, __LINE__);
		//  mol_dyn_code->stout_smeared_fermion_force_device(smeared_gfs);
		//  gf_code->unsmear_gaugefield(hmc_code->get_new_u());
	}
}

void physics::algorithms::calc_detratio_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi_mp, const hardware::System& system)
{
	::calc_detratio_forces(force, gf, phi_mp, system);
}
void physics::algorithms::calc_detratio_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system)
{
	::calc_detratio_forces(force, gf, phi_mp, system);
}
