/** @file
 * Implementation of the metropolis algorithm
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * (c) 2012-2013 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
 */

#include "metropolis.hpp"

#include "solver.hpp"
#include "../lattices/util.hpp"
#include "../../meta/util.hpp"
#include <cmath>

hmc_float physics::algorithms::calc_s_fermion(const physics::lattices::Spinorfield * const phi_inv, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
	using physics::lattices::Spinorfield;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;

	auto params = system.get_inputparameters();

	logger.debug() << "\t\t\tstart solver";

	/** @todo at the moment, we can only put in a cold spinorfield
	 * or a point-source spinorfield as trial-solution
	 */
	const Spinorfield solution(system);
	solution.cold();
	trace_squarenorm("\tinv. field before inversion ", solution);

	if(params.get_solver() == meta::Inputparameters::cg) {
		const QplusQminus fm(kappa, mubar, system);
		const int iterations = cg(&solution, fm, gf, phi, system, params.get_solver_prec());

		logger.debug() << "\t\t\tsolver solved in " << iterations << " iterations!";
		trace_squarenorm("\tinv. field after inversion ", solution);

		const Qminus qminus(kappa, mubar, system);
		qminus(phi_inv, gf, solution);

	} else  {
		const Qplus fm(kappa, mubar, system);
		const int iterations = bicgstab(&solution, fm, gf, phi, system, params.get_solver_prec());

		logger.debug() << "\t\t\tsolver solved in " << iterations << " iterations!";
		trace_squarenorm("\tinv. field after inversion ", solution);

		copyData(phi_inv, solution);
	}
	return squarenorm(*phi_inv);
}

hmc_float physics::algorithms::calc_s_fermion(const physics::lattices::Spinorfield_eo * const phi_inv, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
	using physics::lattices::Spinorfield_eo;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;

	auto params = system.get_inputparameters();

	logger.debug() << "\t\t\tstart solver";

	/** @todo at the moment, we can only put in a cold spinorfield
	 * or a point-source spinorfield as trial-solution
	 */
	const Spinorfield_eo solution(system);

	logger.debug() << "\t\t\tstart solver";

	//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
	if(params.get_solver() == meta::Inputparameters::cg) {
		solution.cold();
		trace_squarenorm("\tinv. field before inversion ", solution);

		const QplusQminus_eo fm(kappa, mubar, system);
		const int iterations = cg(&solution, fm, gf, phi, system, params.get_solver_prec());

		logger.debug() << "\t\t\tsolver solved in " << iterations << " iterations!";
		trace_squarenorm("\tinv. field after inversion ", solution);

		const Qminus_eo qminus(kappa, mubar, system);
		qminus(phi_inv, gf, solution);
	} else {
		solution.zero();
		solution.gamma5();
		trace_squarenorm("\tinv. field before inversion ", solution);
		trace_squarenorm("\tsource before inversion ", phi);

		const Qplus_eo fm(kappa, mubar, system);
		const int iterations = bicgstab(&solution, fm, gf, phi, system, params.get_solver_prec());

		logger.debug() << "\t\t\tsolver solved in " << iterations << " iterations!";
		trace_squarenorm("\tinv. field after inversion ", solution);

		copyData(phi_inv, solution);
	}
	return squarenorm(*phi_inv);
}

hmc_float physics::algorithms::calc_s_fermion_mp(const physics::lattices::Spinorfield * const phi_inv, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system)
{
	//this function essentially performs the same steps as in the non mass-prec case, however, one has to apply one more matrix multiplication
	//  therefore, comments are deleted here...
	//  Furthermore, in the bicgstab-case, the second inversions are not needed

	using physics::lattices::Spinorfield;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;

	logger.trace() << "\tHMC [DH]:\tcalc final fermion energy...";

	auto params = system.get_inputparameters();
	const hmc_float kappa = params.get_kappa();
	const hmc_float mubar = meta::get_mubar(params);

	const Spinorfield tmp(system);
	const Qplus qplus_mp(params.get_kappa_mp(), meta::get_mubar_mp(params), system);
	qplus_mp(&tmp, gf, phi);

	const Spinorfield solution(system);
	solution.cold();
	trace_squarenorm("\tinv. field before inversion ", solution);

	logger.debug() << "\t\t\tstart solver";

	if(params.get_solver() == meta::Inputparameters::cg) {
		const QplusQminus fm(kappa, mubar, system);
		const int iterations = cg(&solution, fm, gf, tmp, system, params.get_solver_prec());

		logger.debug() << "\t\t\tsolver solved in " << iterations << " iterations!";
		trace_squarenorm("\tinv. field after inversion ", solution);

		const Qminus qminus(kappa, mubar, system);
		qminus(phi_inv, gf, solution);
	} else  {
		const Qplus fm(kappa, mubar, system);
		const int iterations = bicgstab(&solution, fm, gf, tmp, system, params.get_solver_prec());

		logger.debug() << "\t\t\tsolver solved in " << iterations << " iterations!";
		trace_squarenorm("\tinv. field after inversion ", solution);

		copyData(phi_inv, solution);
	}
	return squarenorm(*phi_inv);
}

hmc_float physics::algorithms::calc_s_fermion_mp(const physics::lattices::Spinorfield_eo * const phi_inv, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system)
{
	//this function essentially performs the same steps as in the non mass-prec case, however, one has to apply one more matrix multiplication
	//  therefore, comments are deleted here...
	//  Furthermore, in the bicgstab-case, the second inversions are not needed
	using physics::lattices::Spinorfield_eo;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;

	logger.trace() << "\tHMC [DH]:\tcalc final fermion energy...";

	auto params = system.get_inputparameters();
	const hmc_float kappa = params.get_kappa();
	const hmc_float mubar = meta::get_mubar(params);

	logger.debug() << "\t\t\tstart solver";

	//sf_tmp = Qplus(light_mass) phi_mp
	const Spinorfield_eo tmp(system);
	const Qplus_eo qplus_mp(params.get_kappa_mp(), meta::get_mubar_mp(params), system);
	qplus_mp(&tmp, gf, phi);

	const Spinorfield_eo solution(system);

	logger.debug() << "\t\t\tstart solver";

	if(params.get_solver() == meta::Inputparameters::cg) {
		solution.cold();
		trace_squarenorm("\tinv. field before inversion ", solution);

		const QplusQminus_eo fm(kappa, mubar, system);
		const int iterations = cg(&solution, fm, gf, tmp, system, params.get_solver_prec());

		logger.debug() << "\t\t\tsolver solved in " << iterations << " iterations!";
		trace_squarenorm("\tinv. field after inversion ", solution);

		const Qminus_eo qminus(kappa, mubar, system);
		qminus(phi_inv, gf, solution);
	} else {
		solution.zero();
		solution.gamma5();
		trace_squarenorm("\tinv. field before inversion ", solution);
		trace_squarenorm("\tsource before inversion ", tmp);

		const Qplus_eo fm(kappa, mubar, system);
		const int iterations = bicgstab(&solution, fm, gf, tmp, system, params.get_solver_prec());

		logger.debug() << "\t\t\tsolver solved in " << iterations << " iterations!";
		trace_squarenorm("\tinv. field after inversion ", solution);

		copyData(phi_inv, solution);
	}
	return squarenorm(*phi_inv);
}

template <class SPINORFIELD> static hmc_observables metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf, const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p, const physics::lattices::Gaugemomenta& new_p, const SPINORFIELD& phi, const SPINORFIELD * const phi_inv, const hmc_float spinor_energy_init, const hmc_float spinor_energy_mp_init, const hardware::System& system)
{
	using namespace physics::algorithms;

	auto params = system.get_inputparameters();

	//Calc Hamiltonian
	logger.trace() << "\tHMC [DH]:\tCalculate Hamiltonian";
	hmc_float deltaH = 0.;
	hmc_float s_old = 0.;
	hmc_float s_new = 0.;

	//Gauge-Part
	hmc_float tplaq, splaq, plaq;
	hmc_float tplaq_new, splaq_new, plaq_new;
	hmc_float rect_new = 0.;
	hmc_float rect = 0.;
	hmc_complex poly;
	hmc_complex poly_new;
	//In this call, the observables are calculated already with appropiate Weighting factor of 2.0/(VOL4D*NDIM*(NDIM-1)*NC)
	gf.gaugeobservables(&plaq,  &tplaq, &splaq, &poly);
	new_u.gaugeobservables(&plaq_new,  &tplaq_new, &splaq_new, &poly_new);
	//plaq has to be divided by the norm-factor to get s_gauge
	hmc_float factor = 1. / (meta::get_plaq_norm(params));
	if(meta::get_use_rectangles(params) == true) {
		rect = gf.rectangles();
		rect_new = new_u.rectangles();
		hmc_float c0 = meta::get_c0(params);
		hmc_float c1 = meta::get_c1(params);
		deltaH = - beta * ( c0 * (plaq - plaq_new) / factor + c1 * ( rect - rect_new )  );
		s_old = - beta * ( c0 * (plaq) / factor + c1 * ( rect )  );
		s_new = - beta * ( c0 * (plaq_new) / factor + c1 * ( rect_new )  );
	} else {
		/** NOTE: the minus here is introduced to fit tmlqcd!!! */
		deltaH = -(plaq - plaq_new) * beta / factor;
		s_old = -(plaq ) * beta / factor;
		s_new = -(plaq_new) * beta / factor;
	}

	logger.debug() << "\tHMC [DH]:\tS[GF]_0:\t" << std::setprecision(10) << s_old;
	logger.debug() << "\tHMC [DH]:\tS[GF]_1:\t" << std::setprecision(10) << s_new;
	logger.info()  << "\tHMC [DH]:\tdS[GF]: \t" << std::setprecision(10) << deltaH;
	//check on NANs
	if (s_old != s_old || s_new != s_new || deltaH != deltaH) {
		throw Print_Error_Message("NAN occured in HMC! Aborting!", __FILE__, __LINE__);
	}

	//Gaugemomentum-Part
	hmc_float p2 = squarenorm(p);
	hmc_float new_p2 = squarenorm(new_p);
	//the energy is half the squarenorm
	deltaH += 0.5 * (p2 - new_p2);

	logger.debug() << "\tHMC [DH]:\tS[GM]_0:\t" << std::setprecision(10) << 0.5 * p2;
	logger.debug() << "\tHMC [DH]:\tS[GM]_1:\t" << std::setprecision(10) << 0.5 * new_p2;
	logger.info()  << "\tHMC [DH]:\tdS[GM]: \t" << std::setprecision(10) << 0.5 * (p2 - new_p2);
	//check on NANs
	if (p2 != p2 || new_p2 != new_p2 || deltaH != deltaH) {
		throw Print_Error_Message("NAN occured in HMC! Aborting!", __FILE__, __LINE__);
	}

	//Fermion-Part:
	if(! params.get_use_gauge_only() ) {
		if( params.get_use_mp() ) {
			//in this case one has contributions from det(m_light/m_heavy) and det(m_heavy)
			// det(m_heavy)
			hmc_float s_fermion_final;
			//initial energy has been computed in the beginning...
			s_fermion_final = calc_s_fermion(phi_inv, new_u, phi, system, params.get_kappa_mp(),  meta::get_mubar_mp(params));
			deltaH += spinor_energy_init - s_fermion_final;

			logger.debug() << "\tHMC [DH]:\tS[DET]_0:\t" << std::setprecision(10) <<  spinor_energy_init;
			logger.debug() << "\tHMC [DH]:\tS[DET]_1:\t" << std::setprecision(10) << s_fermion_final;
			logger.info() <<  "\tHMC [DH]:\tdS[DET]:\t" << spinor_energy_init - s_fermion_final;
			//check on NANs
			if (spinor_energy_init != spinor_energy_init || s_fermion_final != s_fermion_final || deltaH != deltaH) {
				throw Print_Error_Message("NAN occured in HMC! Aborting!", __FILE__, __LINE__);
			}

			// det(m_light/m_heavy)
			//initial energy has been computed in the beginning...
			hmc_float s_fermion_mp_final = calc_s_fermion_mp(phi_inv, new_u, phi, system);
			deltaH += spinor_energy_mp_init - s_fermion_mp_final;

			logger.debug() << "\tHMC [DH]:\tS[DETRAT]_0:\t" << std::setprecision(10) <<  spinor_energy_mp_init;
			logger.debug() << "\tHMC [DH]:\tS[DETRAT]_1:\t" << std::setprecision(10) << s_fermion_mp_final;
			logger.info() <<  "\tHMC [DH]:\tdS[DETRAT]:\t" << spinor_energy_mp_init - s_fermion_mp_final;
			//check on NANs
			if (spinor_energy_mp_init != spinor_energy_mp_init || s_fermion_mp_final != s_fermion_mp_final || deltaH != deltaH) {
				throw Print_Error_Message("NAN occured in HMC! Aborting!", __FILE__, __LINE__);
			}
		} else {
			hmc_float s_fermion_final = calc_s_fermion(phi_inv, new_u, phi, system);
			deltaH += spinor_energy_init - s_fermion_final;

			logger.debug() << "\tHMC [DH]:\tS[DET]_0:\t" << std::setprecision(10) <<  spinor_energy_init;
			logger.debug() << "\tHMC [DH]:\tS[DET]_1:\t" << std::setprecision(10) << s_fermion_final;
			logger.info() <<  "\tHMC [DH]:\tdS[DET]: \t" << spinor_energy_init - s_fermion_final;
			//check on NANs
			if (spinor_energy_init != spinor_energy_init || s_fermion_final != s_fermion_final || deltaH != deltaH) {
				throw Print_Error_Message("NAN occured in HMC! Aborting!", __FILE__, __LINE__);
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
	logger.info() << "\tHMC [DH]:\tdH:\t\t" << deltaH;
	logger.info() << "\tHMC [MET]:\tAcc-Prop:\t" << compare_prob;
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

hmc_observables physics::algorithms::metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf, const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p, const physics::lattices::Gaugemomenta& new_p, const physics::lattices::Spinorfield& phi, const physics::lattices::Spinorfield * const phi_inv, const hmc_float spinor_energy_init, const hmc_float spinor_energy_mp_init, const hardware::System& system)
{
	return ::metropolis(rnd, beta, gf, new_u, p, new_p, phi, phi_inv, spinor_energy_init, spinor_energy_mp_init, system);
}

hmc_observables physics::algorithms::metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf, const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p, const physics::lattices::Gaugemomenta& new_p, const physics::lattices::Spinorfield_eo& phi, const physics::lattices::Spinorfield_eo * const phi_inv, const hmc_float spinor_energy_init, const hmc_float spinor_energy_mp_init, const hardware::System& system)
{
	return ::metropolis(rnd, beta, gf, new_u, p, new_p, phi, phi_inv, spinor_energy_init, spinor_energy_mp_init, system);
}
