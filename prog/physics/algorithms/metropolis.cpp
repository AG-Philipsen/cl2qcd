/** @file
 * Implementation of the metropolis algorithm
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * (c) 2012-2013 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
 */

#include "metropolis.hpp"

#include "solver.hpp"
#include "../lattices/util.hpp"

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
	solution.cold();

	logger.debug() << "\t\t\tstart solver";

	//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
	if(params.get_solver() == meta::Inputparameters::cg) {
		trace_squarenorm("\tinv. field before inversion ", solution);

		const QplusQminus_eo fm(kappa, mubar, system);
		const int iterations = cg(&solution, fm, gf, phi, system, params.get_solver_prec());

		logger.debug() << "\t\t\tsolver solved in " << iterations << " iterations!";
		trace_squarenorm("\tinv. field after inversion ", solution);

		const Qminus_eo qminus(kappa, mubar, system);
		qminus(phi_inv, gf, solution);
	} else {
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
