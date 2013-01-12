/** @file
 * Implementation of the inversion algorithms
 */

#include "inversion.hpp"
#include "../meta/util.hpp"
#include "solver.hpp"
#include <cassert>
#include "../lattices/util.hpp"

static void invert_M_nf2_upperflavour(const physics::lattices::Spinorfield* result, const physics::lattices::Gaugefield& gaugefield, const physics::lattices::Spinorfield* source, const hardware::System& system);

void physics::algorithms::perform_inversion(const std::vector<const physics::lattices::Spinorfield*> * result, physics::lattices::Gaugefield* gaugefield, const std::vector<const physics::lattices::Spinorfield*>& sources, const hardware::System& system)
{
	int num_sources = sources.size();
	auto params = system.get_inputparameters();

	//apply stout smearing if wanted
	if(params.get_use_smearing())
		gaugefield->smear();

	for(int k = 0; k < num_sources; k++) {
		logger.debug() << "calling solver..";
		invert_M_nf2_upperflavour(result->at(k), *gaugefield, sources[k], system);
	}

	if(params.get_use_smearing())
		gaugefield->unsmear();
}

static void invert_M_nf2_upperflavour(const physics::lattices::Spinorfield* result, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield* source, const hardware::System& system)
{
	using namespace physics::lattices;
	using namespace physics::algorithms::solvers;
	using namespace physics::fermionmatrix;

	/** This solves the sparse-matrix system
	 *  A x = b
	 *  with  x == result
	 *        A == gaugefield
	 *        b == source
	 * using a Krylov-Solver (BiCGStab or CG)
	 */

	// assert a single GPU
	assert(result->get_buffers().size() == 1);

	auto params = system.get_inputparameters();

	int converged;

	if(!params.get_use_eo()) {
		//noneo case
		//Trial solution
		///@todo this should go into a more general function
		result->cold();
		if(params.get_solver() == meta::Inputparameters::cg) {
			Spinorfield tmp(system);
			//to use cg, one needs an hermitian matrix, which is QplusQminus
			//the source must now be gamma5 b, to obtain the desired solution in the end
			copyData(&tmp, source);
			tmp.gamma5();
			QplusQminus f_neo(params.get_kappa(), meta::get_mubar(params), system);
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(result_buf, false, "\tinv. field before inversion ");
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(source_buf, false, "\tsource before inversion ");
			converged = cg(result, f_neo, gf, tmp, system, params.get_solver_prec());
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(result_buf, false, "\tinv. field after inversion ");
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(source_buf, false, "\tsource after inversion ");
			copyData(&tmp, result);
			//now, calc Qminus result_buf to obtain x = A^⁻1 b
			Qminus qminus(params.get_kappa(), meta::get_mubar(params), system);
			qminus(result, gf, tmp);
		} else {
			M f_neo(params.get_kappa(), meta::get_mubar(params), system);
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(result_buf, false, "\tinv. field before inversion ");
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(source_buf, false, "\tsource before inversion ");
			converged = bicgstab(result, f_neo, gf, *source, system, params.get_solver_prec());
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(result_buf, false, "\tinv. field after inversion ");
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(source_buf, false, "\tsource after inversion ");
		}
	} else {
		/**
		 * If even-odd-preconditioning is used, the inversion is split up
		 * into even and odd parts using Schur decomposition, assigning the
		 * non-trivial inversion to the even sites (see DeGran/DeTar p 174ff).
		 */
		//init some helping buffers
		const Spinorfield_eo source_even(system);
		const Spinorfield_eo source_odd(system);
		const Spinorfield_eo tmp1(system);
		const Spinorfield_eo tmp2(system);
		const Spinorfield_eo result_eo(system);
		hmc_complex one = hmc_complex_one;
		hmc_complex mone = { -1., 0.};

		//convert source and input-vector to eoprec-format
		/**
		 * This currently is a workaround connceted to issue #387
		 * the roles of even/odd vectors are interchanged!
		 * @todo: fix
		 * original code:
		 * spinor_code->convert_to_eoprec_device(&clmem_source_even, &clmem_source_odd, source_buf);
		 * workaround:
		 */
		convert_to_eoprec(&source_odd, &source_even, *source);

		// TODO readd if(logger.beDebug()) solver->print_info_inv_field(source_buf, false, "\tsource before inversion ");
		// TODO readd if(logger.beDebug()) solver->print_info_inv_field(&clmem_source_even, true, "\teven source before inversion ");
		// TODO readd if(logger.beDebug()) solver->print_info_inv_field(&clmem_source_odd, true, "\todd source before inversion ");

		//prepare sources
		/**
		 * This changes the even source according to (with A = M + D):
		 *  b_e = b_e - D_eo M_inv b_o
		 */
		if(params.get_fermact() == meta::Inputparameters::wilson) {
			//in this case, the diagonal matrix is just 1 and falls away.
			dslash(&tmp1, gf, source_odd, EVEN);
			saxpy(&source_even, one, source_even, tmp1);
		} else if(params.get_fermact() == meta::Inputparameters::twistedmass) {
			M_tm_inverse_sitediagonal(&tmp1, source_odd);
			dslash(&tmp2, gf, tmp1, EVEN);
			saxpy(&source_even, one, source_even, tmp2);
		}

		//Trial solution
		///@todo this should go into a more general function
		result_eo.cold();
		logger.debug() << "start eoprec-inversion";
		//even solution
		if(params.get_solver() == meta::Inputparameters::cg) {
			//to use cg, one needs an hermitian matrix, which is QplusQminus
			//the source must now be gamma5 b, to obtain the desired solution in the end
			source_even.gamma5();
			QplusQminus_eo f_eo(params.get_kappa(), meta::get_mubar(params), system);
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(&result_buf_eo, true, "\tinv field before inversion ");
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(&clmem_source_even, true, "\tsource before inversion ");
			converged = cg(&result_eo, f_eo, gf, source_even, system, params.get_solver_prec());
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(&result_buf_eo, true, "\tinv field after inversion ");
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(&clmem_source_even, true, "\tsource after inversion ");
			//now, calc Qminus result_buf_eo to obtain x = A^⁻1 b
			//therefore, use source as an intermediate buffer
			Qminus_eo qminus(params.get_kappa(), meta::get_mubar(params), system);
			qminus(&source_even, gf, result_eo);
			//save the result to result_buf
			copyData(&result_eo, source_even);
		} else {
			Aee f_eo(params.get_kappa(), meta::get_mubar(params), system);
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(&result_buf_eo, true, "\tinv field before inversion ");
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(&clmem_source_even, true, "\tsource before inversion ");
			converged = bicgstab(&result_eo, f_eo, gf, source_even, system, params.get_solver_prec());
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(&result_buf_eo, true, "\tinv field after inversion ");
			// TODO readd if(logger.beDebug()) solver->print_info_inv_field(&clmem_source_even, true, "\tsource after inversion ");
		}

		//odd solution
		/** The odd solution is obtained from the even one according to:
		 *  x_o = M_inv D x_e - M_inv b_o
		 * @todo: find out why it must be
		 *  x_o = -(M_inv D x_e + M_inv b_o)
		 */
		if(params.get_fermact() == meta::Inputparameters::wilson) {
			//in this case, the diagonal matrix is just 1 and falls away.
			dslash(&tmp1, gf, result_eo, ODD);
			saxpy(&tmp1, mone, tmp1, source_odd);
			sax(&tmp1, mone, tmp1);
		} else if(params.get_fermact() == meta::Inputparameters::twistedmass) {
			dslash(&tmp2, gf, result_eo, ODD);
			M_tm_inverse_sitediagonal(&tmp1, tmp2);
			M_tm_inverse_sitediagonal(&tmp2, source_odd);
			saxpy(&tmp1, mone, tmp1, tmp2);
			sax(&tmp1, mone, tmp1);
		}

		///CP: whole solution
		//convert source and input-vector to eoprec-format
		/**
		 * This currently is a workaround connceted to issue #387
		 * the roles of even/odd vectors are interchanged!
		 * @todo: fix
		 * original code:
		 * //CP: suppose the even sol is saved in inout_eoprec, the odd one in clmem_tmp_eo_1
		 * spinor_code->convert_from_eoprec_device(&result_buf_eo, &clmem_tmp_eo_1, result_buf);
		 * workaround:
		 */
		//CP: suppose the odd sol is saved in inout_eoprec, the even one in clmem_tmp_eo_1
		convert_from_eoprec(result, tmp1, result_eo);
	}

	// TODO readd if(logger.beDebug()) solver->print_info_inv_field(result_buf, false, "\tsolution ");

	// TODO catch or document exceptions
	logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
}

