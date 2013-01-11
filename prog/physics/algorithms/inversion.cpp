/** @file
 * Implementation of the inversion algorithms
 */

#include "inversion.hpp"
#include "../meta/util.hpp"
#include <cassert>

static void invert_M_nf2_upperflavour(const physics::lattices::Spinorfield* result, const physics::lattices::Gaugefield& gaugefield, const physics::lattices::Spinorfield* source, const meta::Inputparameters& params);

void physics::algorithms::perform_inversion(const std::vector<const physics::lattices::Spinorfield*> * result, physics::lattices::Gaugefield* gaugefield, const std::vector<const physics::lattices::Spinorfield*>& sources, const meta::Inputparameters& params)
{
	int num_sources = params.get_num_sources();

	//apply stout smearing if wanted
	if(params.get_use_smearing())
		gaugefield->smear();

	for(int k = 0; k < num_sources; k++) {
		logger.debug() << "calling solver..";
		invert_M_nf2_upperflavour(result->at(k), *gaugefield, sources.at(k), params);
	}

	if(params.get_use_smearing())
		gaugefield->unsmear();
}

static void invert_M_nf2_upperflavour(const physics::lattices::Spinorfield* result, const physics::lattices::Gaugefield& gaugefield, const physics::lattices::Spinorfield* source, const meta::Inputparameters& params)
{
	/** This solves the sparse-matrix system
	 *  A x = b
	 *  with  x == result
	 *        A == gaugefield
	 *        b == source
	 * using a Krylov-Solver (BiCGStab or CG)
	 */

	// assert a single GPU
	assert(result->get_buffers().size() == 1);

	auto result_buf = result->get_buffers().at(0);
	auto source_buf = source->get_buffers().at(0);
	auto gf_buf = gaugefield.get_buffers().at(0);
	auto device = result_buf->get_device();

	int converged = -1;
	hardware::code::Fermions * solver = device->get_fermion_code();
	auto spinor_code = device->get_spinor_code();

	if(!params.get_use_eo()) {
		//noneo case
		//Trial solution
		///@todo this should go into a more general function
		spinor_code->set_spinorfield_cold_device(result_buf);
		if(params.get_solver() == meta::Inputparameters::cg) {
			//to use cg, one needs an hermitian matrix, which is QplusQminus
			//the source must now be gamma5 b, to obtain the desired solution in the end
			solver->gamma5_device(source_buf);
			hardware::code::QplusQminus f_neo(solver);
			const hardware::buffers::Plain<spinor> clmem_tmp  (meta::get_spinorfieldsize(params), device);
			converged = solver->cg(f_neo, result_buf, source_buf, gf_buf, params.get_solver_prec());
			hardware::buffers::copyData(&clmem_tmp, result_buf);
			//now, calc Qminus result_buf to obtain x = A^⁻1 b
			solver->Qminus(&clmem_tmp, result_buf, gf_buf, params.get_kappa(), meta::get_mubar(params ));
		} else {
			hardware::code::M f_neo(solver);
			converged = solver->bicgstab(f_neo, result_buf, source_buf, gf_buf, params.get_solver_prec());
		}
	} else {
		/**
		 * If even-odd-preconditioning is used, the inversion is split up
		 * into even and odd parts using Schur decomposition, assigning the
		 * non-trivial inversion to the even sites (see DeGran/DeTar p 174ff).
		 */
		//init some helping buffers
		const hardware::buffers::Spinor clmem_source_even  (meta::get_eoprec_spinorfieldsize(params), device);
		const hardware::buffers::Spinor clmem_source_odd  (meta::get_eoprec_spinorfieldsize(params), device);
		const hardware::buffers::Spinor clmem_tmp_eo_1  (meta::get_eoprec_spinorfieldsize(params), device);
		const hardware::buffers::Spinor clmem_tmp_eo_2  (meta::get_eoprec_spinorfieldsize(params), device);
		const hardware::buffers::Plain<hmc_complex> clmem_one (1, device);
		const hardware::buffers::Spinor result_buf_eo  (meta::get_eoprec_spinorfieldsize(params), device);
		const hardware::buffers::Plain<hmc_complex> clmem_mone (1, solver->get_device());
		hmc_complex one = hmc_complex_one;
		hmc_complex mone = {-1.,0.};
		clmem_one.load(&one);
		clmem_mone.load(&mone);

		//convert source and input-vector to eoprec-format
		/**
		 * This currently is a workaround connceted to issue #387
		 * the roles of even/odd vectors are interchanged!
		 * @todo: fix
		 * original code:
		 * spinor_code->convert_to_eoprec_device(&clmem_source_even, &clmem_source_odd, source_buf);
		 * workaround:
		 */
		spinor_code->convert_to_eoprec_device(&clmem_source_odd, &clmem_source_even, source_buf);

		//prepare sources
		/**
		 * This changes the even source according to (with A = M + D):
		 *  b_e = b_e - D_eo M_inv b_o
		 */
		if(params.get_fermact() == meta::Inputparameters::wilson) {
			//in this case, the diagonal matrix is just 1 and falls away.
			solver->dslash_eo_device(&clmem_source_odd, &clmem_tmp_eo_1, gf_buf, EVEN);
			spinor_code->saxpy_eoprec_device(&clmem_source_even, &clmem_tmp_eo_1, &clmem_one, &clmem_source_even);
		} else if(params.get_fermact() == meta::Inputparameters::twistedmass) {
			solver->M_tm_inverse_sitediagonal_device(&clmem_source_odd, &clmem_tmp_eo_1);
			solver->dslash_eo_device(&clmem_tmp_eo_1, &clmem_tmp_eo_2, gf_buf, EVEN);
			spinor_code->saxpy_eoprec_device(&clmem_source_even, &clmem_tmp_eo_2, &clmem_one, &clmem_source_even);
		}

		//Trial solution
		///@todo this should go into a more general function
		spinor_code->set_eoprec_spinorfield_cold_device(&result_buf_eo);
		logger.debug() << "start eoprec-inversion";
		//even solution
		if(params.get_solver() == meta::Inputparameters::cg) {
			//to use cg, one needs an hermitian matrix, which is QplusQminus
			//the source must now be gamma5 b, to obtain the desired solution in the end
			solver->gamma5_eo_device(&clmem_source_even);
			hardware::code::QplusQminus_eo f_eo(solver);
			converged = solver->cg_eo(f_eo, &result_buf_eo, &clmem_source_even, gf_buf, params.get_solver_prec());
			//now, calc Qminus result_buf_eo to obtain x = A^⁻1 b
			//therefore, use source as an intermediate buffer
			solver->Qminus_eo(&result_buf_eo, &clmem_source_even, gf_buf, params.get_kappa(), meta::get_mubar(params));
			//save the result to result_buf
			hardware::buffers::copyData(&result_buf_eo, &clmem_source_even);
		} else {
			hardware::code::Aee f_eo(solver);
			converged = solver->bicgstab_eo(f_eo, &result_buf_eo, &clmem_source_even, gf_buf, params.get_solver_prec());
		}

		//odd solution
		/** The odd solution is obtained from the even one according to:
		 *  x_o = M_inv D x_e - M_inv b_o
		 * @todo: find out why it must be
		 *  x_o = -(M_inv D x_e + M_inv b_o)
		 */
		if(params.get_fermact() == meta::Inputparameters::wilson) {
			//in this case, the diagonal matrix is just 1 and falls away.
			solver->dslash_eo_device(&result_buf_eo, &clmem_tmp_eo_1, gf_buf, ODD);
			spinor_code->saxpy_eoprec_device(&clmem_tmp_eo_1, &clmem_source_odd, &clmem_mone, &clmem_tmp_eo_1);
			spinor_code->sax_eoprec_device(&clmem_tmp_eo_1, &clmem_mone, &clmem_tmp_eo_1);
		} else if(params.get_fermact() == meta::Inputparameters::twistedmass) {
			solver->dslash_eo_device(&result_buf_eo, &clmem_tmp_eo_2, gf_buf, ODD);
			solver->M_tm_inverse_sitediagonal_device(&clmem_tmp_eo_2, &clmem_tmp_eo_1);
			solver->M_tm_inverse_sitediagonal_device(&clmem_source_odd, &clmem_tmp_eo_2);
			spinor_code->saxpy_eoprec_device(&clmem_tmp_eo_1, &clmem_tmp_eo_2, &clmem_mone, &clmem_tmp_eo_1);
			spinor_code->sax_eoprec_device(&clmem_tmp_eo_1, &clmem_mone, &clmem_tmp_eo_1);
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
		spinor_code->convert_from_eoprec_device( &clmem_tmp_eo_1, &result_buf_eo, result_buf);
	}

	if (converged < 0) {
		if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
		else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
	} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
}

