/** @file
 * Implementation of the inversion algorithms
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#include "inversion.hpp"
#include "../../meta/util.hpp"
#include "solvers/solvers.hpp"
#include <cassert>
#include "../lattices/util.hpp"
#include "../lattices/swappable.hpp"

static void invert_M_nf2_upperflavour(const physics::lattices::Spinorfield* result, const physics::lattices::Gaugefield& gaugefield, const physics::lattices::Spinorfield* source, const hardware::System& system);

template<class Spinorfield> static hmc_float print_debug_inv_field(const Spinorfield& in, std::string msg);
template<class Spinorfield> static hmc_float print_debug_inv_field(const Spinorfield* in, std::string msg);

void physics::algorithms::perform_inversion(const std::vector<physics::lattices::Spinorfield*> * result, const physics::lattices::Gaugefield* gaugefield, const std::vector<physics::lattices::Spinorfield*>& sources, const hardware::System& system)
{
	const auto & params = system.get_inputparameters();
	const size_t num_sources = sources.size();

	//apply stout smearing if wanted
	if(params.get_use_smearing())
		gaugefield->smear();

	for(size_t k = 0; k < num_sources; k++) {
		logger.debug() << "calling solver..";

		auto source = sources[k];
		auto res = result->at(k);

		try_swap_in(source);
		try_swap_in(res);

		invert_M_nf2_upperflavour(res, *gaugefield, source, system);

		try_swap_out(source);
		try_swap_out(res);
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

	const auto & params = system.get_inputparameters();

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
			converged = cg(result, f_neo, gf, tmp, system, params.get_solver_prec());
			copyData(&tmp, result);
			//now, calc Qminus result_buf to obtain x = A^⁻1 b
			Qminus qminus(params.get_kappa(), meta::get_mubar(params), system);
			qminus(result, gf, tmp);
		} else {
			M f_neo(params.get_kappa(), meta::get_mubar(params), system);
			converged = bicgstab(result, f_neo, gf, *source, system, params.get_solver_prec());
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
		 * This currently is a workaround connected to issue #387
		 * the roles of even/odd vectors are interchanged!
		 * @todo: fix
		 * original code:
		 * spinor_code->convert_to_eoprec_device(&clmem_source_even, &clmem_source_odd, source_buf);
		 * workaround:
		 */
		convert_to_eoprec(&source_odd, &source_even, *source);

		///@todo: work over these debug outputs once the above issue is settled
		print_debug_inv_field(source, "\tsource before inversion ");
		print_debug_inv_field(&source_even, "\teven source before inversion ");
		print_debug_inv_field(&source_odd, "\todd source before inversion ");

		//prepare sources
		/**
		 * This changes the even source according to (with A = M + D):
		 *  b_e = b_e - D_eo M_inv b_o
		 */
		if(params.get_fermact() == meta::action::wilson) {
			//in this case, the diagonal matrix is just 1 and falls away.
			dslash(&tmp1, gf, source_odd, EVEN);
			saxpy(&source_even, one, source_even, tmp1);
		} else if(params.get_fermact() == meta::action::twistedmass) {
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
			converged = cg(&result_eo, f_eo, gf, source_even, system, params.get_solver_prec());
			//now, calc Qminus result_buf_eo to obtain x = A^⁻1 b
			//therefore, use source as an intermediate buffer
			Qminus_eo qminus(params.get_kappa(), meta::get_mubar(params), system);
			qminus(&source_even, gf, result_eo);
			//save the result to result_buf
			copyData(&result_eo, source_even);
		} else {
			Aee f_eo(params.get_kappa(), meta::get_mubar(params), system);
			converged = bicgstab(&result_eo, f_eo, gf, source_even, system, params.get_solver_prec());
		}

		//odd solution
		/** The odd solution is obtained from the even one according to:
		 *  x_o = M_inv b_o - M_inv D x_e  
		 * @todo: find out why it must be (issue #389)
		 *  x_o = - M_inv b_o - M_inv D x_e
		 *      = -(M_inv D x_e + M_inv b_o)
		 */
		if(params.get_fermact() == meta::action::wilson) {
			//in this case, the diagonal matrix is just 1 and falls away.
			dslash(&tmp1, gf, result_eo, ODD);
			saxpy(&tmp1, mone, tmp1, source_odd);
			sax(&tmp1, mone, tmp1);
		} else if(params.get_fermact() == meta::action::twistedmass) {
			dslash(&tmp2, gf, result_eo, ODD);
			M_tm_inverse_sitediagonal(&tmp1, tmp2);
			M_tm_inverse_sitediagonal(&tmp2, source_odd);
			saxpy(&tmp1, mone, tmp1, tmp2);
			sax(&tmp1, mone, tmp1);
		}

		///CP: whole solution
		//convert source and input-vector to eoprec-format
		/**
		 * This currently is a workaround connected to issue #387
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

	print_debug_inv_field(result, "\tsolution ");

	// TODO catch or document exceptions
	logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
}

template<class Spinorfield> static hmc_float print_debug_inv_field(const Spinorfield* in, std::string msg)
{
	return print_debug_inv_field(*in, msg);
}
template<class Spinorfield> static hmc_float print_debug_inv_field(const Spinorfield& in, std::string msg)
{
	if(logger.beDebug()) {
		hmc_float tmp = squarenorm(in);
		logger.debug() << msg << std::scientific << std::setprecision(10) << tmp;
		return tmp;
	}
	return 0;
}

