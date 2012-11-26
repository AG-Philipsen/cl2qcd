#include "hmc.hpp"

#include "../../logger.hpp"
#include "../../meta/util.hpp"
#include "../device.hpp"
#include <fstream>

using namespace std;

hardware::code::Hmc::Hmc(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module(params, device),
	  clmem_s_fermion_init(1, device),
	  clmem_s_fermion_mp_init(1, device),
	  clmem_p2(1, device),
	  clmem_new_p2(1, device),
	  clmem_s_fermion(1, device),
	  clmem_p(meta::get_vol4d(params) * NDIM, device),
	  clmem_new_p(meta::get_vol4d(params) * NDIM, device),
	  new_u(device->get_gaugefield_code()->get_gaugefield()->get_elements(), device),
	  clmem_force(meta::get_vol4d(params) * NDIM, device),
	  clmem_phi_inv(meta::get_spinorfieldsize(params), device),
	  clmem_phi_inv_eo(meta::get_eoprec_spinorfieldsize(params), device),
	  clmem_phi(meta::get_spinorfieldsize(params), device),
	  clmem_phi_mp(meta::get_spinorfieldsize(params), device),
	  clmem_phi_eo(meta::get_eoprec_spinorfieldsize(params), device),
	  clmem_phi_mp_eo(meta::get_eoprec_spinorfieldsize(params), device)
{
}

hardware::code::Hmc::~Hmc()
{
}

void hardware::code::Hmc::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
}

size_t hardware::code::Hmc::get_read_write_size(const std::string& in) const
{
}

uint64_t hardware::code::Hmc::get_flop_size(const std::string& in) const
{
}

void hardware::code::Hmc::print_profiling(const std::string& filename, int number) const
{
}

////////////////////////////////////////////////////
//Access to members

const hardware::buffers::Gaugemomentum * hardware::code::Hmc::get_clmem_p()
{
	return &clmem_p;
}

const hardware::buffers::Gaugemomentum * hardware::code::Hmc::get_clmem_new_p()
{
	return &clmem_new_p;
}

const hardware::buffers::SU3 * hardware::code::Hmc::get_new_u()
{
	return &new_u;
}

const hardware::buffers::Plain<spinor> * hardware::code::Hmc::get_clmem_phi()
{
	return &clmem_phi;
}

const hardware::buffers::Plain<spinor> * hardware::code::Hmc::get_clmem_phi_inv()
{
	return &clmem_phi_inv;
}

const hardware::buffers::Spinor * hardware::code::Hmc::get_clmem_phi_eo()
{
	return &clmem_phi_eo;
}

const hardware::buffers::Spinor * hardware::code::Hmc::get_clmem_phi_inv_eo()
{
	return &clmem_phi_inv_eo;
}

const hardware::buffers::Plain<spinor> * hardware::code::Hmc::get_clmem_phi_mp()
{
	return &clmem_phi_mp;
}

const hardware::buffers::Spinor * hardware::code::Hmc::get_clmem_phi_mp_eo()
{
	return &clmem_phi_mp_eo;
}

hardware::buffers::Plain<hmc_float> * hardware::code::Hmc::get_clmem_s_fermion_init()
{
	return &clmem_s_fermion_init;
}

hardware::buffers::Plain<hmc_float> * hardware::code::Hmc::get_clmem_s_fermion_mp_init()
{
	return &clmem_s_fermion_mp_init;
}

////////////////////////////////////////////////////
//Methods needed for the HMC-algorithm
void hardware::code::Hmc::generate_spinorfield_gaussian(const hardware::buffers::PRNGBuffer * prng)
{
	auto spinor_code = get_device()->get_spinor_code();
	if(get_parameters().get_use_eo() == true) {
		spinor_code->generate_gaussian_spinorfield_eo_device(get_clmem_phi_inv_eo(), prng);
	} else {
		spinor_code->generate_gaussian_spinorfield_device(get_clmem_phi_inv(), prng);
	}
}

void hardware::code::Hmc::md_update_spinorfield(const hardware::buffers::SU3 * gaugefield, hmc_float kappa, hmc_float mubar)
{
	auto fermion_code = get_device()->get_fermion_code();

	//suppose the initial gaussian field is saved in clmem_phi_inv (see above).
	//  then the "phi" = Dpsi from the algorithm is stored in clmem_phi
	//  which then has to be the source of the inversion
	if(get_parameters().get_use_eo() == true) {
		fermion_code->Qplus_eo (&clmem_phi_inv_eo, &clmem_phi_eo , gaugefield, kappa, mubar);
		if(logger.beDebug()) fermion_code->print_info_inv_field(&clmem_phi_eo, true, "\tinit field after update ");
	} else {
		fermion_code->Qplus(&clmem_phi_inv, &clmem_phi , gaugefield, kappa, mubar);
		if(logger.beDebug()) fermion_code->print_info_inv_field(&clmem_phi, false, "\tinit field after update ");
	}
}

void hardware::code::Hmc::md_update_spinorfield_mp(usetimer * solvertimer, const hardware::buffers::SU3 * gaugefield)
{
	using namespace hardware::buffers;

	auto fermion_code = get_device()->get_fermion_code();
	auto spinor_code = get_device()->get_spinor_code();

	///@todo solvertimer is not used here yet...
	//suppose the initial gaussian field is saved in clmem_phi_inv (see above).
	//  then the "phi" = Dpsi from the algorithm is stored in clmem_phi
	//  which then has to be the source of the inversion
	//in the mass preconditioning case, this is a bit more complicated and involves an inversion
	if(get_parameters().get_use_eo() == true) {
		//CP: Init tmp spinorfield
		Spinor sf_eo_tmp(clmem_phi_inv_eo.get_elements(), get_device());

		//sf_eo_tmp = Qplus_eo(light_mass) phi_inv_eo
		fermion_code->Qplus_eo (&clmem_phi_inv_eo, &sf_eo_tmp , gaugefield);

		//Now one needs ( Qplus_eo )^-1 (heavy_mass) using sf_eo_tmp as source to get phi_mp_eo
		///@todo Do not always use bicgstab here
		logger.debug() << "\t\t\tstart solver";

		/** @todo at the moment, we can only put in a cold spinorfield
		 * or a point-source spinorfield as trial-solution
		 */
		spinor_code->set_zero_spinorfield_eoprec_device(get_clmem_phi_mp_eo());

		fermion_code->gamma5_eo_device(get_clmem_phi_mp_eo());

		int converged = -1;
		if(logger.beDebug()) fermion_code->print_info_inv_field(get_clmem_phi_mp_eo(), true, "\tinv. field before inversion ");
		if(logger.beDebug()) fermion_code->print_info_inv_field(&sf_eo_tmp, true, "\tsource before inversion ");
		converged = fermion_code->bicgstab_eo(hardware::code::Qplus_eo(fermion_code), this->get_clmem_phi_mp_eo(), &sf_eo_tmp, &new_u, get_parameters().get_solver_prec(), get_parameters().get_kappa_mp(), meta::get_mubar_mp(get_parameters()));
		if (converged < 0) {
			if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
			else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
		} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
		if(logger.beDebug()) fermion_code->print_info_inv_field(get_clmem_phi_mp_eo(), true, "\tinv. field after inversion ");
		if(logger.beDebug()) fermion_code->print_info_inv_field(&clmem_phi_mp_eo, true, "\tinit field after update ");
	} else {
		//CP: Init tmp spinorfield
		Plain<spinor> sf_tmp(clmem_phi_inv.get_elements(), get_device());

		//sf_tmp = Qplus(light_mass) phi_inv
		fermion_code->Qplus (&clmem_phi_inv, &sf_tmp , gaugefield);

		//Now one needs ( Qplus )^-1 (heavy_mass) using sf_tmp as source to get phi_mp
		//use always bicgstab here
		logger.debug() << "\t\t\tstart solver";

		/** @todo at the moment, we can only put in a cold spinorfield
		 * or a point-source spinorfield as trial-solution
		 */
		spinor_code->set_zero_spinorfield_device(get_clmem_phi_mp());
		fermion_code->gamma5_device(get_clmem_phi_mp());

		int converged = -1;
		if(logger.beDebug()) fermion_code->print_info_inv_field(get_clmem_phi_mp(), false, "\tinv. field before inversion ");
		if(logger.beDebug()) fermion_code->print_info_inv_field(&sf_tmp, false, "\tsource before inversion ");
		converged = fermion_code->bicgstab(hardware::code::Qplus(fermion_code), this->get_clmem_phi_mp(), &sf_tmp, &new_u, get_parameters().get_solver_prec(), get_parameters().get_kappa_mp(), get_mubar_mp(get_parameters()));
		if (converged < 0) {
			if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
			else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
		} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
		if(logger.beDebug()) fermion_code->print_info_inv_field(get_clmem_phi_mp(), false, "\tinv. field after inversion ");
		if(logger.beDebug()) fermion_code->print_info_inv_field(&clmem_phi_mp, false, "\tinit field after update ");
	}
}

//this function takes to args kappa and mubar because one has to use it with different masses when mass-prec is used and when not
void hardware::code::Hmc::calc_fermion_force(usetimer * solvertimer, hmc_float kappa, hmc_float mubar)
{
	logger.debug() << "\t\tcalc fermion_force...";
	auto fermion_code = get_device()->get_fermion_code();
	auto spinor_code = get_device()->get_spinor_code();
	int converged = -1;
	if(get_parameters().get_use_eo() == true) {
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
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
			spinor_code->set_eoprec_spinorfield_cold_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field before inversion ");
			converged = fermion_code->cg_eo(hardware::code::QplusQminus_eo(fermion_code), fermion_code->get_inout_eo(), this->get_clmem_phi_eo(), &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field after inversion ");

			/**
			 * Y_even is now just
			 *  Y_even = (Qminus_eo) X_even = (Qminus_eo) (Qplusminus_eo)^-1 psi =
			 *    = (Qplus_eo)^-1 psi
			 */
			fermion_code->Qminus_eo(fermion_code->get_inout_eo(), &clmem_phi_inv_eo, &new_u, kappa, mubar);
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
			spinor_code->set_zero_spinorfield_eoprec_device(fermion_code->get_inout_eo());
			fermion_code->gamma5_eo_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field before inversion ");
			if(logger.beDebug()) fermion_code->print_info_inv_field(get_clmem_phi_eo(), true, "\t\t\tsource before inversion ");
			converged = fermion_code->bicgstab_eo(hardware::code::Qplus_eo(fermion_code), fermion_code->get_inout_eo(), this->get_clmem_phi_eo(), &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv_eo, fermion_code->get_inout_eo());

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
			hardware::buffers::copyData(fermion_code->get_source_even(), fermion_code->get_inout_eo());

			//this sets clmem_inout cold as trial-solution
			spinor_code->set_eoprec_spinorfield_cold_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field before inversion ");
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_source_even(), true, "\t\t\tsource before inversion ");
			converged = fermion_code->bicgstab_eo(hardware::code::Qminus_eo(fermion_code), fermion_code->get_inout_eo(), fermion_code->get_source_even(), &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field after inversion ");
		}
		/**
		 * At this point, one has calculated X_odd and Y_odd.
		 * If one has a fermionmatrix
		 *  M = R + D
		 * these are:
		 *  X_odd = -R(-mu)_inv D X_even
		 *  Y_odd = -R(mu)_inv D Y_even
		 */

		/** @fixme below usages of dslash should work, but only because we use bicgstab above
		    proper implementation needs to make sure this is always the case */

		///@NOTE the following calculations could also go in a new function for convenience
		//calculate X_odd
		//therefore, clmem_tmp_eo_1 is used as intermediate state. The result is saved in clmem_inout, since
		//  this is used as a default in the force-function.
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			fermion_code->dslash_eo_device(fermion_code->get_inout_eo(), fermion_code->get_tmp_eo_1(), &new_u, ODD, kappa);
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_1(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			fermion_code->dslash_eo_device(fermion_code->get_inout_eo(), fermion_code->get_tmp_eo_1(), &new_u, ODD, kappa);
			fermion_code->M_tm_inverse_sitediagonal_minus_device(fermion_code->get_tmp_eo_1(), fermion_code->get_tmp_eo_2(), mubar);
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_2(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		}
		//logger.debug() << "\t\tcalc eo fermion_force F(Y_even, X_odd)...";
		//Calc F(Y_even, X_odd) = F(clmem_phi_inv_eo, clmem_tmp_eo_1)
		fermion_force_eo_device(&clmem_phi_inv_eo,  fermion_code->get_tmp_eo_1(), EVEN, kappa);

		//calculate Y_odd
		//therefore, clmem_tmp_eo_1 is used as intermediate state. The result is saved in clmem_phi_inv, since
		//  this is used as a default in the force-function.
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			fermion_code->dslash_eo_device(&clmem_phi_inv_eo, fermion_code->get_tmp_eo_1(), &new_u, ODD, kappa);
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_1(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			fermion_code->dslash_eo_device(&clmem_phi_inv_eo, fermion_code->get_tmp_eo_1(), &new_u, ODD, kappa);
			fermion_code->M_tm_inverse_sitediagonal_device(fermion_code->get_tmp_eo_1(), fermion_code->get_tmp_eo_2(), mubar);
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_2(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		}

		//logger.debug() << "\t\tcalc eoprec fermion_force F(Y_odd, X_even)...";
		//Calc F(Y_odd, X_even) = F(clmem_tmp_eo_1, clmem_inout_eo)
		fermion_force_eo_device(fermion_code->get_tmp_eo_1(), fermion_code->get_inout_eo(), ODD, kappa);
	} else {
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
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
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			converged = fermion_code->cg(hardware::code::QplusQminus(fermion_code), fermion_code->get_inout(), get_clmem_phi(), &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			/**
			 * Y is now just
			 *  Y = (Qminus) X = (Qminus) (Qplusminus)^-1 psi =
			 *    = (Qplus)^-1 psi
			 */
			fermion_code->Qminus(fermion_code->get_inout(), &clmem_phi_inv, &new_u, kappa, mubar);

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
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			converged = fermion_code->bicgstab(hardware::code::Qplus(fermion_code), fermion_code->get_inout(), get_clmem_phi(), &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv, fermion_code->get_inout());

			/**
			 * Now, one has to calculate
			 * X = (Qminus)^-1 Y = (Qminus)^-1 (Qplus)^-1 psi = (QplusQminus)^-1 psi ??
			 * out of
			 * Qminus clmem_inout = clmem_phi_inv
			 * Therefore, when calculating the final energy of the spinorfield,
			 *  one can just take clmem_phi_inv (see also above)!!
			 */

			//copy former solution to clmem_source
			hardware::buffers::copyData(fermion_code->get_source(), fermion_code->get_inout());
			logger.debug() << "\t\t\tstart solver";

			//this sets clmem_inout cold as trial-solution
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			converged = fermion_code->bicgstab(hardware::code::Qminus(fermion_code), fermion_code->get_inout(), fermion_code->get_source(), &new_u, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

		}
		if(logger.beDebug()) {
			fermion_code->print_info_inv_field(&clmem_phi_inv, false, "\tY ");
			fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tX ");
		}
		logger.debug() << "\t\tcalc fermion_force...";
		fermion_force_device(&clmem_phi_inv, fermion_code->get_inout(), kappa);
	}
}

void hardware::code::Hmc::calc_fermion_force_detratio(usetimer * solvertimer, const hardware::buffers::SU3 * gaugefield)
{
	logger.debug() << "\t\tcalc fermion_force_detratio...";
	auto fermion_code = get_device()->get_fermion_code();
	auto spinor_code = get_device()->get_spinor_code();

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
	int converged = -1;
	hmc_float kappa = get_parameters().get_kappa();
	hmc_float mubar = get_mubar(get_parameters());
	hmc_float kappa2 = get_parameters().get_kappa_mp();
	hmc_float mubar2 = meta::get_mubar_mp(get_parameters());
	if(get_parameters().get_use_eo() == true) {
		//CP: Init tmp spinorfield
		hardware::buffers::Spinor sf_eo_tmp(clmem_phi_eo.get_elements(), get_device());
		//the source is now Q_2^+ phi = sf_eo_tmp
		fermion_code->Qplus_eo (get_clmem_phi_eo(), &sf_eo_tmp , gaugefield, kappa2, mubar2);
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
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
			spinor_code->set_eoprec_spinorfield_cold_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field before inversion ");
			converged = fermion_code->cg_eo(hardware::code::QplusQminus_eo(fermion_code), fermion_code->get_inout_eo(), &sf_eo_tmp, gaugefield, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\t\t\tinv. field after inversion ");

			/**
			 * Y_even is now just
			 *  Y_even = (Qminus_eo) X_even = (Qminus_eo) (Qplusminus_eo)^-1 sf_eo_tmp =
			 *    = (Qplus_eo)^-1 Q_2^+ psi
			 */
			fermion_code->Qminus_eo(fermion_code->get_inout_eo(), &clmem_phi_inv_eo, gaugefield, kappa, mubar);
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
			spinor_code->set_zero_spinorfield_eoprec_device(fermion_code->get_inout_eo());
			fermion_code->gamma5_eo_device(fermion_code->get_inout_eo());

			//if(logger.beDebug()) fermion_code->print_info_inv_field(get_inout_eo(), true, "\t\t\tinv. field before inversion ");
			//if(logger.beDebug()) print_info_inv_field(get_clmem_phi_eo(), true, "\t\t\tsource before inversion ");
			converged = fermion_code->bicgstab_eo(hardware::code::Qplus_eo(fermion_code), fermion_code->get_inout_eo(), &sf_eo_tmp, gaugefield, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			//if(logger.beDebug()) print_info_inv_field(get_inout_eo(), true, "\t\t\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv_eo, fermion_code->get_inout_eo());

			/**
			 * Now, one has to calculate
			 * X = (Qminus_eo)^-1 Y = (Qminus_eo)^-1 (Qplus_eo)^-1 sf_eo_tmp = (QplusQminus_eo)^-1 sf_eo_tmp = (QplusQminus_eo)^-1 Q_2^+ psi
			 * out of
			 * Qminus_eo clmem_inout_eo = clmem_phi_inv_eo
			 */

			logger.debug() << "\t\tcalc X_even...";
			//copy former solution to clmem_source
			hardware::buffers::copyData(fermion_code->get_source_even(), fermion_code->get_inout_eo());

			//this sets clmem_inout cold as trial-solution
			spinor_code->set_eoprec_spinorfield_cold_device(fermion_code->get_inout_eo());

			//if(logger.beDebug()) print_info_inv_field(get_inout_eo(), true, "\t\t\tinv. field before inversion ");
			//if(logger.beDebug()) print_info_inv_field(get_source_even(), true, "\t\t\tsource before inversion ");
			converged = fermion_code->bicgstab_eo(hardware::code::Qminus_eo(fermion_code), fermion_code->get_inout_eo(), fermion_code->get_source_even(), gaugefield, get_parameters().get_force_prec(), kappa, mubar);
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			//if(logger.beDebug()) print_info_inv_field(get_inout_eo(), true, "\t\t\tinv. field after inversion ");
		}
		/**
		 * At this point, one has to calculate X_odd and Y_odd.
		 * If one has a fermionmatrix
		 *  M = R + D
		 * these are:
		 *  X_odd = -R(-mu)_inv D X_even
		 *  Y_odd = -R(mu)_inv D Y_even
		 */

		/** @fixme below usages of dslash should work, but only because we use bicgstab above
		proper implementation needs to make sure this is always the case 
		CP: tmp1 and tmp2 are initialized independently of bicgstab as far as I know! */

		///@NOTE the following calculations could also go in a new function for convenience
		//calculate X_odd
		//therefore, sf_eo_tmp is used as intermediate state. 
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			fermion_code->dslash_eo_device(fermion_code->get_inout_eo(), fermion_code->get_tmp_eo_1(), gaugefield, ODD);
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_1(), fermion_code->get_clmem_minusone(), &sf_eo_tmp);
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			fermion_code->dslash_eo_device(fermion_code->get_inout_eo(), fermion_code->get_tmp_eo_1(), gaugefield, ODD);
			fermion_code->M_tm_inverse_sitediagonal_minus_device(fermion_code->get_tmp_eo_1(), fermion_code->get_tmp_eo_2());
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_2(), fermion_code->get_clmem_minusone(), &sf_eo_tmp);
		}

		//logger.debug() << "\t\tcalc eo fermion_force F(Y_even, X_odd)...";
		//Calc F(Y_even, X_odd) = F(clmem_phi_inv_eo, sf_eo_tmp)
		fermion_force_eo_device(&clmem_phi_inv_eo,  &sf_eo_tmp, EVEN, kappa);

		//calculate Y_odd
		//therefore, clmem_tmp_eo_1 is used as intermediate state. 
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
			fermion_code->dslash_eo_device(&clmem_phi_inv_eo, fermion_code->get_tmp_eo_1(), gaugefield, ODD);
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_1(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
			fermion_code->dslash_eo_device(&clmem_phi_inv_eo, fermion_code->get_tmp_eo_1(), gaugefield, ODD);
			fermion_code->M_tm_inverse_sitediagonal_device(fermion_code->get_tmp_eo_1(), fermion_code->get_tmp_eo_2());
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_2(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		}

		//logger.debug() << "\t\tcalc eoprec fermion_force F(Y_odd, X_even)...";
		//Calc F(Y_odd, X_even) = F(clmem_tmp_eo_1, clmem_inout_eo)
		fermion_force_eo_device(fermion_code->get_tmp_eo_1(), fermion_code->get_inout_eo(), ODD, kappa);

		/**
		 *Now, one has the additional term - phi^+ deriv(Q_2) X
		 *This works exactly the same as above, except replacing Y -> -phi and different mass parameters
		 */

		//Y is not needed anymore, therefore use clmem_phi_inv_eo to store -phi
		spinor_code->sax_eoprec_device(get_clmem_phi_eo(), fermion_code->get_clmem_minusone(), &clmem_phi_inv_eo);

		//logger.debug() << "\t\tcalc eo fermion_force F(Y_even, X_odd)...";
		//Calc F(Y_even, X_odd) = F(clmem_phi_inv_eo, sf_eo_tmp)
		fermion_force_eo_device(&clmem_phi_inv_eo,  &sf_eo_tmp, EVEN, kappa2);

		//calculate phi_odd
		//this works in the same way as with Y above, since -phi_even is saved in the same buffer as Y_even
		//therefore, clmem_tmp_eo_1 is used as intermediate state. 
		if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
		  fermion_code->dslash_eo_device(&clmem_phi_inv_eo, fermion_code->get_tmp_eo_1(), gaugefield, ODD, kappa2);
			spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_1(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		} else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
		  fermion_code->dslash_eo_device(&clmem_phi_inv_eo, fermion_code->get_tmp_eo_1(), gaugefield, ODD, kappa2);
		  fermion_code->M_tm_inverse_sitediagonal_device(fermion_code->get_tmp_eo_1(), fermion_code->get_tmp_eo_2(), mubar2);
		  spinor_code->sax_eoprec_device(fermion_code->get_tmp_eo_2(), fermion_code->get_clmem_minusone(), fermion_code->get_tmp_eo_1());
		}

		//logger.debug() << "\t\tcalc eoprec fermion_force F(Y_odd, X_even)...";
		//Calc F(Y_odd, X_even) = F(clmem_tmp_eo_1, clmem_inout_eo)
		fermion_force_eo_device(fermion_code->get_tmp_eo_1(), fermion_code->get_inout_eo(), ODD, kappa2);
	} else {
		//CP: Init tmp spinorfield
		hardware::buffers::Plain<spinor> sf_tmp(clmem_phi.get_elements(), get_device());
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		//the source is now Q_2^+ phi = sf_tmp
		fermion_code->Qplus (get_clmem_phi(), &sf_tmp , gaugefield, kappa2, mubar2);
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
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
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			converged = fermion_code->cg(hardware::code::QplusQminus(fermion_code), fermion_code->get_inout(), &sf_tmp, gaugefield, get_parameters().get_force_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			/**
			 * Y is now just
			 *  Y = (Qminus) X = (Qminus) (Qplusminus)^-1 sf_tmp =
			 *    = (Qplus)^-1 sf_tmp
			 */
			fermion_code->Qminus(fermion_code->get_inout(), &clmem_phi_inv, gaugefield);

		} else  {
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
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			//here, the "normal" solver can be used since the inversion is of the same structure as in the inverter
			converged = fermion_code->bicgstab(hardware::code::Qplus(fermion_code), fermion_code->get_inout(), &sf_tmp, gaugefield, get_parameters().get_force_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv, fermion_code->get_inout());

			/**
			 * Now, one has to calculate
			 * X = (Qminus)^-1 Y = (Qminus)^-1 (Qplus)^-1 psi = (QplusQminus)^-1 psi ??
			 * out of
			 * Qminus clmem_inout = clmem_phi_inv
			 * Therefore, when calculating the final energy of the spinorfield,
			 *  one can just take clmem_phi_inv (see also above)!!
			 */

			//copy former solution to clmem_source
			hardware::buffers::copyData(fermion_code->get_source(), fermion_code->get_inout());
			logger.debug() << "\t\t\tstart solver";

			//this sets clmem_inout cold as trial-solution
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			converged = fermion_code->bicgstab(hardware::code::Qminus(fermion_code), fermion_code->get_inout(), fermion_code->get_source(), gaugefield, get_parameters().get_force_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

		}
		if(logger.beDebug()) {
			fermion_code->print_info_inv_field(&clmem_phi_inv, false, "\tY ");
			fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tX ");
		}
		logger.debug() << "\t\tcalc fermion_force...";
		fermion_force_device(&clmem_phi_inv, fermion_code->get_inout(), kappa);

		/**
		 *Now, one has the additional term - phi^+ deriv(Q_2) X
		 *This works exactly the same as above, except replacing Y -> -phi and different mass parameters
		 */

		//Y is not needed anymore, therefore use clmem_phi_inv_eo to store -phi
		spinor_code->sax_device(get_clmem_phi(), fermion_code->get_clmem_minusone(), &clmem_phi_inv);

		fermion_force_device(&clmem_phi_inv, fermion_code->get_inout(), kappa2);
	}
}

void hardware::code::Hmc::calc_gauge_force()
{
	logger.debug() << "\t\tcalc gauge_force...";
	gauge_force_device();
	if(meta::get_use_rectangles(get_parameters()) == true) {
		logger.debug() << "\t\tcalc rect gauge_force...";
		gauge_force_tlsym_device();
	}
}

hmc_float hardware::code::Hmc::calc_s_fermion()
{
	auto fermion_code = get_device()->get_fermion_code();
	auto spinor_code = get_device()->get_spinor_code();

	logger.debug() << "HMC:\tcalc final fermion energy...";
	//this function essentially performs the same steps as in the force-calculation, but with higher precision.
	//  therefore, comments are deleted here...
	//  Furthermore, in the bicgstab-case, the second inversions are not needed
	int converged = -1;
	if(get_parameters().get_use_eo() == true) {
		//the source is already set, it is Dpsi, where psi is the initial gaussian spinorfield
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			logger.debug() << "\t\t\tstart solver";

			spinor_code->set_eoprec_spinorfield_cold_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field before inversion ");
			converged = fermion_code->cg_eo(hardware::code::QplusQminus_eo(fermion_code), fermion_code->get_inout_eo(), get_clmem_phi_eo(), &new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field after inversion ");

			fermion_code->Qminus_eo(fermion_code->get_inout_eo(), &clmem_phi_inv_eo, &new_u);
		} else {
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			spinor_code->set_zero_spinorfield_eoprec_device(fermion_code->get_inout_eo());
			fermion_code->gamma5_eo_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field before inversion ");
			if(logger.beDebug()) fermion_code->print_info_inv_field(get_clmem_phi_eo(), true, "\tsource before inversion ");
			converged = fermion_code->bicgstab_eo(hardware::code::Qplus_eo(fermion_code), fermion_code->get_inout_eo(), get_clmem_phi_eo(), &new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv_eo, fermion_code->get_inout_eo());
		}
	} else {
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			converged = fermion_code->cg(hardware::code::QplusQminus(fermion_code), fermion_code->get_inout(), get_clmem_phi(), &new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			fermion_code->Qminus(fermion_code->get_inout(), &clmem_phi_inv, &new_u);

		} else  {

			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			  * or a point-source spinorfield as trial-solution
			  */
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			converged = fermion_code->bicgstab(hardware::code::Qplus(fermion_code), fermion_code->get_inout(), get_clmem_phi(), &new_u, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv, fermion_code->get_inout());
		}
	}
	///@todo: this can be moved in the ifs above!!
	if(get_parameters().get_use_eo() == true) {
		get_device()->get_spinor_code()->set_float_to_global_squarenorm_eoprec_device(&clmem_phi_inv_eo, &clmem_s_fermion);
	} else {
		get_device()->get_spinor_code()->set_float_to_global_squarenorm_device(&clmem_phi_inv, &clmem_s_fermion);
	}
	hmc_float tmp;
	clmem_s_fermion.dump(&tmp);
	return tmp;
}

hmc_float hardware::code::Hmc::calc_s_fermion_mp(const hardware::buffers::SU3 * gaugefield)
{
	auto fermion_code = get_device()->get_fermion_code();
	auto spinor_code = get_device()->get_spinor_code();

	logger.debug() << "HMC:\tcalc final fermion energy...";
	//this function essentially performs the same steps as in the non mass-prec case, however, one has to apply one more matrix multiplication
	//  therefore, comments are deleted here...
	//  Furthermore, in the bicgstab-case, the second inversions are not needed
	int converged = -1;
	if(get_parameters().get_use_eo() == true) {
		//CP: Init tmp spinorfield
		hardware::buffers::Spinor sf_eo_tmp(clmem_phi_mp_eo.get_elements(), get_device());

		//sf_eo_tmp = Qplus_eo(heavy_mass) phi_mp_eo
		fermion_code->Qplus_eo (get_clmem_phi_mp_eo(), &sf_eo_tmp , gaugefield, get_parameters().get_kappa_mp(), meta::get_mubar_mp(get_parameters()));
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			logger.debug() << "\t\t\tstart solver";
			spinor_code->set_eoprec_spinorfield_cold_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field before inversion ");
			converged = fermion_code->cg_eo(hardware::code::QplusQminus_eo(fermion_code), fermion_code->get_inout_eo(), &sf_eo_tmp, gaugefield, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field after inversion ");

			fermion_code->Qminus_eo(fermion_code->get_inout_eo(), &clmem_phi_inv_eo, gaugefield);
		} else {
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			spinor_code->set_zero_spinorfield_eoprec_device(fermion_code->get_inout_eo());
			fermion_code->gamma5_eo_device(fermion_code->get_inout_eo());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field before inversion ");
			if(logger.beDebug()) fermion_code->print_info_inv_field(&sf_eo_tmp, true, "\tsource before inversion ");
			converged = fermion_code->bicgstab_eo(hardware::code::Qplus_eo(fermion_code), fermion_code->get_inout_eo(), &sf_eo_tmp, gaugefield, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout_eo(), true, "\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv_eo, fermion_code->get_inout_eo());
		}
	} else {
		//CP: Init tmp spinorfield
		hardware::buffers::Plain<spinor> sf_tmp(clmem_phi_mp.get_elements(), get_device());

		//sf_eo_tmp = Qplus(light_mass) phi_mp
		fermion_code->Qplus (get_clmem_phi_mp(), &sf_tmp , gaugefield, get_parameters().get_kappa_mp(), meta::get_mubar_mp(get_parameters()));
		if(get_parameters().get_solver() == meta::Inputparameters::cg) {
			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			converged = fermion_code->cg(hardware::code::QplusQminus(fermion_code), fermion_code->get_inout(), &sf_tmp, gaugefield, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			fermion_code->Qminus(fermion_code->get_inout(), &clmem_phi_inv, gaugefield);

		} else  {

			logger.debug() << "\t\t\tstart solver";

			/** @todo at the moment, we can only put in a cold spinorfield
			 * or a point-source spinorfield as trial-solution
			 */
			spinor_code->set_spinorfield_cold_device(fermion_code->get_inout());

			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field before inversion ");
			converged = fermion_code->bicgstab(hardware::code::Qplus(fermion_code), fermion_code->get_inout(), &sf_tmp, gaugefield, get_parameters().get_solver_prec());
			if (converged < 0) {
				if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
				else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
			} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
			if(logger.beDebug()) fermion_code->print_info_inv_field(fermion_code->get_inout(), false, "\tinv. field after inversion ");

			//store this result in clmem_phi_inv
			hardware::buffers::copyData(&clmem_phi_inv, fermion_code->get_inout());
		}
	}
	///@todo: this can be moved in the ifs above!!
	if(get_parameters().get_use_eo() == true) {
		get_device()->get_spinor_code()->set_float_to_global_squarenorm_eoprec_device(&clmem_phi_inv_eo, &clmem_s_fermion);
	} else {
		get_device()->get_spinor_code()->set_float_to_global_squarenorm_device(&clmem_phi_inv, &clmem_s_fermion);
	}
	hmc_float tmp;
	clmem_s_fermion.dump(&tmp);
	return tmp;
}

hmc_observables hardware::code::Hmc::metropolis(hmc_float rnd, hmc_float beta, const hardware::buffers::SU3 * gaugefield)
{
	auto gf_code = get_device()->get_gaugefield_code();
	auto fermion_code = get_device()->get_fermion_code();
	auto gm_code = get_device()->get_gaugemomentum_code();

	//Calc Hamiltonian
	logger.debug() << "HMC:\tCalculate Hamiltonian";
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
	gf_code->gaugeobservables(gaugefield, &plaq,  &tplaq, &splaq, &poly);
	gf_code->gaugeobservables(&new_u, &plaq_new,  &tplaq_new, &splaq_new, &poly_new);
	//plaq has to be divided by the norm-factor to get s_gauge
	hmc_float factor = 1. / (meta::get_plaq_norm(get_parameters()));
	if(meta::get_use_rectangles(get_parameters()) == true) {
		gf_code->gaugeobservables_rectangles(gaugefield, &rect);
		gf_code->gaugeobservables_rectangles(&new_u, &rect_new);
		hmc_float c0 = meta::get_c0(get_parameters());
		hmc_float c1 = meta::get_c1(get_parameters());
		deltaH = - beta * ( c0 * (plaq - plaq_new) / factor + c1 * ( rect - rect_new )  );
		s_old = - beta * ( c0 * (plaq) / factor + c1 * ( rect )  );
		s_new = - beta * ( c0 * (plaq_new) / factor + c1 * ( rect_new )  );

	} else {
		/** NOTE: the minus here is introduced to fit tmlqcd!!! */
		deltaH = -(plaq - plaq_new) * beta / factor;
		s_old = -(plaq ) * beta / factor;
		s_new = -(plaq_new) * beta / factor;
	}

	logger.debug() << "HMC:\tS_gauge(old field) = " << setprecision(10) << s_old;
	logger.debug() << "HMC:\tS_gauge(new field) = " << setprecision(10) << s_new;
	logger.info()  << "HMC:\tdeltaS_gauge = " << setprecision(10) << deltaH;

	//Gaugemomentum-Part
	hmc_float p2, new_p2;
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&clmem_p, &clmem_p2);
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&clmem_new_p, &clmem_new_p2);
	clmem_p2.dump(&p2);
	clmem_new_p2.dump(&new_p2);
	//the energy is half the squarenorm
	deltaH += 0.5 * (p2 - new_p2);

	logger.debug() << "HMC:\tS_gaugemom(old field) = " << setprecision(10) << 0.5 * p2;
	logger.debug() << "HMC:\tS_gaugemom(new field) = " << setprecision(10) << 0.5 * new_p2;
	logger.info()  << "HMC:\tdeltaS_gaugemom = " << setprecision(10) << 0.5 * (p2 - new_p2);

	//Fermion-Part:
	if(! get_parameters().get_use_gauge_only() ) {
		hmc_float spinor_energy_init, s_fermion_final;
		//initial energy has been computed in the beginning...
		clmem_s_fermion_init.dump(&spinor_energy_init);
		// sum_links phi*_i (M^+M)_ij^-1 phi_j
		s_fermion_final = calc_s_fermion();
		deltaH += spinor_energy_init - s_fermion_final;

		logger.debug() << "HMC:\tS_ferm(old field) = " << setprecision(10) <<  spinor_energy_init;
		logger.debug() << "HMC:\tS_ferm(new field) = " << setprecision(10) << s_fermion_final;
		logger.info() <<  "HMC:\tdeltaS_ferm = " << spinor_energy_init - s_fermion_final;
		if( get_parameters().get_use_mp() ) {
			hmc_float spinor_energy_mp_init, s_fermion_mp_final;
			//initial energy has been computed in the beginning...
			clmem_s_fermion_mp_init.dump(&spinor_energy_mp_init);
			// sum_links phi*_i (M^+M)_ij^-1 phi_j
			s_fermion_mp_final = calc_s_fermion_mp(&new_u);
			deltaH += spinor_energy_mp_init - s_fermion_mp_final;

			logger.debug() << "HMC:\tS_ferm_mp(old field) = " << setprecision(10) <<  spinor_energy_mp_init;
			logger.debug() << "HMC:\tS_ferm_mp(new field) = " << setprecision(10) << s_fermion_mp_final;
			logger.info() <<  "HMC:\tdeltaS_ferm_mp = " << spinor_energy_mp_init - s_fermion_mp_final;
		}
	}
	//Metropolis-Part
	hmc_float compare_prob;
	if(deltaH < 0) {
		compare_prob = exp(deltaH);
	} else {
		compare_prob = 1.0;
	}
	logger.info() << "HMC:\tdeltaH = " << deltaH << "\tAcc-Prop = " << compare_prob;
	hmc_observables tmp;
	if(rnd <= compare_prob) {
		tmp.accept = 1;
		tmp.plaq = plaq_new;
		tmp.tplaq = tplaq_new;
		tmp.splaq = splaq_new;
		tmp.poly = poly_new;
		tmp.deltaH = deltaH;
		tmp.prob = compare_prob;
		if(meta::get_use_rectangles(get_parameters()) ) tmp.rectangles = rect_new / get_rect_norm(get_parameters());
	} else {
		tmp.accept = 0;
		tmp.plaq = plaq;
		tmp.tplaq = tplaq;
		tmp.splaq = splaq;
		tmp.poly = poly;
		tmp.deltaH = deltaH;
		tmp.prob = compare_prob;
		if(meta::get_use_rectangles(get_parameters()) ) tmp.rectangles = rect / get_rect_norm(get_parameters());
	}

	return tmp;
}

void hardware::code::Hmc::calc_spinorfield_init_energy(hardware::buffers::Plain<hmc_float> * dest)
{
	auto spinor_code = get_device()->get_spinor_code();

	//Suppose the initial spinorfield is saved in phi_inv
	//  it is created in generate_gaussian_spinorfield_device
	if(get_parameters().get_use_eo() == true) {
		spinor_code->set_float_to_global_squarenorm_eoprec_device(&clmem_phi_inv_eo, dest);
	} else {
		spinor_code->set_float_to_global_squarenorm_device(&clmem_phi_inv, dest);
	}
}

void hardware::code::Hmc::md_update_gaugemomentum_device(hmc_float eps)
{
	auto mol_dyn_code = get_device()->get_molecular_dynamics_code();
	mol_dyn_code->md_update_gaugemomentum_device(&clmem_force, &clmem_new_p, eps);
}

void hardware::code::Hmc::md_update_gaugefield_device(hmc_float eps)
{
	auto mol_dyn_code = get_device()->get_molecular_dynamics_code();
	mol_dyn_code->md_update_gaugefield_device(&clmem_new_p, &new_u, eps);
}

void hardware::code::Hmc::set_zero_clmem_force_device()
{
	auto gm_code = get_device()->get_gaugemomentum_code();
	gm_code->set_zero_gaugemomentum(&clmem_force);
}

void hardware::code::Hmc::gauge_force_device()
{
	auto mol_dyn_code = get_device()->get_molecular_dynamics_code();
	mol_dyn_code->gauge_force_device(&new_u, &clmem_force);
}
void hardware::code::Hmc::gauge_force_tlsym_device()
{
	auto mol_dyn_code = get_device()->get_molecular_dynamics_code();
	mol_dyn_code->gauge_force_tlsym_device(&new_u, &clmem_force);
}
void hardware::code::Hmc::fermion_force_device(const hardware::buffers::Plain<spinor> * Y, const hardware::buffers::Plain<spinor> * X, hmc_float kappa)
{
	using namespace hardware::buffers;
	auto mol_dyn_code = get_device()->get_molecular_dynamics_code();
	mol_dyn_code->fermion_force_device(Y, X, &new_u, &clmem_force, kappa);
}

//the argument kappa is set to ARG_DEF as default
void hardware::code::Hmc::fermion_force_eo_device(const hardware::buffers::Spinor * Y, const hardware::buffers::Spinor * X, int evenodd, hmc_float kappa)
{
	using namespace hardware::buffers;
	auto mol_dyn_code = get_device()->get_molecular_dynamics_code();

	mol_dyn_code->fermion_force_eo_device(Y, X, &new_u, &clmem_force, evenodd, kappa);
}
