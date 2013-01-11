/** @file
 * Implementation of explicit fermionamtrix operations
 */

#include "fermionmatrix.hpp"

void physics::fermionmatrix::M(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa, hmc_float mubar, const hardware::System& system)
{
	switch(system.get_inputparameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			//in the pure Wilson case there is just one fermionmatrix
			M_wilson(out, gf, in, kappa);
			break;
		case meta::Inputparameters::twistedmass:
			M_tm_plus(out, gf, in, kappa, mubar);
			break;
	}
}

void physics::fermionmatrix::Qplus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa, hmc_float mubar, const hardware::System& system)
{
	switch(system.get_inputparameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			//in the pure Wilson case there is just one fermionmatrix
			M_wilson(out, gf, in, kappa);
			break;
		case meta::Inputparameters::twistedmass:
			M_tm_plus(out, gf, in, kappa, mubar);
			break;
	}
	out->gamma5();
}

void physics::fermionmatrix::Qminus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa, hmc_float mubar, const hardware::System& system)
{
	switch(system.get_inputparameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			//in the pure Wilson case there is just one fermionmatrix
			M_wilson(out, gf, in, kappa);
			break;
		case meta::Inputparameters::twistedmass:
			M_tm_minus(out, gf, in, kappa, mubar);
			break;
	}
	out->gamma5();
}

void physics::fermionmatrix::QplusQminus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa, hmc_float mubar, const hardware::System& system)
{
	/** @todo one could save one field here if an additional copying would be included in the end...
	 * or the field should be created in here, local */
	/** @todo The local creation of the temporary field is known to cause performance problems... */
	physics::lattices::Spinorfield tmp(system);
	Qminus(&tmp, gf, in, kappa, mubar, system);
	Qplus(out, gf, tmp, kappa, mubar, system);
}

void physics::fermionmatrix::Qplus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, hmc_float kappa, hmc_float mubar, const hardware::System& system)
{
	if(system.get_inputparameters().get_use_merge_kernels_fermion() == false) {
		Aee(out, gf, in, kappa, mubar, system);
		out->gamma5();
	} else {
		throw Print_Error_Message("kernel merging not implemented for Qplus", __FILE__, __LINE__);
		//Aee_AND_gamma5_eo(in, out, gf, kappa, mubar);
	}
}

void physics::fermionmatrix::Qminus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, hmc_float kappa, hmc_float mubar, const hardware::System& system)
{
	if(system.get_inputparameters().get_use_merge_kernels_fermion() == false) {
		Aee_minus(out, gf, in, kappa, mubar, system);
		out->gamma5();
	} else {
		throw Print_Error_Message("kernel merging not implemented for Qminus", __FILE__, __LINE__);
		//Aee_minus_AND_gamma5_eo(in, out, gf, kappa, mubar);
	}
}

void physics::fermionmatrix::QplusQminus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, hmc_float kappa, hmc_float mubar, const hardware::System& system)
{
	//CP: this should be an eoprec-sized field. However, this induces problems in the CG algorithm!!!
	//MB: This is because of padding, the eoprec buffer size shoulw always be queried from Opencl_Module_Spinor
	/// @todo use a buffer from a pool
	/** @todo The local creation of the temporary field is known to cause performance problems... */
	physics::lattices::Spinorfield_eo tmp(system); // FIXME somehow pool this...
	Qminus(&tmp, gf, in, kappa, mubar, system);
	Qplus(out, gf, tmp, kappa, mubar, system);
}

void physics::fermionmatrix::Aee(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, hmc_float kappa, hmc_float mubar, const hardware::System& system)
{
	/**
	 * This is the even-odd preconditioned fermion matrix with the
	 * non-trivial inversion on the even sites (see DeGran/DeTar p. 174).
	 * If one has fermionmatrix
	 *  M = R + D,
	 * then Aee is:
	 * Aee = R_e - D_eo R_o_inv D_oe
	 */

	/** @todo The local creation of the temporary field is known to cause performance problems... */
	physics::lattices::Spinorfield_eo tmp(system); // FIXME somehow pool this...

	switch(system.get_inputparameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			//in this case, the diagonal matrix is just 1 and falls away.
			dslash(&tmp, gf, in, ODD, kappa);
			dslash(out, gf, tmp, EVEN, kappa);
			saxpy(out, {1., 0.}, *out, in);
			break;
		case meta::Inputparameters::twistedmass: {
			physics::lattices::Spinorfield_eo tmp2(system); // FIXME somehow pool this...
			dslash(&tmp, gf, in, ODD, kappa);
			M_tm_inverse_sitediagonal(&tmp2, tmp, mubar);
			dslash(out, gf, tmp2, EVEN, kappa);
			M_tm_sitediagonal(&tmp, in, mubar);
			saxpy(out, {1., 0.}, *out, tmp);
		}
		break;
	}
}

/**
 *  This is the equivalent of the above function, but for the lower
 *  flavour, which essentially means mu -> -mu in the tm-case and
 *  no changes in the meta::Inputparameters::wilson case.
 */
void physics::fermionmatrix::Aee_minus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, hmc_float kappa, hmc_float mubar, const hardware::System& system)
{
	/**
	 * This is the even-odd preconditioned fermion matrix with the
	 * non-trivial inversion on the even sites (see DeGran/DeTar p. 174).
	 * If one has fermionmatrix
	 *  M = R + D,
	 * then Aee is:
	 * Aee = R_e - D_eo R_o_inv D_oe
	 * and Aee_minus is:
	 * Aee_minus = R_e(-mu) - D_eo R_o(-mu)_inv D_oe
	 */
	/** @todo The local creation of the temporary field is known to cause performance problems... */
	physics::lattices::Spinorfield_eo tmp(system); // FIXME somehow pool this...

	switch(system.get_inputparameters().get_fermact()) {
		case meta::Inputparameters::wilson:
			//in this case, the diagonal matrix is just 1 and falls away.
			dslash(&tmp, gf, in, ODD, kappa);
			dslash(out, gf, tmp, EVEN, kappa);
			saxpy(out, {1., 0.}, *out, in);
			break;
		case meta::Inputparameters::twistedmass: {
			physics::lattices::Spinorfield_eo tmp2(system); // FIXME somehow pool this...
			dslash(&tmp, gf, in, ODD, kappa);
			M_tm_inverse_sitediagonal_minus(&tmp2, tmp, mubar);
			dslash(out, gf, tmp2, EVEN, kappa);
			M_tm_sitediagonal_minus(&tmp, in, mubar);
			saxpy(out, {1., 0.}, *out, tmp);
		}
		break;
	}
}
