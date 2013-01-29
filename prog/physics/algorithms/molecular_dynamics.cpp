/** @file
 * Implementation of the molecular dynamics algorithms
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * (c) 2012-2013 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
 */

#include "molecular_dynamics.hpp"

#include <stdexcept>
#include "../fermionmatrix/fermionmatrix.hpp"
#include "../../meta/util.hpp"
#include "solver.hpp"
#include "fermion_force.hpp"

static void md_update_gaugemomenta(const physics::lattices::Gaugemomenta * dest, const physics::lattices::Gaugemomenta& src, hmc_float eps);
static void md_update_gaugefield(const physics::lattices::Gaugefield * gf, const physics::lattices::Gaugemomenta& , hmc_float eps);
static void md_update_spinorfield(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& orig, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
static void md_update_spinorfield(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& orig, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
static void md_update_spinorfield_mp(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& orig, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
static void md_update_spinorfield_mp(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& orig, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);


static void md_update_gaugemomenta(const physics::lattices::Gaugemomenta * const dest, const physics::lattices::Gaugemomenta& src, const hmc_float eps)
{
	auto dest_bufs = dest->get_buffers();
	auto src_bufs = src.get_buffers();
	size_t num_bufs = dest_bufs.size();
	if(num_bufs != src_bufs.size()) {
		throw std::invalid_argument("The source and destination gaugemomenta use different devices.");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto dest_buf = dest_bufs[i];
		auto src_buf = src_bufs[i];
		auto code = dest_buf->get_device()->get_molecular_dynamics_code();
		code->md_update_gaugemomentum_device(src_buf, dest_buf, eps);
	}
}

static void md_update_gaugefield(const physics::lattices::Gaugefield * const gf, const physics::lattices::Gaugemomenta& gm, const hmc_float eps)
{
	auto gf_bufs = gf->get_buffers();
	auto gm_bufs = gm.get_buffers();
	size_t num_bufs = gf_bufs.size();
	if(num_bufs != gm_bufs.size()) {
		throw std::invalid_argument("The gaugefield and the gaugemomenta use different devices.");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto gf_buf = gf_bufs[i];
		auto gm_buf = gm_bufs[i];
		auto code = gf_buf->get_device()->get_molecular_dynamics_code();
		code->md_update_gaugefield_device(gm_buf, gf_buf, eps);
	}
}

void physics::algorithms::gauge_force(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield& gf)
{
	auto gm_bufs = gm->get_buffers();
	auto gf_bufs = gf.get_buffers();
	size_t num_bufs = gm_bufs.size();
	if(num_bufs != 1 || num_bufs != gf_bufs.size()) {
		throw Print_Error_Message(std::string(__func__) + " is only implemented for a single device.", __FILE__, __LINE__);
	}

	auto gm_buf = gm_bufs[0];
	auto gf_buf = gf_bufs[0];
	auto code = gm_buf->get_device()->get_molecular_dynamics_code();
	code->gauge_force_device(gf_buf, gm_buf);
}

void physics::algorithms::gauge_force_tlsym(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield& gf)
{
	auto gm_bufs = gm->get_buffers();
	auto gf_bufs = gf.get_buffers();
	size_t num_bufs = gm_bufs.size();
	if(num_bufs != 1 || num_bufs != gf_bufs.size()) {
		throw Print_Error_Message(std::string(__func__) + " is only implemented for a single device.", __FILE__, __LINE__);
	}

	auto gm_buf = gm_bufs[0];
	auto gf_buf = gf_bufs[0];
	auto code = gm_buf->get_device()->get_molecular_dynamics_code();
	code->gauge_force_tlsym_device(gf_buf, gm_buf);
}

void physics::algorithms::fermion_force(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Spinorfield& Y, const physics::lattices::Spinorfield& X, const physics::lattices::Gaugefield& gf, const hmc_float kappa)
{
	auto gm_bufs = gm->get_buffers();
	auto Y_bufs = Y.get_buffers();
	auto X_bufs = X.get_buffers();
	auto gf_bufs = gf.get_buffers();
	size_t num_bufs = gm_bufs.size();
	if(num_bufs != 1 || num_bufs != Y_bufs.size() || num_bufs != X_bufs.size() || num_bufs != gf_bufs.size()) {
		throw Print_Error_Message(std::string(__func__) + " is only implemented for a single device.", __FILE__, __LINE__);
	}

	auto gm_buf = gm_bufs[0];
	auto Y_buf = Y_bufs[0];
	auto X_buf = X_bufs[0];
	auto gf_buf = gf_bufs[0];
	auto code = gm_buf->get_device()->get_molecular_dynamics_code();
	code->fermion_force_device(Y_buf, X_buf, gf_buf, gm_buf, kappa);
}

void physics::algorithms::fermion_force(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Spinorfield_eo& Y, const physics::lattices::Spinorfield_eo& X, const int evenodd, const physics::lattices::Gaugefield& gf, const hmc_float kappa)
{
	auto gm_bufs = gm->get_buffers();
	auto Y_bufs = Y.get_buffers();
	auto X_bufs = X.get_buffers();
	auto gf_bufs = gf.get_buffers();
	size_t num_bufs = gm_bufs.size();
	if(num_bufs != 1 || num_bufs != Y_bufs.size() || num_bufs != X_bufs.size() || num_bufs != gf_bufs.size()) {
		throw Print_Error_Message(std::string(__func__) + " is only implemented for a single device.", __FILE__, __LINE__);
	}

	auto gm_buf = gm_bufs[0];
	auto Y_buf = Y_bufs[0];
	auto X_buf = X_bufs[0];
	auto gf_buf = gf_bufs[0];
	auto code = gm_buf->get_device()->get_molecular_dynamics_code();
	code->fermion_force_eo_device(Y_buf, X_buf, gf_buf, gm_buf, evenodd, kappa);
}

static void md_update_spinorfield(const physics::lattices::Spinorfield * const out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& orig, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
	physics::fermionmatrix::Qplus qplus(kappa, mubar, system);
	qplus(out, gf, orig);
	trace_squarenorm("Spinorfield after update", *out);
}

static void md_update_spinorfield(const physics::lattices::Spinorfield_eo * const out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& orig, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
	physics::fermionmatrix::Qplus_eo qplus(kappa, mubar, system);
	qplus(out, gf, orig);
	trace_squarenorm("Spinorfield after update", *out);
}

template<class FERMIONMATRIX, class SPINORFIELD> static void md_update_spinorfield_mp(const SPINORFIELD * const out, const physics::lattices::Gaugefield& gf, const SPINORFIELD& orig, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
	using physics::algorithms::solvers::bicgstab;

	SPINORFIELD tmp(system);
	FERMIONMATRIX qplus(kappa, mubar, system);

	trace_squarenorm("Spinorfield before update", orig);

	qplus(&tmp, gf, orig);

	//Now one needs ( Qplus )^-1 (heavy_mass) using tmp as source to get phi_mp
	//use always bicgstab here
	logger.debug() << "\t\t\tstart solver";

	/** @todo at the moment, we can only put in a cold spinorfield
	 * or a point-source spinorfield as trial-solution
	 */
	out->zero();
	out->gamma5();

	trace_squarenorm("\tinv. field before inversion.", *out);
	trace_squarenorm("\tsource before inversion.", tmp);

	auto params = system.get_inputparameters();
	FERMIONMATRIX qplus_mp(params.get_kappa_mp(), meta::get_mubar_mp(params), system);
	const int iterations = bicgstab(out, qplus_mp, gf, tmp, system, params.get_solver_prec());
	logger.debug() << "\t\t\tsolver solved in " << iterations << " iterations!";

	trace_squarenorm("\tinv. field after inversion ", *out);
}

static void md_update_spinorfield_mp(const physics::lattices::Spinorfield * const out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& orig, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
	using physics::fermionmatrix::Qplus;

	::md_update_spinorfield_mp<Qplus>(out, gf, orig, system, kappa, mubar);
}

static void md_update_spinorfield_mp(const physics::lattices::Spinorfield_eo * const out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& orig, const hardware::System& system, const hmc_float kappa, const hmc_float mubar)
{
	using physics::fermionmatrix::Qplus_eo;

	::md_update_spinorfield_mp<Qplus_eo>(out, gf, orig, system, kappa, mubar);
}

template<class SPINORFIELD> static void md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const SPINORFIELD& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	using namespace physics::algorithms;

	physics::lattices::Gaugemomenta delta_p(system);
	delta_p.zero();
	calc_total_force(&delta_p, gf, phi, system, kappa, mubar);

	logger.debug() << "\tHMC [UP]:\tupdate GM [" << eps << "]";
	md_update_gaugemomenta(inout, delta_p, -1.*eps);
}
void physics::algorithms::md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	::md_update_gaugemomentum(inout, eps, gf, phi, system, kappa, mubar);
}
void physics::algorithms::md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	::md_update_gaugemomentum(inout, eps, gf, phi, system, kappa, mubar);
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

void physics::algorithms::calc_gauge_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf, const hardware::System& system)
{
	gauge_force(gm, gf);
	if(meta::get_use_rectangles(system.get_inputparameters())) {
		gauge_force_tlsym(gm, gf);
	}
}

template<class SPINORFIELD> static void calc_total_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const SPINORFIELD& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	using namespace physics::algorithms;

	force->zero();
	if(!system.get_inputparameters().get_use_gauge_only() )
		calc_fermion_forces(force, gf, phi, system, kappa, mubar);
	calc_gauge_force(force, gf, system);
}

void physics::algorithms::calc_total_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	::calc_total_force(force, gf, phi, system, kappa, mubar);
}
void physics::algorithms::calc_total_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar)
{
	::calc_total_force(force, gf, phi, system, kappa, mubar);
}
