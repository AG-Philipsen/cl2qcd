/** @file
 * Implementation of the force algorithms
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * (c) 2012-2013 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
 */

#include "forces.hpp"

#include "../../meta/util.hpp"
#include "../../hardware/device.hpp"
#include "../../hardware/code/molecular_dynamics.hpp"

void physics::algorithms::gauge_force(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield& gf)
{
	auto gm_bufs = gm->get_buffers();
	auto gf_bufs = gf.get_buffers();
	size_t num_bufs = gm_bufs.size();
	if(num_bufs != gf_bufs.size()) {
		throw Print_Error_Message("Input buffers seem to use different devices.", __FILE__, __LINE__);
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto gm_buf = gm_bufs[i];
		auto gf_buf = gf_bufs[i];
		auto code = gm_buf->get_device()->get_molecular_dynamics_code();
		code->gauge_force_device(gf_buf, gm_buf);
	}
	gm->update_halo();
}

void physics::algorithms::gauge_force_tlsym(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield& gf)
{
	auto gm_bufs = gm->get_buffers();
	auto gf_bufs = gf.get_buffers();
	size_t num_bufs = gm_bufs.size();
	if(num_bufs != gf_bufs.size()) {
		throw Print_Error_Message("Input buffers seem to use different devices.", __FILE__, __LINE__);
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto gm_buf = gm_bufs[i];
		auto gf_buf = gf_bufs[i];
		auto code = gm_buf->get_device()->get_molecular_dynamics_code();
		code->gauge_force_tlsym_device(gf_buf, gm_buf);
	}
	gm->update_halo();
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
