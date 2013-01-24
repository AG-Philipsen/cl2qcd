/** @file
 * Declaration of the molecular dynamics algorithms
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "molecular_dynamics.hpp"

#include <stdexcept>

void physics::algorithms::md_update_gaugemomenta(const physics::lattices::Gaugemomenta * const dest, const physics::lattices::Gaugemomenta& src, const hmc_float eps)
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

void physics::algorithms::md_update_gaugefield(const physics::lattices::Gaugefield * const gf, const physics::lattices::Gaugemomenta& gm, const hmc_float eps)
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
