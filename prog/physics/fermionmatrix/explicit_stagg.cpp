/** @file
 * Implementation of explicit fermionamtrix operations
 */

#include "fermionmatrix.hpp"


#if 0
void physics::fermionmatrix::M_wilson(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa)
{
	auto out_bufs = out->get_buffers();
	auto gf_bufs = gf.get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != gf_bufs.size() || num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->get_fermion_code();
		fermion_code->M_wilson_device(in_bufs[i], out_bufs[i], gf_bufs[i], kappa);
	}
	out->update_halo();
}

void physics::fermionmatrix::M_tm_plus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa, hmc_float mubar)
{
	auto out_bufs = out->get_buffers();
	auto gf_bufs = gf.get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != gf_bufs.size() || num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->get_fermion_code();
		fermion_code->M_tm_plus_device(in_bufs[i], out_bufs[i], gf_bufs[i], kappa, mubar);
	}
	out->update_halo();
}

void physics::fermionmatrix::M_tm_minus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa, hmc_float mubar)
{
	auto out_bufs = out->get_buffers();
	auto gf_bufs = gf.get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != gf_bufs.size() || num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->get_fermion_code();
		fermion_code->M_tm_minus_device(in_bufs[i], out_bufs[i], gf_bufs[i], kappa, mubar);
	}
	out->update_halo();
}

void physics::fermionmatrix::M_tm_inverse_sitediagonal(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar)
{
	auto out_bufs = out->get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->get_fermion_code();
		fermion_code->M_tm_inverse_sitediagonal_device(in_bufs[i], out_bufs[i], mubar);
	}
}

void physics::fermionmatrix::M_tm_sitediagonal(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar)
{
	auto out_bufs = out->get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->get_fermion_code();
		fermion_code->M_tm_sitediagonal_device(in_bufs[i], out_bufs[i], mubar);
	}
}

void physics::fermionmatrix::M_tm_inverse_sitediagonal_minus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar)
{
	auto out_bufs = out->get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->get_fermion_code();
		fermion_code->M_tm_inverse_sitediagonal_minus_device(in_bufs[i], out_bufs[i], mubar);
	}
}

void physics::fermionmatrix::M_tm_sitediagonal_minus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar)
{
	auto out_bufs = out->get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->get_fermion_code();
		fermion_code->M_tm_sitediagonal_minus_device(in_bufs[i], out_bufs[i], mubar);
	}
}

void physics::fermionmatrix::dslash(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, int evenodd, hmc_float kappa)
{
	auto out_bufs = out->get_buffers();
	auto gf_bufs = gf.get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != gf_bufs.size() || num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->get_fermion_code();
		fermion_code->dslash_eo_device(in_bufs[i], out_bufs[i], gf_bufs[i], evenodd, kappa);
	}
	out->update_halo();
}

#endif

