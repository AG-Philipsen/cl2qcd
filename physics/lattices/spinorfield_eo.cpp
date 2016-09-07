/** @file
 * Implementation of the physics::lattices::Spinorfield_eo class
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

#include "spinorfield_eo.hpp"
#include <cassert>
#include "../../hardware/code/spinors.hpp"
#include "../../hardware/code/fermions.hpp"
#include "../../hardware/buffers/halo_update.hpp"
#include "util.hpp"
#include <algorithm>

physics::lattices::Spinorfield_eo::Spinorfield_eo(const hardware::System& system, const SpinorfieldEoParametersInterface& spinorfieldEoParametersInterface)
	: system(system), spinorfieldEo(system), spinorfieldEoParametersInterface(spinorfieldEoParametersInterface)
#ifdef LAZY_HALO_UPDATES
	  , valid_halo_width(0)
#endif
{
}

physics::lattices::Spinorfield_eo::~Spinorfield_eo()
{
}

const std::vector<const hardware::buffers::Spinor *> physics::lattices::Spinorfield_eo::get_buffers() const noexcept
{
	return spinorfieldEo.get_buffers();
}

hmc_complex physics::lattices::scalar_product(const Spinorfield_eo& left, const Spinorfield_eo& right)
{
	const Scalar<hmc_complex> res(left.system);
	scalar_product(&res, left, right);
	return res.get();
}

hmc_complex physics::lattices::scalar_product(const Spinorfield_eo& left, const Spinorfield_eo& right, const Scalar<hmc_complex>* res)
{
	auto res_buffers = res->get_buffers();
	auto left_buffers = left.get_buffers();
	auto right_buffers = right.get_buffers();
	size_t num_buffers = res_buffers.size();

	if(num_buffers != left_buffers.size() || num_buffers != right_buffers.size()) {
		throw std::invalid_argument("The given lattices do not use the same number of devices.");
	}

	for(size_t i = 0; i < num_buffers; ++i) {
		auto res_buf = res_buffers[i];
		auto left_buf = left_buffers[i];
		auto right_buf = right_buffers[i];
		auto device = res_buf->get_device();
		auto spinor_code = device->getSpinorCode();

		spinor_code->set_complex_to_scalar_product_eoprec_device(left_buf, right_buf, res_buf);
	}
	return res->get_sum();
}

void physics::lattices::scalar_product(const Scalar<hmc_complex>* res, const Spinorfield_eo& left, const Spinorfield_eo& right)
{
	auto res_buffers = res->get_buffers();
	auto left_buffers = left.get_buffers();
	auto right_buffers = right.get_buffers();
	size_t num_buffers = res_buffers.size();

	if(num_buffers != left_buffers.size() || num_buffers != right_buffers.size()) {
		throw std::invalid_argument("The given lattices do not use the same number of devices.");
	}

	for(size_t i = 0; i < num_buffers; ++i) {
		auto res_buf = res_buffers[i];
		auto left_buf = left_buffers[i];
		auto right_buf = right_buffers[i];
		auto device = res_buf->get_device();
		auto spinor_code = device->getSpinorCode();

		spinor_code->set_complex_to_scalar_product_eoprec_device(left_buf, right_buf, res_buf);
	}
	res->sum();
}

hmc_float physics::lattices::squarenorm(const Spinorfield_eo& field)
{
	const Scalar<hmc_float> res(field.system);
	squarenorm(&res, field);
	return res.get();
}

hmc_float physics::lattices::squarenorm(const Spinorfield_eo& field, const Scalar<hmc_float>* res)
{
	auto field_buffers = field.get_buffers();
	auto res_buffers = res->get_buffers();
	size_t num_buffers = field_buffers.size();

	if(num_buffers != res_buffers.size()) {
		throw std::invalid_argument("The given lattices do not use the same number of devices.");
	}

	for(size_t i = 0; i < num_buffers; ++i) {
		auto field_buf = field_buffers[i];
		auto res_buf = res_buffers[i];
		auto device = field_buf->get_device();
		auto spinor_code = device->getSpinorCode();

		spinor_code->set_float_to_global_squarenorm_eoprec_device(field_buf, res_buf);
	}
	return res->get_sum();
}

void physics::lattices::squarenorm(const Scalar<hmc_float>* res, const Spinorfield_eo& field)
{
	auto field_buffers = field.get_buffers();
	auto res_buffers = res->get_buffers();
	size_t num_buffers = field_buffers.size();

	if(num_buffers != res_buffers.size()) {
		throw std::invalid_argument("The given lattices do not use the same number of devices.");
	}

	for(size_t i = 0; i < num_buffers; ++i) {
		auto field_buf = field_buffers[i];
		auto res_buf = res_buffers[i];
		auto device = field_buf->get_device();
		auto spinor_code = device->getSpinorCode();

		spinor_code->set_float_to_global_squarenorm_eoprec_device(field_buf, res_buf);
	}
	res->sum();
}

void physics::lattices::Spinorfield_eo::zero() const
{
for(auto buffer: spinorfieldEo.get_buffers()) {
		auto spinor_code = buffer->get_device()->getSpinorCode();
		spinor_code->set_zero_spinorfield_eoprec_device(buffer);
	}
	mark_halo_clean();
}

void physics::lattices::Spinorfield_eo::cold() const
{
for(auto buffer: spinorfieldEo.get_buffers()) {
		auto spinor_code = buffer->get_device()->getSpinorCode();
		spinor_code->set_eoprec_spinorfield_cold_device(buffer);
	}
	mark_halo_clean();
}

void physics::lattices::Spinorfield_eo::gaussian(const physics::PRNG& prng) const
{
	auto prng_bufs = prng.get_buffers();

	if(spinorfieldEo.get_buffers().size() != prng_bufs.size()) {
		throw std::invalid_argument("PRNG does not use same devices as spinorfield");
	}

	for(size_t i = 0; i < spinorfieldEo.get_buffers().size(); ++i) {
		auto spin_buf = spinorfieldEo.get_buffers()[i];
		auto prng_buf = prng_bufs[i];
		spin_buf->get_device()->getSpinorCode()->generate_gaussian_spinorfield_eo_device(spin_buf, prng_buf);
	}
	mark_halo_dirty();
}

void physics::lattices::saxpy(const Spinorfield_eo* out, const hmc_complex alpha, const Spinorfield_eo& x, const Spinorfield_eo& y)
{
	auto out_bufs = out->get_buffers();
	auto x_bufs = x.get_buffers();
	auto y_bufs = y.get_buffers();

	if(out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->getSpinorCode()->saxpy_eoprec_device(x_bufs[i], y_bufs[i], alpha, out_buf);
	}

	auto const valid_halo_width = std::min(x.get_valid_halo_width(), y.get_valid_halo_width());
	if(valid_halo_width) {
		out->mark_halo_clean(valid_halo_width);
	} else {
		out->mark_halo_dirty();
	}
}

void physics::lattices::saxpy_AND_squarenorm(const Spinorfield_eo* out, const Scalar<hmc_complex>& alpha, const Spinorfield_eo& x, const Spinorfield_eo& y, const Scalar<hmc_complex>& squarenorm)
{
	auto out_bufs = out->get_buffers();
	auto x_bufs = x.get_buffers();
	auto y_bufs = y.get_buffers();
	auto alpha_bufs = alpha.get_buffers();
	auto squarenorm_bufs = squarenorm.get_buffers();

	if(out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->getSpinorCode()->saxpy_AND_squarenorm_eo_device(x_bufs[i], y_bufs[i], alpha_bufs[i], out_buf, squarenorm_bufs[i]);
	}

	auto const valid_halo_width = std::min(x.get_valid_halo_width(), y.get_valid_halo_width());
	if(valid_halo_width) {
		out->mark_halo_clean(valid_halo_width);
	} else {
		out->mark_halo_dirty();
	}
}

void physics::lattices::saxpy(const Spinorfield_eo* out, const Scalar<hmc_complex>& alpha, const Spinorfield_eo& x, const Spinorfield_eo& y)
{
	auto out_bufs = out->get_buffers();
	auto alpha_bufs = alpha.get_buffers();
	auto x_bufs = x.get_buffers();
	auto y_bufs = y.get_buffers();

	if(out_bufs.size() != alpha_bufs.size() || out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->getSpinorCode()->saxpy_eoprec_device(x_bufs[i], y_bufs[i], alpha_bufs[i], out_buf);
	}

	auto const valid_halo_width = std::min(x.get_valid_halo_width(), y.get_valid_halo_width());
	if(valid_halo_width) {
		out->mark_halo_clean(valid_halo_width);
	} else {
		out->mark_halo_dirty();
	}
}

void physics::lattices::saxpy_AND_gamma5_eo(const Spinorfield_eo* out, const hmc_complex alpha, const Spinorfield_eo& x, const Spinorfield_eo& y )
{
	logger.trace() << "calling saxpy and g5";
	auto x_bufs = x.get_buffers();
	auto y_bufs = y.get_buffers();
	auto out_bufs = out->get_buffers();

	if(out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->getFermionCode()->saxpy_AND_gamma5_eo_device(x_bufs[i], y_bufs[i], alpha, out_buf);
	}

	auto const valid_halo_width = std::min(x.get_valid_halo_width(), y.get_valid_halo_width());
	if(valid_halo_width) {
		out->mark_halo_clean(valid_halo_width);
	} else {
		out->mark_halo_dirty();
	}
}

void physics::lattices::sax(const Spinorfield_eo* out, const hmc_complex alpha, const Spinorfield_eo& x)
{
	const Scalar<hmc_complex> alpha_buf(out->system);
	alpha_buf.store(alpha);
	sax(out, alpha_buf, x);
}

void physics::lattices::sax(const Spinorfield_eo* out, const Scalar<hmc_complex>& alpha, const Spinorfield_eo& x)
{
	auto out_bufs = out->get_buffers();
	auto alpha_bufs = alpha.get_buffers();
	auto x_bufs = x.get_buffers();

	if(out_bufs.size() != alpha_bufs.size() || out_bufs.size() != x_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->getSpinorCode()->sax_eoprec_device(x_bufs[i], alpha_bufs[i], out_buf);
	}

	auto const valid_halo_width = x.get_valid_halo_width();
	if(valid_halo_width) {
		out->mark_halo_clean(valid_halo_width);
	} else {
		out->mark_halo_dirty();
	}
}

void physics::lattices::saxsbypz(const Spinorfield_eo* out, const hmc_complex alpha, const Spinorfield_eo& x, const hmc_complex beta, const Spinorfield_eo& y, const Spinorfield_eo& z)
{
	const Scalar<hmc_complex> alpha_buf(out->system);
	const Scalar<hmc_complex> beta_buf(out->system);
	alpha_buf.store(alpha);
	beta_buf.store(beta);
	saxsbypz(out, alpha_buf, x, beta_buf, y, z);
}

void physics::lattices::saxsbypz(const Spinorfield_eo* out, const Scalar<hmc_complex>& alpha, const Spinorfield_eo& x, const Scalar<hmc_complex>& beta, const Spinorfield_eo& y, const Spinorfield_eo& z)
{
	auto out_bufs = out->get_buffers();
	auto alpha_bufs = alpha.get_buffers();
	auto x_bufs = x.get_buffers();
	auto beta_bufs = beta.get_buffers();
	auto y_bufs = y.get_buffers();
	auto z_bufs = z.get_buffers();

	if(out_bufs.size() != alpha_bufs.size() || out_bufs.size() != beta_bufs.size() || out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size() || out_bufs.size() != z_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}


	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->getSpinorCode()->saxsbypz_eoprec_device(x_bufs[i], y_bufs[i], z_bufs[i], alpha_bufs[i], beta_bufs[i], out_buf);
	}

	auto const valid_halo_width = std::min({x.get_valid_halo_width(), y.get_valid_halo_width(), z.get_valid_halo_width()});
	if(valid_halo_width) {
		out->mark_halo_clean(valid_halo_width);
	} else {
		out->mark_halo_dirty();
	}
}

void physics::lattices::convert_to_eoprec(const Spinorfield_eo* even, const Spinorfield_eo* odd, const Spinorfield& in)
{
	auto even_bufs = even->get_buffers();
	auto odd_bufs = odd->get_buffers();
	auto in_bufs = in.get_buffers();
	size_t num_bufs = in_bufs.size();

	if(even_bufs.size() != num_bufs || odd_bufs.size() != num_bufs) {
		throw std::invalid_argument("Output buffers do not use the same devices as input buffers");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto spinor_code = in_bufs[i]->get_device()->getSpinorCode();
		spinor_code->convert_to_eoprec_device(even_bufs[i], odd_bufs[i], in_bufs[i]);
	}
	even->mark_halo_clean();
	odd->mark_halo_clean();
}

void physics::lattices::convert_from_eoprec(const Spinorfield* merged, const Spinorfield_eo& even, const Spinorfield_eo& odd)
{
	even.require_halo();
	odd.require_halo();
	auto merged_bufs = merged->get_buffers();
	auto even_bufs = even.get_buffers();
	auto odd_bufs = odd.get_buffers();
	size_t num_bufs = merged_bufs.size();

	if(even_bufs.size() != num_bufs || odd_bufs.size() != num_bufs) {
		throw std::invalid_argument("Output buffers do not use the same devices as input buffers");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto spinor_code = merged_bufs[i]->get_device()->getSpinorCode();
		spinor_code->convert_from_eoprec_device(even_bufs[i], odd_bufs[i], merged_bufs[i]);
	}
}

void physics::lattices::Spinorfield_eo::gamma5() const
{
for(auto buffer: spinorfieldEo.get_buffers()) {
		auto fermion_code = buffer->get_device()->getFermionCode();
		fermion_code->gamma5_eo_device(buffer);
	}
}

template<> size_t physics::lattices::get_flops<physics::lattices::Spinorfield_eo, physics::lattices::scalar_product>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->getSpinorCode();
	return spinor_code->get_flop_size("scalar_product_eoprec");
}
template<> size_t physics::lattices::get_read_write_size<physics::lattices::Spinorfield_eo, physics::lattices::scalar_product>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->getSpinorCode();
	return spinor_code->get_read_write_size("scalar_product_eoprec");
}

template<> size_t physics::lattices::get_flops<physics::lattices::Spinorfield_eo, physics::lattices::squarenorm>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->getSpinorCode();
	return spinor_code->get_flop_size("global_squarenorm_eoprec");
}
template<> size_t physics::lattices::get_read_write_size<physics::lattices::Spinorfield_eo, physics::lattices::squarenorm>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->getSpinorCode();
	return spinor_code->get_read_write_size("global_squarenorm_eoprec");
}

template<> size_t physics::lattices::get_flops<physics::lattices::Spinorfield_eo, physics::lattices::sax>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->getSpinorCode();
	return spinor_code->get_read_write_size("sax_eoprec");
}
template<> size_t physics::lattices::get_read_write_size<physics::lattices::Spinorfield_eo, physics::lattices::sax>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->getSpinorCode();
	return spinor_code->get_read_write_size("sax_eoprec");
}

template<> size_t physics::lattices::get_flops<physics::lattices::Spinorfield_eo, physics::lattices::saxpy>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->getSpinorCode();
	return spinor_code->get_flop_size("saxpy_eoprec");
}
template<> size_t physics::lattices::get_read_write_size<physics::lattices::Spinorfield_eo, physics::lattices::saxpy>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->getSpinorCode();
	return spinor_code->get_read_write_size("saxpy_eoprec");
}

template<> size_t physics::lattices::get_flops<physics::lattices::Spinorfield_eo, physics::lattices::saxsbypz>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->getSpinorCode();
	return spinor_code->get_flop_size("saxsbypz_eoprec");
}
template<> size_t physics::lattices::get_read_write_size<physics::lattices::Spinorfield_eo, physics::lattices::saxsbypz>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->getSpinorCode();
	return spinor_code->get_read_write_size("saxsbypz_eoprec");
}

void physics::lattices::log_squarenorm(const std::string& msg, const physics::lattices::Spinorfield_eo& x)
{
	if(logger.beDebug()) {
		hmc_float tmp = squarenorm(x);
		logger.debug() << msg << std::scientific << std::setprecision(10) << tmp;
	}
}

void physics::lattices::Spinorfield_eo::mark_halo_dirty() const
{
	spinorfieldEo.mark_halo_dirty();
}

void physics::lattices::Spinorfield_eo::require_halo(unsigned reqd_width) const
{
	spinorfieldEo.require_halo(reqd_width);
}

hardware::lattices::Spinorfield_eoHaloUpdate physics::lattices::Spinorfield_eo::require_halo_async(unsigned reqd_width) const
{
	return spinorfieldEo.require_halo_async(reqd_width);
}

void physics::lattices::Spinorfield_eo::mark_halo_clean(unsigned width) const
{
	spinorfieldEo.mark_halo_clean(width);
}

unsigned physics::lattices::Spinorfield_eo::get_valid_halo_width() const
{
	return valid_halo_width;
}

void physics::lattices::copyData(const physics::lattices::Spinorfield_eo* to, const physics::lattices::Spinorfield_eo& from)
{
	copyData<Spinorfield_eo>(to, from);

	auto const valid_halo_width = from.get_valid_halo_width();
	if(valid_halo_width) {
		to->mark_halo_clean(valid_halo_width);
	} else {
		to->mark_halo_dirty();
	}
}

void physics::lattices::copyData(const physics::lattices::Spinorfield_eo* to, const physics::lattices::Spinorfield_eo* from)
{
	copyData(to, *from);
}
