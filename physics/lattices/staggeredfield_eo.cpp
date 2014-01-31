/** @file
 * Implementation of the physics::lattices::Staggeredfield_eo class
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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


#include "staggeredfield_eo.hpp"
#include "../../meta/util.hpp"
#include <cassert>
#include "../../hardware/code/spinors_staggered.hpp"
//#include "../../hardware/code/fermions_staggered.hpp"
#include "../../meta/type_ops.hpp"
#include "../../hardware/buffers/halo_update.hpp"
//For hardware::code::get_eoprec_spinorfieldsize()
#include "../../hardware/code/spinors.hpp"

static std::vector<const hardware::buffers::SU3vec *> allocate_buffers(const hardware::System& system);

/*To be added...
static void update_halo_soa(const std::vector<const hardware::buffers::Spinor *> buffers, const meta::Inputparameters& params);
static void update_halo_aos(const std::vector<const hardware::buffers::Spinor *> buffers, const meta::Inputparameters& params);
static void extract_boundary(char* host, const hardware::buffers::Spinor * buffer, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS);
static void send_halo(const hardware::buffers::Spinor * buffer, const char* host, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS);
*/

physics::lattices::Staggeredfield_eo::Staggeredfield_eo(const hardware::System& system)
	: system(system), buffers(allocate_buffers(system))
{
}

static std::vector<const hardware::buffers::SU3vec *> allocate_buffers(const hardware::System& system)
{
	using hardware::buffers::SU3vec;

	auto devices = system.get_devices();
	std::vector<const SU3vec*> buffers;
	buffers.reserve(devices.size());
	for(auto device: devices) {
		buffers.push_back(new SU3vec(hardware::code::get_eoprec_spinorfieldsize(device->get_mem_lattice_size()), device));
	}
	return buffers;
}

physics::lattices::Staggeredfield_eo::~Staggeredfield_eo()
{
for(auto buffer: buffers) {
		delete buffer;
	}
}

const std::vector<const hardware::buffers::SU3vec *> physics::lattices::Staggeredfield_eo::get_buffers() const noexcept
{
	return buffers;
}

hmc_complex physics::lattices::scalar_product(const Staggeredfield_eo& left, const Staggeredfield_eo& right)
{
	const Scalar<hmc_complex> res(left.system);
	scalar_product(&res, left, right);
	return res.get();
}

void physics::lattices::scalar_product(const Scalar<hmc_complex>* res, const Staggeredfield_eo& left, const Staggeredfield_eo& right)
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
		auto spinor_code = device->get_spinor_staggered_code();

		spinor_code->set_complex_to_scalar_product_eoprec_device(left_buf, right_buf, res_buf);
	}
	res->sum();
}

hmc_float physics::lattices::squarenorm(const Staggeredfield_eo& field)
{
	const Scalar<hmc_float> res(field.system);
	squarenorm(&res, field);
	return res.get();
}

void physics::lattices::squarenorm(const Scalar<hmc_float>* res, const Staggeredfield_eo& field)
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
		auto spinor_code = device->get_spinor_staggered_code();

		spinor_code->set_float_to_global_squarenorm_eoprec_device(field_buf, res_buf);
	}
	res->sum();
}

void physics::lattices::Staggeredfield_eo::set_zero() const
{
for(auto buffer: buffers) {
		auto spinor_code = buffer->get_device()->get_spinor_staggered_code();
		spinor_code->set_zero_spinorfield_eoprec_device(buffer);
	}
}

void physics::lattices::Staggeredfield_eo::set_cold() const
{
for(auto buffer: buffers) {
		auto spinor_code = buffer->get_device()->get_spinor_staggered_code();
		spinor_code->set_cold_spinorfield_eoprec_device(buffer);
	}
}

void physics::lattices::Staggeredfield_eo::set_gaussian(const physics::PRNG& prng) const
{
	auto prng_bufs = prng.get_buffers();

	if(buffers.size() != prng_bufs.size()) {
		throw std::invalid_argument("PRNG does not use same devices as spinorfield");
	}

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto spin_buf = buffers[i];
		auto prng_buf = prng_bufs[i];
		spin_buf->get_device()->get_spinor_staggered_code()->set_gaussian_spinorfield_eoprec_device(spin_buf, prng_buf);
	}
	if((system.get_devices()).size()!=1)
		update_halo();
}

void physics::lattices::sax(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x)
{
	auto out_bufs = out->get_buffers();
	auto x_bufs = x.get_buffers();

	if(out_bufs.size() != x_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->get_spinor_staggered_code()->sax_eoprec_device(x_bufs[i], alpha, out_buf);
	}
}

void physics::lattices::sax(const Staggeredfield_eo* out, const Scalar<hmc_complex>& alpha, const Staggeredfield_eo& x)
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
		device->get_spinor_staggered_code()->sax_eoprec_device(x_bufs[i], alpha_bufs[i], out_buf);
	}
}

void physics::lattices::saxpy(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y)
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
		device->get_spinor_staggered_code()->saxpy_eoprec_device(x_bufs[i], y_bufs[i], alpha, out_buf);
	}
}

void physics::lattices::saxpy(const Staggeredfield_eo* out, const Scalar<hmc_complex>& alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y)
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
		device->get_spinor_staggered_code()->saxpy_eoprec_device(x_bufs[i], y_bufs[i], alpha_bufs[i], out_buf);
	}
}

void physics::lattices::saxpy(const Staggeredfield_eo* out, const hmc_float alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y)
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
		device->get_spinor_staggered_code()->saxpy_eoprec_device(x_bufs[i], y_bufs[i], alpha, out_buf);
	}
}

void physics::lattices::saxpy(const Staggeredfield_eo* out, const Scalar<hmc_float>& alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y)
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
		device->get_spinor_staggered_code()->saxpy_eoprec_device(x_bufs[i], y_bufs[i], alpha_bufs[i], out_buf);
	}
}

void physics::lattices::saxpby(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x, const hmc_complex beta, const Staggeredfield_eo& y)
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
		device->get_spinor_staggered_code()->saxpby_eoprec_device(x_bufs[i], y_bufs[i], alpha, beta, out_buf);
	}
}

void physics::lattices::saxpby(const Staggeredfield_eo* out, const Scalar<hmc_complex>& alpha, const Staggeredfield_eo& x, const Scalar<hmc_complex>& beta, const Staggeredfield_eo& y)
{
	auto out_bufs = out->get_buffers();
	auto alpha_bufs = alpha.get_buffers();
	auto x_bufs = x.get_buffers();
	auto beta_bufs = beta.get_buffers();
	auto y_bufs = y.get_buffers();
	
	if(out_bufs.size() != alpha_bufs.size() || out_bufs.size() != beta_bufs.size() || out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->get_spinor_staggered_code()->saxpby_eoprec_device(x_bufs[i], y_bufs[i], alpha_bufs[i], beta_bufs[i], out_buf);
	}
}

void physics::lattices::saxpby(const Staggeredfield_eo* out, const hmc_float alpha, const Staggeredfield_eo& x, const hmc_float beta, const Staggeredfield_eo& y)
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
		device->get_spinor_staggered_code()->saxpby_eoprec_device(x_bufs[i], y_bufs[i], alpha, beta, out_buf);
	}
}

void physics::lattices::saxpby(const Staggeredfield_eo* out, const Scalar<hmc_float>& alpha, const Staggeredfield_eo& x, const Scalar<hmc_float>& beta, const Staggeredfield_eo& y)
{
	auto out_bufs = out->get_buffers();
	auto alpha_bufs = alpha.get_buffers();
	auto x_bufs = x.get_buffers();
	auto beta_bufs = beta.get_buffers();
	auto y_bufs = y.get_buffers();
	
	if(out_bufs.size() != alpha_bufs.size() || out_bufs.size() != beta_bufs.size() || out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->get_spinor_staggered_code()->saxpby_eoprec_device(x_bufs[i], y_bufs[i], alpha_bufs[i], beta_bufs[i], out_buf);
	}
}

void physics::lattices::saxpbypz(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x, const hmc_complex beta, const Staggeredfield_eo& y, const Staggeredfield_eo& z)
{
	auto out_bufs = out->get_buffers();
	auto x_bufs = x.get_buffers();
	auto y_bufs = y.get_buffers();
	auto z_bufs = z.get_buffers();

	if(out_bufs.size() != x_bufs.size() || out_bufs.size() != y_bufs.size() || out_bufs.size() != z_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		device->get_spinor_staggered_code()->saxpbypz_eoprec_device(x_bufs[i], y_bufs[i], z_bufs[i], alpha, beta, out_buf);
	}
}

void physics::lattices::saxpbypz(const Staggeredfield_eo* out, const Scalar<hmc_complex>& alpha, const Staggeredfield_eo& x, const Scalar<hmc_complex>& beta, const Staggeredfield_eo& y, const Staggeredfield_eo& z)
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
		device->get_spinor_staggered_code()->saxpbypz_eoprec_device(x_bufs[i], y_bufs[i], z_bufs[i], alpha_bufs[i], beta_bufs[i], out_buf);
	}
}

void physics::lattices::sax_vec_and_squarenorm(const Vector<hmc_float>* res, const Vector<hmc_float>& alpha, const Staggeredfield_eo& x)
{
	auto x_buffers = x.get_buffers();
	auto alpha_buffers = alpha.get_buffers();
	auto res_buffers = res->get_buffers();
	size_t num_buffers = x_buffers.size();

	if(num_buffers != res_buffers.size() || num_buffers != alpha_buffers.size()) {
		throw std::invalid_argument("The given lattices do not use the same number of devices.");
	}

	for(size_t i = 0; i < num_buffers; ++i) {
		auto x_buf = x_buffers[i];
		auto alpha_buf = alpha_buffers[i];
		auto res_buf = res_buffers[i];
		auto device = x_buf->get_device();
		auto spinor_code = device->get_spinor_staggered_code();

		spinor_code->sax_vectorized_and_squarenorm_eoprec_device(x_buf, alpha_buf, alpha.get_vector_size(), res_buf);
	}
	//res->sum();
}


template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, physics::lattices::scalar_product>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("scalar_product_staggered_eoprec");
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, physics::lattices::squarenorm>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("global_squarenorm_staggered_eoprec");
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, physics::lattices::sax>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("sax_cplx_stagg_eoprec"); //the same as sax_cplx_arg_stagg_eoprec
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_complex, physics::lattices::saxpy>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("saxpy_cplx_stagg_eoprec"); //the same as saxpy_cplx_arg_stagg_eoprec
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::saxpy>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("saxpy_real_stagg_eoprec"); //the same as saxpy_real_arg_stagg_eoprec
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_complex, physics::lattices::saxpby>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("saxpby_cplx_stagg_eoprec"); //the same as saxpby_cplx_arg_stagg_eoprec
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::saxpby>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("saxpby_real_stagg_eoprec"); //the same as saxpby_cplx_arg_stagg_eoprec
}

template<> size_t physics::lattices::get_flops<physics::lattices::Staggeredfield_eo, physics::lattices::saxpbypz>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto spinor_code = devices[0]->get_spinor_staggered_code();
	return spinor_code->get_flop_size("saxpbypz_cplx_stagg_eoprec");//the same as saxpbypz_cplx_arg_stagg_eoprec
}

void physics::lattices::log_squarenorm(const std::string& msg, const physics::lattices::Staggeredfield_eo& x)
{
	if(logger.beDebug()) {
		hmc_float tmp = squarenorm(x);
		logger.debug() << msg << std::scientific << std::setprecision(10) << tmp;
	}
}

void physics::lattices::Staggeredfield_eo::update_halo() const
{
	//So far, this is simple a throw
	throw Print_Error_Message("Halo update is not implemented, yet.", __FILE__, __LINE__);
	/*
	if(buffers.size() > 1) { // for a single device this will be a noop
		// currently either all or none of the buffers must be SOA
		if(buffers[0]->is_soa()) {
			update_halo_soa(buffers, system.get_inputparameters());
		} else {
			update_halo_aos(buffers, system.get_inputparameters());
		}
	}
	*/
}

//This code is useful for using the template pseudo_randomize defined in lattices/util.hpp
void physics::lattices::Staggeredfield_eo::import(const su3vec * const host) const
{
	logger.trace() << "importing staggeredfield_eo";
	if(buffers.size() == 1) {
		buffers[0]->load(host);
	} else {
		throw Print_Error_Message("Import not implemented for multi device staggeredfield_eo!", __FILE__, __LINE__);
		/*
		auto params = system.get_inputparameters();
		auto const _device = buffers.at(0)->get_device();
		auto const local_size = _device->get_local_lattice_size();
		size_4 const halo_size(local_size.x, local_size.y, local_size.z, _device->get_halo_size());
		auto const grid_size = _device->get_grid_size();
		if(grid_size.x != 1 || grid_size.y != 1 || grid_size.z != 1) {
			throw Print_Error_Message("Not implemented!", __FILE__, __LINE__);
		}
		for(auto const buffer: buffers) {
			auto device = buffer->get_device();

			size_4 offset(0, 0, 0, device->get_grid_pos().t * local_size.t);
			logger.debug() << offset;
			const size_t local_volume = get_vol4d(local_size);
			buffer->load(&host[get_global_pos(offset, params)], local_volume);

			const size_t halo_volume = get_vol4d(halo_size);
			size_4 halo_offset(0, 0, 0, (offset.t + local_size.t) % params.get_ntime());
			logger.debug() << halo_offset;
			logger.trace() << get_global_pos(halo_offset, params);
			logger.trace() << halo_volume;
			logger.trace() << get_elements();
			assert(get_global_pos(halo_offset, params) + halo_volume <= get_elements());
			buffer->load(&host[get_global_pos(halo_offset, params)], halo_volume, local_volume);

			halo_offset = size_4(0, 0, 0, (offset.t + params.get_ntime() - halo_size.t) % params.get_ntime());
			logger.debug() << halo_offset;
			assert(get_global_pos(halo_offset, params) + halo_volume <= get_elements());
			buffer->load(&host[get_global_pos(halo_offset, params)], halo_volume, local_volume + halo_volume);
		}
		*/
	}
	logger.trace() << "import complete";
}

unsigned physics::lattices::Staggeredfield_eo::get_elements() const noexcept
{
	return hardware::code::get_spinorfieldsize(system.get_inputparameters());
}



//To be added (so far only single GPU and only EO preconditioning)...
#if 0

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
		auto spinor_code = in_bufs[i]->get_device()->get_spinor_code();
		spinor_code->convert_to_eoprec_device(even_bufs[i], odd_bufs[i], in_bufs[i]);
	}
}

void physics::lattices::convert_from_eoprec(const Spinorfield* merged, const Spinorfield_eo& even, const Spinorfield_eo& odd)
{
	auto merged_bufs = merged->get_buffers();
	auto even_bufs = even.get_buffers();
	auto odd_bufs = odd.get_buffers();
	size_t num_bufs = merged_bufs.size();

	if(even_bufs.size() != num_bufs || odd_bufs.size() != num_bufs) {
		throw std::invalid_argument("Output buffers do not use the same devices as input buffers");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto spinor_code = merged_bufs[i]->get_device()->get_spinor_code();
		spinor_code->convert_from_eoprec_device(even_bufs[i], odd_bufs[i], merged_bufs[i]);
	}
}

static void update_halo_aos(const std::vector<const hardware::buffers::Spinor *> buffers, const meta::Inputparameters& params)
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo<spinor>(buffers, params, .5 /* only even or odd sites */ );
}

static void update_halo_soa(const std::vector<const hardware::buffers::Spinor *> buffers, const meta::Inputparameters& params)
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(!buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo_soa<spinor>(buffers, params, .5 /* only even or odd sites */ );
}
#endif
