/** @file
 * Implementation of the physics::lattices::Gaugemomenta class
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
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

#include "gaugemomenta.hpp"
#include "../../meta/util.hpp"
#include <stdexcept>
#include "../../hardware/device.hpp"
#include "../../hardware/code/gaugemomentum.hpp"
#include "../../hardware/buffers/halo_update.hpp"
#include "../../geometry/index.hpp"
#include <cstring>

static std::vector<const hardware::buffers::Gaugemomentum *> allocate_buffers(const hardware::System& system);
static void update_halo_aos(const std::vector<const hardware::buffers::Gaugemomentum *> buffers, const hardware::System& system);
static void update_halo_soa(const std::vector<const hardware::buffers::Gaugemomentum *> buffers, const hardware::System& system);

physics::lattices::Gaugemomenta::Gaugemomenta(const hardware::System& system, const GaugemomentaParametersInterface& parametersInterface)
	: system(system), gaugemomentaParametersInterface(parametersInterface),
	  buffers(allocate_buffers(system))
{
}

static  std::vector<const hardware::buffers::Gaugemomentum *> allocate_buffers(const hardware::System& system)
{
	using hardware::buffers::Gaugemomentum;

	// only use device 0 for now
	auto devices = system.get_devices();
	std::vector<const Gaugemomentum*> buffers;
	buffers.reserve(devices.size());
	for(auto device: devices) {
		buffers.push_back(new Gaugemomentum(NDIM * get_vol4d(device->get_mem_lattice_size()), device));
	}
	return buffers;
}

physics::lattices::Gaugemomenta::~Gaugemomenta()
{
    for(auto buffer: buffers) {
		delete buffer;
	}
}

const std::vector<const hardware::buffers::Gaugemomentum *> physics::lattices::Gaugemomenta::get_buffers() const noexcept
{
	return buffers;
}

void physics::lattices::Gaugemomenta::zero() const
{
for(auto buffer: buffers) {
		buffer->get_device()->getGaugemomentumCode()->set_zero_gaugemomentum(buffer);
	}
}

void physics::lattices::Gaugemomenta::gaussian(const physics::PRNG& prng) const
{
	size_t num_bufs = buffers.size();
	auto prng_bufs = prng.get_buffers();
	if(num_bufs != prng_bufs.size()) {
		throw std::invalid_argument("The PRNG is using different devices than the gaugemomenta");
	}
	for(size_t i = 0; i < num_bufs; ++i) {
		auto buf = buffers[i];
		buf->get_device()->getGaugemomentumCode()->generate_gaussian_gaugemomenta_device(buf, prng_bufs[i]);
	}
	update_halo();
}

void physics::lattices::saxpy(const Gaugemomenta* out, const hmc_float alpha, const Gaugemomenta& x)
{
	saxpy(out, alpha, x, *out);
}

void physics::lattices::saxpy(const Gaugemomenta* out, const Scalar<hmc_float>& alpha, const Gaugemomenta& x)
{
	saxpy(out, alpha, x, *out);
}

void physics::lattices::saxpy(const Gaugemomenta* out, const hmc_float alpha, const Gaugemomenta& x, const Gaugemomenta& y)
{
	const Scalar<hmc_float> alpha_buf(out->system);
	alpha_buf.store(alpha);
	saxpy(out, alpha_buf, x, y);
}

void physics::lattices::saxpy(const Gaugemomenta* out, const Scalar<hmc_float>& alpha, const Gaugemomenta& x, const Gaugemomenta& y)
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
		device->getGaugemomentumCode()->saxpy_device(x_bufs[i], y_bufs[i], alpha_bufs[i], out_buf);
	}
}

hmc_float physics::lattices::squarenorm(const Gaugemomenta& field)
{
	const Scalar<hmc_float> res(field.system);
	squarenorm(&res, field);
	return res.get();
}

void physics::lattices::squarenorm(const Scalar<hmc_float>* res, const Gaugemomenta& field)
{
	auto field_buffers = field.get_buffers();
	auto res_buffers = res->get_buffers();
	size_t num_buffers = field_buffers.size();

	if(num_buffers != res_buffers.size()) {
		throw std::invalid_argument("The given lattices do not sue the same number of devices.");
	}

	for(size_t i = 0; i < num_buffers; ++i) {
		auto field_buf = field_buffers[i];
		auto res_buf = res_buffers[i];
		auto device = field_buf->get_device();
		auto code = device->getGaugemomentumCode();

		code->set_float_to_gaugemomentum_squarenorm_device(field_buf, res_buf);
	}

	res->sum();
}

template<> size_t physics::lattices::get_flops<physics::lattices::Gaugemomenta, physics::lattices::saxpy>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto gaugemomentum_code = devices[0]->getGaugemomentumCode();
	return gaugemomentum_code->get_flop_size("gaugemomentum_saxpy");
}

template<> size_t physics::lattices::get_flops<physics::lattices::Gaugemomenta, physics::lattices::squarenorm>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto gaugemomentum_code = devices[0]->getGaugemomentumCode();
	return gaugemomentum_code->get_flop_size("gaugemomentum_squarenorm");
}


void physics::lattices::log_squarenorm(const std::string& msg, const physics::lattices::Gaugemomenta& x)
{
	if(logger.beDebug()) {
		hmc_float tmp = squarenorm(x);
		logger.debug() << msg << std::scientific << std::setprecision(10) << tmp;
	}
}

void physics::lattices::Gaugemomenta::update_halo() const
{
	if(buffers.size() > 1) { // for a single device this will be a noop
		// currently either all or none of the buffers must be SOA
		if(buffers[0]->is_soa()) {
			update_halo_soa(buffers, system);
		} else {
			update_halo_aos(buffers, system);
		}
	}
}

static void update_halo_aos(const std::vector<const hardware::buffers::Gaugemomentum *> buffers, const hardware::System& system)
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo<ae>(buffers, system, NDIM);
}

static void update_halo_soa(const std::vector<const hardware::buffers::Gaugemomentum *> buffers, const hardware::System& system)
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(!buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo_soa<ae>(buffers, system, .5, 2 * NDIM);
}

unsigned physics::lattices::Gaugemomenta::get_elements() const noexcept
{
    return gaugemomentaParametersInterface.getNumberOfElements();
}

void physics::lattices::Gaugemomenta::import(const ae * const host) const
{
	logger.trace() << "importing gaugemomenta";
	if(buffers.size() == 1) {
		auto device = buffers[0]->get_device();
		device->getGaugemomentumCode()->importGaugemomentumBuffer(buffers[0], host);
	} else {
		auto const _device = buffers.at(0)->get_device();
		auto const local_size = _device->get_local_lattice_size();
		size_4 const halo_size(local_size.x, local_size.y, local_size.z, _device->get_halo_size());
		auto const grid_size = _device->getGridSize();
		if(grid_size.x != 1 || grid_size.y != 1 || grid_size.z != 1) {
			throw Print_Error_Message("Not implemented!", __FILE__, __LINE__);
		}
		for(auto const buffer: buffers) {
			auto device = buffer->get_device();
			ae * mem_host = new ae[buffer->get_elements()];

			size_4 offset(0, 0, 0, device->getGridPos().t * local_size.t);
			logger.debug() << offset;
			const size_t local_volume = get_vol4d(local_size) * NDIM;
			memcpy(mem_host, &host[uint(LinkIndex(Index(offset, LatticeExtents(gaugemomentaParametersInterface.getNs(),gaugemomentaParametersInterface.getNt())),TDIR))], local_volume * sizeof(ae));

			const size_t halo_volume = get_vol4d(halo_size) * NDIM;
			size_4 halo_offset(0, 0, 0, (offset.t + local_size.t) % gaugemomentaParametersInterface.getNt());
			logger.debug() << halo_offset;
			memcpy(&mem_host[local_volume], &host[uint(LinkIndex(Index(halo_offset, LatticeExtents(gaugemomentaParametersInterface.getNs(),gaugemomentaParametersInterface.getNt())),TDIR))], halo_volume * sizeof(ae));

			halo_offset = size_4(0, 0, 0, (offset.t + gaugemomentaParametersInterface.getNt() - halo_size.t) % gaugemomentaParametersInterface.getNt());
			logger.debug() << halo_offset;
			memcpy(&mem_host[local_volume + halo_volume], &host[uint(LinkIndex(Index(halo_offset, LatticeExtents(gaugemomentaParametersInterface.getNs(),gaugemomentaParametersInterface.getNt())),TDIR))], halo_volume * sizeof(ae));

			device->getGaugemomentumCode()->importGaugemomentumBuffer(buffer, mem_host);

			delete[] mem_host;
		}
	}
	logger.trace() << "import complete";
}
