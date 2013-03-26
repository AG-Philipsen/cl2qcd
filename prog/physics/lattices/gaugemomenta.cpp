/** @file
 * Implementation of the physics::lattices::Gaugemomenta class
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "gaugemomenta.hpp"
#include "../../meta/util.hpp"
#include <stdexcept>
#include "../../hardware/device.hpp"
#include "../../hardware/code/gaugemomentum.hpp"
#include "../../hardware/buffers/halo_update.hpp"

static std::vector<const hardware::buffers::Gaugemomentum *> allocate_buffers(const hardware::System& system);
static void update_halo_aos(const std::vector<const hardware::buffers::Gaugemomentum *> buffers, const meta::Inputparameters& params);
static void update_halo_soa(const std::vector<const hardware::buffers::Gaugemomentum *> buffers, const meta::Inputparameters& params);

physics::lattices::Gaugemomenta::Gaugemomenta(const hardware::System& system)
	: system(system), buffers(allocate_buffers(system))
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
		buffer->get_device()->get_gaugemomentum_code()->set_zero_gaugemomentum(buffer);
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
		buf->get_device()->get_gaugemomentum_code()->generate_gaussian_gaugemomenta_device(buf, prng_bufs[i]);
	}
	update_halo();
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
		auto code = device->get_gaugemomentum_code();

		code->set_float_to_gaugemomentum_squarenorm_device(field_buf, res_buf);
	}

	res->sum();
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
			update_halo_soa(buffers, system.get_inputparameters());
		} else {
			update_halo_aos(buffers, system.get_inputparameters());
		}
	}
}

static void update_halo_aos(const std::vector<const hardware::buffers::Gaugemomentum *> buffers, const meta::Inputparameters& params)
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo<ae>(buffers, params, NDIM);
}

static void update_halo_soa(const std::vector<const hardware::buffers::Gaugemomentum *> buffers, const meta::Inputparameters& params)
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(!buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo_soa<ae>(buffers, params, .5, 2 * NDIM);
}
