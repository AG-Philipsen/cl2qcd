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

static std::vector<const hardware::buffers::Gaugemomentum *> allocate_buffers(const hardware::System& system);

physics::lattices::Gaugemomenta::Gaugemomenta(const hardware::System& system)
	: system(system), buffers(allocate_buffers(system))
{
}

static  std::vector<const hardware::buffers::Gaugemomentum *> allocate_buffers(const hardware::System& system)
{
	using hardware::buffers::Gaugemomentum;

	// only use device 0 for now
	hardware::Device * device = system.get_devices().at(0);
	std::vector<const Gaugemomentum*> buffers;
	buffers.push_back(new Gaugemomentum(NDIM * get_vol4d(device->get_mem_lattice_size()), device));
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

	// TODO implemente for more than one device
	if(num_buffers != 1) {
		throw Print_Error_Message("physics::lattices::squarenorm(const Gaugemomenta&) is not implemented for multiple devices", __FILE__, __LINE__);
	}
	if(num_buffers != res_buffers.size()) {
		throw std::invalid_argument("The given lattices do not sue the same number of devices.");
	}

	auto field_buf = field_buffers[0];
	auto res_buf = res_buffers[0];
	auto device = field_buf->get_device();
	auto code = device->get_gaugemomentum_code();

	code->set_float_to_gaugemomentum_squarenorm_device(field_buf, res_buf);
}

void physics::lattices::log_squarenorm(const std::string& msg, const physics::lattices::Gaugemomenta& x)
{
	if(logger.beDebug()) {
		hmc_float tmp = squarenorm(x);
		logger.debug() << msg << std::scientific << std::setprecision(10) << tmp;
	}
}
