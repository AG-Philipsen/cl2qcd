/** @file
 * Implementation of the physics::lattices::Spinorfield class
 */

#include "spinorfield.hpp"
#include "../../meta/util.hpp"
#include <cassert>

static std::vector<const hardware::buffers::Plain<spinor> *> allocate_buffers(const hardware::System& system);

physics::lattices::Spinorfield::Spinorfield(const hardware::System& system)
	: system(system), buffers(allocate_buffers(system))
{
}

static  std::vector<const hardware::buffers::Plain<spinor> *> allocate_buffers(const hardware::System& system)
{
	using hardware::buffers::Plain;

	// only use device 0 for now
	hardware::Device * device = system.get_devices().at(0);
	std::vector<const Plain<spinor>*> buffers;
	buffers.push_back(new Plain<spinor>(meta::get_spinorfieldsize(system.get_inputparameters()), device));
	return buffers;
}

physics::lattices::Spinorfield::~Spinorfield()
{
for(auto buffer: buffers) {
		delete buffer;
	}
}

const std::vector<const physics::lattices::Spinorfield *> physics::lattices::create_spinorfields(const hardware::System& system, const size_t n)
{
	std::vector<const Spinorfield *> fields;
	fields.reserve(n);

	for(size_t i = 0; i < n; ++i) {
		fields.push_back(new Spinorfield(system));
	}

	return fields;
}

void physics::lattices::release_spinorfields(const std::vector<const physics::lattices::Spinorfield *> fields)
{
for(auto field: fields) {
		delete field;
	}
}

const std::vector<const hardware::buffers::Plain<spinor> *> physics::lattices::Spinorfield::get_buffers() const noexcept
{
	return buffers;
}

void physics::lattices::Spinorfield::gamma5() const
{
for(auto buffer: buffers) {
		auto fermion_code = buffer->get_device()->get_fermion_code();
		fermion_code->gamma5_device(buffer);
	}
}

hmc_complex physics::lattices::scalar_product(const Spinorfield& left, const Spinorfield& right)
{
	auto left_buffers = left.get_buffers();
	auto right_buffers = right.get_buffers();

	// TODO implemente for more than one device
	if(left_buffers.size() > 1 || right_buffers.size() > 1) {
		throw Print_Error_Message("physics::lattices::scalar_product(const Spinorfield&, const Spinorfield&) is not implemented for multiple devices", __FILE__, __LINE__);
	}

	auto left_buf = left_buffers[0];
	auto right_buf = right_buffers[0];
	auto device = left_buf->get_device();
	hardware::buffers::Plain<hmc_complex> result_buf(1, device);
	auto spinor_code = device->get_spinor_code();

	spinor_code->set_complex_to_scalar_product_device(left_buf, right_buf, &result_buf);

	hmc_complex result;
	result_buf.dump(&result);
	return result;
}
