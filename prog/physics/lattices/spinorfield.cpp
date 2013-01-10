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

hmc_float physics::lattices::squarenorm(const Spinorfield& field)
{
	auto field_buffers = field.get_buffers();

	// TODO implemente for more than one device
	if(field_buffers.size() > 1) {
		throw Print_Error_Message("physics::lattices::squarenorm(const Spinorfield&) is not implemented for multiple devices", __FILE__, __LINE__);
	}

	auto field_buf = field_buffers[0];
	auto device = field_buf->get_device();
	hardware::buffers::Plain<hmc_float> result_buf(1, device);
	auto spinor_code = device->get_spinor_code();

	spinor_code->set_float_to_global_squarenorm_device(field_buf, &result_buf);

	hmc_float result;
	result_buf.dump(&result);
	return result;
}

void physics::lattices::Spinorfield::zero() const
{
for(auto buffer: buffers) {
		auto spinor_code = buffer->get_device()->get_spinor_code();
		spinor_code->set_zero_spinorfield_device(buffer);
	}
}

void physics::lattices::Spinorfield::cold() const
{
for(auto buffer: buffers) {
		auto spinor_code = buffer->get_device()->get_spinor_code();
		spinor_code->set_spinorfield_cold_device(buffer);
	}
}

void physics::lattices::Spinorfield::gaussian(const physics::PRNG& prng) const
{
	auto prng_bufs = prng.get_buffers();

	if(buffers.size() != prng_bufs.size()) {
		throw std::invalid_argument("PRNG does not use same devices as spinorfield");
	}

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto spin_buf = buffers[i];
		auto prng_buf = prng_bufs[i];
		spin_buf->get_device()->get_spinor_code()->generate_gaussian_spinorfield_device(spin_buf, prng_buf);
	}
}

void physics::lattices::saxpy(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x, const Spinorfield& y)
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
		hardware::buffers::Plain<hmc_complex> alpha_buf(1, device);
		alpha_buf.load(&alpha);
		device->get_spinor_code()->saxpy_device(x_bufs[i], y_bufs[i], &alpha_buf, out_buf);
	}
}

void physics::lattices::sax(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x)
{
	auto out_bufs = out->get_buffers();
	auto x_bufs = x.get_buffers();

	if(out_bufs.size() != x_bufs.size()) {
		throw std::invalid_argument("Output buffers does not use same devices as input buffers");
	}

	for(size_t i = 0; i < out_bufs.size(); ++i) {
		auto out_buf = out_bufs[i];
		auto device = out_buf->get_device();
		hardware::buffers::Plain<hmc_complex> alpha_buf(1, device);
		alpha_buf.load(&alpha);
		device->get_spinor_code()->sax_device(x_bufs[i], &alpha_buf, out_buf);
	}
}

void physics::lattices::saxsbypz(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x, const hmc_complex beta, const Spinorfield& y, const Spinorfield& z)
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
		hardware::buffers::Plain<hmc_complex> alpha_buf(1, device);
		hardware::buffers::Plain<hmc_complex> beta_buf(1, device);
		alpha_buf.load(&alpha);
		beta_buf.load(&beta);
		device->get_spinor_code()->saxsbypz_device(x_bufs[i], y_bufs[i], z_bufs[i], &alpha_buf, &beta_buf, out_buf);
	}
}
