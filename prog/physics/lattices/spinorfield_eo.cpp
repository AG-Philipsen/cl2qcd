/** @file
 * Implementation of the physics::lattices::Spinorfield_eo class
 */

#include "spinorfield_eo.hpp"
#include "../../meta/util.hpp"
#include <cassert>
#include "../../hardware/code/spinors.hpp"
#include "../../hardware/code/fermions.hpp"
#include "../../meta/type_ops.hpp"
#include "../../hardware/buffers/halo_update.hpp"

static std::vector<const hardware::buffers::Spinor *> allocate_buffers(const hardware::System& system);
static void update_halo_soa(const std::vector<const hardware::buffers::Spinor *> buffers, const meta::Inputparameters& params);
static void update_halo_aos(const std::vector<const hardware::buffers::Spinor *> buffers, const meta::Inputparameters& params);
static void extract_boundary(char* host, const hardware::buffers::Spinor * buffer, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS);
static void send_halo(const hardware::buffers::Spinor * buffer, const char* host, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS);

physics::lattices::Spinorfield_eo::Spinorfield_eo(const hardware::System& system)
	: system(system), buffers(allocate_buffers(system))
{
}

static  std::vector<const hardware::buffers::Spinor *> allocate_buffers(const hardware::System& system)
{
	using hardware::buffers::Spinor;

	auto devices = system.get_devices();
	std::vector<const Spinor*> buffers;
	buffers.reserve(devices.size());
	for(auto device: devices) {
		buffers.push_back(new Spinor(hardware::code::get_eoprec_spinorfieldsize(device->get_mem_lattice_size()), device));
	}
	return buffers;
}

physics::lattices::Spinorfield_eo::~Spinorfield_eo()
{
for(auto buffer: buffers) {
		delete buffer;
	}
}

const std::vector<const hardware::buffers::Spinor *> physics::lattices::Spinorfield_eo::get_buffers() const noexcept
{
	return buffers;
}

hmc_complex physics::lattices::scalar_product(const Spinorfield_eo& left, const Spinorfield_eo& right)
{
	const Scalar<hmc_complex> res(left.system);
	scalar_product(&res, left, right);
	return res.get();
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
		auto spinor_code = device->get_spinor_code();

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
		auto spinor_code = device->get_spinor_code();

		spinor_code->set_float_to_global_squarenorm_eoprec_device(field_buf, res_buf);
	}
	res->sum();
}

void physics::lattices::Spinorfield_eo::zero() const
{
for(auto buffer: buffers) {
		auto spinor_code = buffer->get_device()->get_spinor_code();
		spinor_code->set_zero_spinorfield_eoprec_device(buffer);
	}
}

void physics::lattices::Spinorfield_eo::cold() const
{
for(auto buffer: buffers) {
		auto spinor_code = buffer->get_device()->get_spinor_code();
		spinor_code->set_eoprec_spinorfield_cold_device(buffer);
	}
}

void physics::lattices::Spinorfield_eo::gaussian(const physics::PRNG& prng) const
{
	auto prng_bufs = prng.get_buffers();

	if(buffers.size() != prng_bufs.size()) {
		throw std::invalid_argument("PRNG does not use same devices as spinorfield");
	}

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto spin_buf = buffers[i];
		auto prng_buf = prng_bufs[i];
		spin_buf->get_device()->get_spinor_code()->generate_gaussian_spinorfield_eo_device(spin_buf, prng_buf);
	}
	update_halo();
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
		device->get_spinor_code()->saxpy_eoprec_device(x_bufs[i], y_bufs[i], alpha, out_buf);
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
		device->get_spinor_code()->saxpy_eoprec_device(x_bufs[i], y_bufs[i], alpha_bufs[i], out_buf);
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
		device->get_spinor_code()->sax_eoprec_device(x_bufs[i], alpha_bufs[i], out_buf);
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
		device->get_spinor_code()->saxsbypz_eoprec_device(x_bufs[i], y_bufs[i], z_bufs[i], alpha_bufs[i], beta_bufs[i], out_buf);
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

void physics::lattices::Spinorfield_eo::gamma5() const
{
for(auto buffer: buffers) {
		auto fermion_code = buffer->get_device()->get_fermion_code();
		fermion_code->gamma5_eo_device(buffer);
	}
}

template<> size_t physics::lattices::get_flops<physics::lattices::Spinorfield_eo, physics::lattices::scalar_product>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	if(devices.size() != 1) {
		throw Print_Error_Message("Currently only supported on a single device", __FILE__, __LINE__);
	}

	auto spinor_code = devices[0]->get_spinor_code();
	return spinor_code->get_flop_size("scalar_product_eoprec");
}

template<> size_t physics::lattices::get_flops<physics::lattices::Spinorfield_eo, physics::lattices::squarenorm>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	if(devices.size() != 1) {
		throw Print_Error_Message("Currently only supported on a single device", __FILE__, __LINE__);
	}

	auto spinor_code = devices[0]->get_spinor_code();
	return spinor_code->get_flop_size("global_squarenorm_eoprec");
}

template<> size_t physics::lattices::get_flops<physics::lattices::Spinorfield_eo, physics::lattices::sax>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	if(devices.size() != 1) {
		throw Print_Error_Message("Currently only supported on a single device", __FILE__, __LINE__);
	}

	auto spinor_code = devices[0]->get_spinor_code();
	return spinor_code->get_flop_size("sax_eoprec");
}
template<> size_t physics::lattices::get_flops<physics::lattices::Spinorfield_eo, physics::lattices::saxpy>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	if(devices.size() != 1) {
		throw Print_Error_Message("Currently only supported on a single device", __FILE__, __LINE__);
	}

	auto spinor_code = devices[0]->get_spinor_code();
	return spinor_code->get_flop_size("saxpy_eoprec");
}
template<> size_t physics::lattices::get_flops<physics::lattices::Spinorfield_eo, physics::lattices::saxsbypz>(const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	if(devices.size() != 1) {
		throw Print_Error_Message("Currently only supported on a single device", __FILE__, __LINE__);
	}

	auto spinor_code = devices[0]->get_spinor_code();
	return spinor_code->get_flop_size("saxsbypz_eoprec");
}

void physics::lattices::log_squarenorm(const std::string& msg, const physics::lattices::Spinorfield_eo& x)
{
	if(logger.beDebug()) {
		hmc_float tmp = squarenorm(x);
		logger.debug() << msg << std::scientific << std::setprecision(10) << tmp;
	}
}

void physics::lattices::Spinorfield_eo::update_halo() const
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

	const auto main_device = buffers[0]->get_device();
	const size_4 grid_dims = main_device->get_grid_size();
	if(grid_dims.x != 1 || grid_dims.y != 1 || grid_dims.z != 1) {
		throw Print_Error_Message("Only the time-direction can be parallelized");
	}
	const unsigned GRID_SIZE = grid_dims.t;
	const unsigned HALO_SIZE = main_device->get_halo_size();
	const unsigned VOLSPACE = meta::get_volspace(params) / 2;
	const unsigned HALO_ELEMS = HALO_SIZE * VOLSPACE;
	const unsigned HALO_CHUNK_ELEMS = HALO_SIZE * VOLSPACE;
	const unsigned VOL4D_LOCAL = get_vol4d(main_device->get_local_lattice_size()) / 2;
	const size_t num_buffers = buffers.size();

	std::vector<char*> upper_boundaries; upper_boundaries.reserve(num_buffers);
	std::vector<char*> lower_boundaries; lower_boundaries.reserve(num_buffers);
	for(size_t i = 0; i < num_buffers; ++i) {
		upper_boundaries.push_back(new char[HALO_ELEMS * sizeof(spinor)]);
		lower_boundaries.push_back(new char[HALO_ELEMS * sizeof(spinor)]);
	}

	// copy inside of boundaries to host
	for(size_t i = 0; i < num_buffers; ++i) {
		const auto buffer = buffers[i];
		logger.debug() << "Extracting data from buffer " << i;
		extract_boundary(upper_boundaries[i], buffer, VOL4D_LOCAL - HALO_CHUNK_ELEMS, HALO_CHUNK_ELEMS);
		extract_boundary(lower_boundaries[i], buffer, 0, HALO_CHUNK_ELEMS);
	}

	// copy data from host to halo (of coure getting what the neighbour stored
	for(size_t i = 0; i < num_buffers; ++i) {
		const auto buffer = buffers[i];
		logger.debug() << "Sending data to buffer " << i;
		// our lower halo is the upper bounary of our lower neighbour
		// its storage location is wrapped around to be the last chunk of data in our buffer, that is after local data and upper halo
		send_halo(buffer, upper_boundaries[lower_grid_neighbour(i, GRID_SIZE)], VOL4D_LOCAL  + HALO_CHUNK_ELEMS, HALO_CHUNK_ELEMS);
		// our upper halo is the lower bounary of our upper neighbour, it's stored right after our local data
		send_halo(buffer, lower_boundaries[upper_grid_neighbour(i, GRID_SIZE)], VOL4D_LOCAL, HALO_CHUNK_ELEMS);
	}

	// clean up host
	for(size_t i = 0; i < num_buffers; ++i) {
		delete[] upper_boundaries[i];
		delete[] lower_boundaries[i];
	}
}

static void extract_boundary(char* host, const hardware::buffers::Spinor * buffer, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS)
{
	logger.debug() << "Extracting boundary. Offset: " << in_lane_offset << " - Elements per chunk: " << HALO_CHUNK_ELEMS;
	const unsigned NUM_LANES = buffer->get_lane_count();
	const unsigned STORAGE_TYPE_SIZE = buffer->get_storage_type_size();
	const unsigned VOL4D_MEM = get_vol4d(buffer->get_device()->get_mem_lattice_size()) / 2;
	const unsigned CHUNKS_PER_LANE = 1;
	const unsigned CHUNK_STRIDE = VOL4D_MEM;

	for(size_t lane = 0; lane < NUM_LANES; ++lane) {
		logger.trace() << "Reading lane " << lane;
		size_t lane_offset = lane * buffer->get_lane_stride();
		for(size_t chunk = 0; chunk < CHUNKS_PER_LANE; ++chunk) {
			size_t host_offset = (lane * CHUNKS_PER_LANE + chunk) * HALO_CHUNK_ELEMS;
			size_t dev_offset = lane_offset + in_lane_offset + chunk * CHUNK_STRIDE;
			buffer->dump_raw(&host[host_offset * STORAGE_TYPE_SIZE], HALO_CHUNK_ELEMS * STORAGE_TYPE_SIZE, dev_offset * STORAGE_TYPE_SIZE);
		}
	}
}

static void send_halo(const hardware::buffers::Spinor * buffer, const char* host, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS)
{
	logger.debug() << "Sending Halo. Offset: " << in_lane_offset << " - Elements per chunk: " << HALO_CHUNK_ELEMS;
	const unsigned NUM_LANES = buffer->get_lane_count();
	const unsigned STORAGE_TYPE_SIZE = buffer->get_storage_type_size();
	const unsigned VOL4D_MEM = get_vol4d(buffer->get_device()->get_mem_lattice_size()) / 2;
	const unsigned CHUNKS_PER_LANE = 1;
	const unsigned CHUNK_STRIDE = VOL4D_MEM;

	for(size_t lane = 0; lane < NUM_LANES; ++lane) {
		logger.trace() << "Sending lane " << lane;
		size_t lane_offset = lane * buffer->get_lane_stride();
		for(size_t chunk = 0; chunk < CHUNKS_PER_LANE; ++chunk) {
			size_t host_offset = (lane * CHUNKS_PER_LANE + chunk) * HALO_CHUNK_ELEMS;
			size_t dev_offset = lane_offset + in_lane_offset + chunk * CHUNK_STRIDE;
			buffer->load_raw(&host[host_offset * STORAGE_TYPE_SIZE], HALO_CHUNK_ELEMS * STORAGE_TYPE_SIZE, dev_offset * STORAGE_TYPE_SIZE);
		}
	}
}
