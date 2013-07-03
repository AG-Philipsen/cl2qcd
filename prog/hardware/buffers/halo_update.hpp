#ifndef _HARDWARE_BUFFERS_HALO_UPDATE_
#define _HARDWARE_BUFFERS_HALO_UPDATE_

#include "../../meta/inputparameters.hpp"
#include "../../meta/size_4.hpp"
#include "../../meta/util.hpp"
#include "plain.hpp"
#include "../../exceptions.h"
#include "../device.hpp"
#include <vector>
#include <map>
#include <sstream>
#include <stdexcept>
#include "../system.hpp"

namespace hardware {

namespace buffers {

template <typename T, class BUFFER> void update_halo(std::vector<BUFFER*> buffers, const meta::Inputparameters& params, const float ELEMS_PER_SITE = 1.);
template <typename T, class BUFFER> void update_halo_soa(std::vector<BUFFER*> buffers, const hardware::System& system, const float ELEMS_PER_SITE = 1., const unsigned CHUNKS_PER_LANE = 1, unsigned reqd_width = 0);

template<typename BUFFER, typename PROXY> static hardware::SynchronizationEvent extract_boundary(const PROXY* proxy, const BUFFER * buffer, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE, const hardware::SynchronizationEvent& event);
template<typename BUFFER, typename PROXY> static hardware::SynchronizationEvent send_halo(const BUFFER * buffer, const PROXY* proxy, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE, const hardware::SynchronizationEvent& event);

}

}


static unsigned lower_grid_neighbour(const unsigned idx, const unsigned GRID_SIZE);
static unsigned upper_grid_neighbour(const unsigned idx, const unsigned GRID_SIZE);
template <typename T, class BUFFER> void hardware::buffers::update_halo(std::vector<BUFFER*> buffers, const meta::Inputparameters& params, const float ELEMS_PER_SITE)
{
	// no-op on single device
	size_t num_buffers = buffers.size();
	if(num_buffers > 1) {
		const auto main_device = buffers[0]->get_device();
		const size_4 grid_dims = main_device->get_grid_size();
		if(grid_dims.x != 1 || grid_dims.y != 1 || grid_dims.z != 1) {
			throw Print_Error_Message("Only the time-direction can be parallelized");
		}
		const unsigned GRID_SIZE = grid_dims.t;
		const unsigned HALO_SIZE = main_device->get_halo_size();
		const unsigned VOLSPACE = meta::get_volspace(params) * ELEMS_PER_SITE;
		const unsigned HALO_ELEMS = HALO_SIZE * VOLSPACE;
		const unsigned VOL4D_LOCAL = get_vol4d(main_device->get_local_lattice_size()) * ELEMS_PER_SITE;
		const size_t num_buffers = buffers.size();

		// host buffers for intermediate data storage (no direct device to device copy)
		std::vector<T*> upper_boundaries;
		upper_boundaries.reserve(num_buffers);
		std::vector<T*> lower_boundaries;
		lower_boundaries.reserve(num_buffers);
		for(size_t i = 0; i < num_buffers; ++i) {
			upper_boundaries.push_back(new T[HALO_ELEMS]);
			lower_boundaries.push_back(new T[HALO_ELEMS]);
		}

		// copy inside of boundaries to host
		for(size_t i = 0; i < num_buffers; ++i) {
			const auto buffer = buffers[i];
			buffer->dump(upper_boundaries[i], HALO_ELEMS, VOL4D_LOCAL - HALO_ELEMS);
			buffer->dump(lower_boundaries[i], HALO_ELEMS, 0);
		}

		// copy data from host to halo (of coure getting what the neighbour stored
		for(size_t i = 0; i < num_buffers; ++i) {
			const auto buffer = buffers[i];
			// our lower halo is the upper bounary of our lower neighbour
			// its storage location is wrapped around to be the last chunk of data in our buffer, that is after local data and upper halo
			buffer->load(upper_boundaries[lower_grid_neighbour(i, GRID_SIZE)], HALO_ELEMS, VOL4D_LOCAL + HALO_ELEMS);
			// our upper halo is the lower bounary of our upper neighbour, it's stored right after our local data
			buffer->load(lower_boundaries[upper_grid_neighbour(i, GRID_SIZE)], HALO_ELEMS, VOL4D_LOCAL);
		}

		// clean up host
		for(size_t i = 0; i < num_buffers; ++i) {
			delete[] upper_boundaries[i];
			delete[] lower_boundaries[i];
		}
	}
}

static unsigned upper_grid_neighbour(const unsigned idx, const unsigned GRID_SIZE)
{
	return (idx + 1) % GRID_SIZE;
}

static unsigned lower_grid_neighbour(const unsigned idx, const unsigned GRID_SIZE)
{
	return (idx + GRID_SIZE - 1) % GRID_SIZE;
}

template<typename BUFFER, typename PROXY> static hardware::SynchronizationEvent hardware::buffers::extract_boundary(const PROXY* proxy, const BUFFER * buffer, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE, const hardware::SynchronizationEvent& event)
{
	logger.debug() << "Extracting boundary. Offset: " << in_lane_offset << " - Elements per chunk: " << HALO_CHUNK_ELEMS;
	const unsigned NUM_LANES = buffer->get_lane_count();
	const unsigned STORAGE_TYPE_SIZE = buffer->get_storage_type_size();
	const unsigned CHUNK_STRIDE = get_vol4d(buffer->get_device()->get_mem_lattice_size()) * ELEMS_PER_SITE;

	const size_t buffer_origin[] = { in_lane_offset * STORAGE_TYPE_SIZE, 0, 0 };
	const size_t proxy_origin[] = { 0, 0, 0 };

	const size_t region[] = { HALO_CHUNK_ELEMS * STORAGE_TYPE_SIZE, CHUNKS_PER_LANE, NUM_LANES };

	const size_t buffer_row_pitch = CHUNK_STRIDE * STORAGE_TYPE_SIZE;
	const size_t proxy_row_pitch = 0; // automatically consecutive
	const size_t buffer_slice_pitch = buffer->get_lane_stride() * STORAGE_TYPE_SIZE;
	const size_t proxy_slice_pitch = 0; // automatically consecutive

	return copyDataRect(buffer->get_device(), proxy, buffer, proxy_origin, buffer_origin, region, proxy_row_pitch, proxy_slice_pitch, buffer_row_pitch, buffer_slice_pitch, {event});

// Original code kept for documentation purposes:
//	for(size_t lane = 0; lane < NUM_LANES; ++lane) {
//		logger.trace() << "Reading lane " << lane;
//		size_t lane_offset = lane * buffer->get_lane_stride();
//		for(size_t chunk = 0; chunk < CHUNKS_PER_LANE; ++chunk) {
//			size_t host_offset = (lane * CHUNKS_PER_LANE + chunk) * HALO_CHUNK_ELEMS;
//			size_t dev_offset = lane_offset + in_lane_offset + chunk * CHUNK_STRIDE;
//			logger.trace() << "Chunk " << chunk << " - host offset: " << host_offset << " - device offset: " << dev_offset;
//			buffer->dump_raw(&host[host_offset * STORAGE_TYPE_SIZE], HALO_CHUNK_ELEMS * STORAGE_TYPE_SIZE, dev_offset * STORAGE_TYPE_SIZE);
//		}
//	}
}

template<typename BUFFER, typename PROXY> static hardware::SynchronizationEvent hardware::buffers::send_halo(const BUFFER * buffer, const PROXY* proxy, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE, const hardware::SynchronizationEvent& event)
{
	logger.debug() << "Sending Halo. Offset: " << in_lane_offset << " - Elements per chunk: " << HALO_CHUNK_ELEMS;
	const unsigned NUM_LANES = buffer->get_lane_count();
	const unsigned STORAGE_TYPE_SIZE = buffer->get_storage_type_size();
	const unsigned CHUNK_STRIDE = get_vol4d(buffer->get_device()->get_mem_lattice_size()) * ELEMS_PER_SITE;

	const size_t buffer_origin[] = { in_lane_offset * STORAGE_TYPE_SIZE, 0, 0 };
	const size_t proxy_origin[] = { 0, 0, 0 };

	const size_t region[] = { HALO_CHUNK_ELEMS * STORAGE_TYPE_SIZE, CHUNKS_PER_LANE, NUM_LANES };

	const size_t buffer_row_pitch = CHUNK_STRIDE * STORAGE_TYPE_SIZE;
	const size_t proxy_row_pitch = 0; // automatically consecutive
	const size_t buffer_slice_pitch = buffer->get_lane_stride() * STORAGE_TYPE_SIZE;
	const size_t proxy_slice_pitch = 0; // automatically consecutive

	return copyDataRect(buffer->get_device(), buffer, proxy, buffer_origin, proxy_origin, region, buffer_row_pitch, buffer_slice_pitch, proxy_row_pitch, proxy_slice_pitch, {event});

// Original code kept for documentation purposes:
//	for(size_t lane = 0; lane < NUM_LANES; ++lane) {
//		logger.trace() << "Sending lane " << lane;
//		size_t lane_offset = lane * buffer->get_lane_stride();
//		for(size_t chunk = 0; chunk < CHUNKS_PER_LANE; ++chunk) {
//			size_t host_offset = (lane * CHUNKS_PER_LANE + chunk) * HALO_CHUNK_ELEMS;
//			size_t dev_offset = lane_offset + in_lane_offset + chunk * CHUNK_STRIDE;
//			logger.trace() << "Chunk " << chunk << " - host offset: " << host_offset << " - device offset: " << dev_offset;
//			buffer->load_raw(&host[host_offset * STORAGE_TYPE_SIZE], HALO_CHUNK_ELEMS * STORAGE_TYPE_SIZE, dev_offset * STORAGE_TYPE_SIZE);
//		}
//	}
}

namespace hardware {
namespace buffers {
class DeviceAccessibleMemory : public Buffer {
	public:
		DeviceAccessibleMemory(const size_t bytes, hardware::Device * device);
};

class ProxyBufferCache {
	private:
		ProxyBufferCache();
		~ProxyBufferCache();
		std::map<std::pair<cl_context,std::pair<size_t,size_t> >, std::vector<DeviceAccessibleMemory*>> cache;
	public:
		static ProxyBufferCache& getInstance();
		const std::vector<DeviceAccessibleMemory*>& getBuffers(size_t rows, size_t bytes, const std::vector<hardware::Device*>& devices);
};

}
}

template <typename T, class BUFFER> void hardware::buffers::update_halo_soa(std::vector<BUFFER*> buffers, const hardware::System& system, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE, unsigned reqd_width = 0)
{
	auto const & params = system.get_inputparameters();
	const size_t num_buffers = buffers.size();

	if(num_buffers > 1) {
		const auto main_device = buffers[0]->get_device();
		const size_4 grid_dims = main_device->get_grid_size();
		if(grid_dims.x != 1 || grid_dims.y != 1 || grid_dims.z != 1) {
			throw Print_Error_Message("Only the time-direction can be parallelized");
		}
		const unsigned GRID_SIZE = grid_dims.t;
		const unsigned HALO_SIZE = main_device->get_halo_size();
		const unsigned VOLSPACE = meta::get_volspace(params) * ELEMS_PER_SITE;
		if(reqd_width > HALO_SIZE) {
			std::ostringstream tmp;
			tmp << "Requested halo width is " << reqd_width << ", but maximum halo width is " << HALO_SIZE;
			throw std::invalid_argument(tmp.str());
		}
		if(!reqd_width) {
			reqd_width = HALO_SIZE;
		}
		const unsigned TOTAL_HALO_ELEMS = HALO_SIZE * VOLSPACE * CHUNKS_PER_LANE;
		const unsigned TOTAL_HALO_CHUNK_ELEMS = HALO_SIZE * VOLSPACE;
		const unsigned REQD_HALO_ELEMS = reqd_width * VOLSPACE * CHUNKS_PER_LANE;
		const unsigned REQD_HALO_CHUNK_ELEMS = reqd_width * VOLSPACE;
		const unsigned VOL4D_LOCAL = get_vol4d(main_device->get_local_lattice_size()) * ELEMS_PER_SITE;

		logger.trace() << "GRID_SIZE: " << GRID_SIZE;
		logger.trace() << "HALO_SIZE: " << HALO_SIZE;
		logger.trace() << "reqd_width: " << reqd_width;
		logger.trace() << "eff. VOLSPACE: " << VOLSPACE;
		logger.trace() << "TOTAL HALO ELEMS: " << TOTAL_HALO_ELEMS;
		logger.trace() << "TOTAL HALO ELEMS per CHUNK: " << TOTAL_HALO_CHUNK_ELEMS;
		logger.trace() << "REQD HALO ELEMS: " << REQD_HALO_ELEMS;
		logger.trace() << "REQD HALO ELEMS per CHUNK: " << REQD_HALO_CHUNK_ELEMS;
		logger.trace() << "eff. VOL4D_LOCAL: " << VOL4D_LOCAL;
		logger.trace() << "CHUNKS_PER_LANE: " << CHUNKS_PER_LANE;

		std::vector<hardware::Device*> devices(num_buffers);
		for(size_t i = 0; i < num_buffers; ++i) {
			devices[i] = buffers[i]->get_device();
		}

		auto& cache = ProxyBufferCache::getInstance();
		auto& host_buffers = cache.getBuffers(2, REQD_HALO_ELEMS * sizeof(T), devices);

		std::vector<hardware::SynchronizationEvent> extract_events(2 * num_buffers);

		// copy inside of boundaries to host
		for(size_t i = 0; i < num_buffers; ++i) {
			const auto buffer = buffers[i];
			logger.debug() << "Extracting data from buffer " << i;
			extract_events[i] = extract_boundary(host_buffers[i], buffer, VOL4D_LOCAL - REQD_HALO_CHUNK_ELEMS, REQD_HALO_CHUNK_ELEMS, ELEMS_PER_SITE, CHUNKS_PER_LANE, SynchronizationEvent());
			extract_events[num_buffers + i] = extract_boundary(host_buffers[num_buffers + i], buffer, 0, REQD_HALO_CHUNK_ELEMS, ELEMS_PER_SITE, CHUNKS_PER_LANE, SynchronizationEvent());
		}

		for(auto device: devices) {
			device->flush();
		}

		std::vector<hardware::SynchronizationEvent> send_events(2 * num_buffers);

		// copy data from host to halo (of coure getting what the neighbour stored
		for(size_t i = 0; i < num_buffers; ++i) {
			const auto buffer = buffers[i];
			logger.debug() << "Sending data to buffer " << i;
			// our lower halo is the upper bounary of our lower neighbour
			// its storage location is wrapped around to be the last chunk of data in our buffer, that is after local data and upper halo
			const size_t lower_i = lower_grid_neighbour(i, GRID_SIZE);
			send_events[lower_i] = send_halo(buffer, host_buffers[lower_i], VOL4D_LOCAL + 2 * TOTAL_HALO_CHUNK_ELEMS - REQD_HALO_CHUNK_ELEMS, REQD_HALO_CHUNK_ELEMS, ELEMS_PER_SITE, CHUNKS_PER_LANE, extract_events[lower_i]);
			// our upper halo is the lower bounary of our upper neighbour, it's stored right after our local data
			const size_t upper_i = upper_grid_neighbour(i, GRID_SIZE);
			send_events[num_buffers + upper_i] = send_halo(buffer, host_buffers[num_buffers + upper_i], VOL4D_LOCAL, REQD_HALO_CHUNK_ELEMS, ELEMS_PER_SITE, CHUNKS_PER_LANE, extract_events[num_buffers + upper_i]);
		}

		// ensure that command queue are blocked until corresponding proxy buffer has been completely read by neighbouring device
#ifdef CL_VERSION_1_2
		// and send buffer back to the correct device while we are at it, dumping their content as we we don't need it anymore, anyways
		std::vector<hardware::SynchronizationEvent> events(2);
		for(size_t i = 0; i < num_buffers; ++i) {
			events[0] = send_events[i];
			events[1] = send_events[num_buffers + i];
			host_buffers[i]->migrate(devices[i], events, CL_MIGRATE_MEM_OBJECT_CONTENT_UNDEFINED);
		}
#else
		for(size_t i = 0; i < num_buffers; ++i) {
			auto device = devices[i];
			device->enqueue_barrier(send_events[i], send_events[num_buffers + i]);
		}
#endif
	}
}

#endif /* _HARDWARE_BUFFERS_HALO_UPDATE_ */
