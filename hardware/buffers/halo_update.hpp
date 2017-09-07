/*
 * Copyright 2013 Matthias Bach
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

#ifndef _HARDWARE_BUFFERS_HALO_UPDATE_
#define _HARDWARE_BUFFERS_HALO_UPDATE_

#include "../size_4.hpp"
#include "plain.hpp"
#include "../device.hpp"
#include <vector>
#include <map>
#include <sstream>
#include <stdexcept>
#include "../system.hpp"
#include "../transfer/transfer.hpp"

unsigned const UP_TRANSFER = 0x1;
unsigned const DOWN_TRANSFER = 0x2;

namespace hardware {

namespace buffers {

template <typename T, class BUFFER> void update_halo(std::vector<BUFFER*> buffers, const hardware::System& system, const float ELEMS_PER_SITE = 1.);
/**
 * Update the halo of the given buffers.
 */
template <typename T, class BUFFER> void update_halo_soa(std::vector<BUFFER*> buffers, const hardware::System& system, const float ELEMS_PER_SITE = 1., const unsigned CHUNKS_PER_LANE = 1, unsigned reqd_width = 0);
/**
 * Initialize the halo update of the given buffers.
 *
 * Extraction of the halo data happens synchroneous to the default queue of each buffer.
 * Therefore commands using the default queue of this buffer can safely operate on the non-halo / own elements of this buffer if queued after this call.
 *
 * Note that depending on the exact transfer method used this might already transfer data between devices!
 */
template <typename T, class BUFFER> void initialize_update_halo_soa(std::vector<BUFFER*> buffers, const hardware::System& system, const float ELEMS_PER_SITE = 1., const unsigned CHUNKS_PER_LANE = 1, unsigned reqd_width = 0);
/**
 * Finish up a previously started halo update.
 *
 * Access the the halo data happens synchroneous to the default queue of each buffer.
 * Therefore commands using the default queue of this buffer can safely operate the halo elements if queued after this call.
 *
 * Note that depending on the exact transfer methods this might cause data transfer between devices.
 */
template <typename T, class BUFFER> void finalize_update_halo_soa(std::vector<BUFFER*> buffers, const hardware::System& system, const float ELEMS_PER_SITE = 1., const unsigned CHUNKS_PER_LANE = 1, unsigned reqd_width = 0);

template<typename BUFFER> static hardware::SynchronizationEvent extract_boundary(hardware::Transfer * transfer, const BUFFER * buffer, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE, const hardware::SynchronizationEvent& event);
template<typename BUFFER> static hardware::SynchronizationEvent send_halo(hardware::Transfer * transfer, const BUFFER * buffer, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE, const hardware::SynchronizationEvent& event);

}

}

//todo: move to cpp!
static unsigned lower_grid_neighbour(const unsigned idx, const unsigned GRID_SIZE);
static unsigned upper_grid_neighbour(const unsigned idx, const unsigned GRID_SIZE);
template <typename T, class BUFFER> void hardware::buffers::update_halo(std::vector<BUFFER*> buffers, const hardware::System& system, const float ELEMS_PER_SITE)
{
	// no-op on single device
	size_t num_buffers = buffers.size();
	if(num_buffers > 1) {
		const auto main_device = buffers[0]->get_device();
		LatticeGrid lG(main_device->getGridSize());
		const unsigned GRID_SIZE = lG.tExtent;
		const unsigned HALO_SIZE = main_device->getHaloExtent();
		const unsigned VOLSPACE = system.getHardwareParameters()->getSpatialLatticeVolume() * ELEMS_PER_SITE;
		const unsigned HALO_ELEMS = HALO_SIZE * VOLSPACE;
		const unsigned VOL4D_LOCAL = get_vol4d(main_device->getLocalLatticeExtents()) * ELEMS_PER_SITE;
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

template<typename BUFFER> static hardware::SynchronizationEvent hardware::buffers::extract_boundary(hardware::Transfer * transfer, const BUFFER * buffer, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE, const hardware::SynchronizationEvent& event)
{
	logger.debug() << "Extracting boundary. Offset: " << in_lane_offset << " - Elements per chunk: " << HALO_CHUNK_ELEMS;
	const unsigned NUM_LANES = buffer->get_lane_count();
	const unsigned STORAGE_TYPE_SIZE = buffer->get_storage_type_size();
	const unsigned CHUNK_STRIDE = get_vol4d(buffer->get_device()->getLocalLatticeMemoryExtents()) * ELEMS_PER_SITE;

	const size_t buffer_origin[] = { in_lane_offset * STORAGE_TYPE_SIZE, 0, 0 };

	const size_t region[] = { HALO_CHUNK_ELEMS * STORAGE_TYPE_SIZE, CHUNKS_PER_LANE, NUM_LANES };

	const size_t buffer_row_pitch = CHUNK_STRIDE * STORAGE_TYPE_SIZE;
	const size_t buffer_slice_pitch = buffer->get_lane_stride() * STORAGE_TYPE_SIZE;

	return transfer->load(buffer, buffer_origin, region, buffer_row_pitch, buffer_slice_pitch, event);

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

template<typename BUFFER> static hardware::SynchronizationEvent hardware::buffers::send_halo(hardware::Transfer* transfer, const BUFFER * buffer, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE, const hardware::SynchronizationEvent& event)
{
	logger.debug() << "Sending Halo. Offset: " << in_lane_offset << " - Elements per chunk: " << HALO_CHUNK_ELEMS;
	const unsigned NUM_LANES = buffer->get_lane_count();
	const unsigned STORAGE_TYPE_SIZE = buffer->get_storage_type_size();
	const unsigned CHUNK_STRIDE = get_vol4d(buffer->get_device()->getLocalLatticeMemoryExtents()) * ELEMS_PER_SITE;

	const size_t buffer_origin[] = { in_lane_offset * STORAGE_TYPE_SIZE, 0, 0 };

	const size_t region[] = { HALO_CHUNK_ELEMS * STORAGE_TYPE_SIZE, CHUNKS_PER_LANE, NUM_LANES };

	const size_t buffer_row_pitch = CHUNK_STRIDE * STORAGE_TYPE_SIZE;
	const size_t buffer_slice_pitch = buffer->get_lane_stride() * STORAGE_TYPE_SIZE;

	return transfer->dump(buffer, buffer_origin, region, buffer_row_pitch, buffer_slice_pitch, event);

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

template<class BUFFER> struct UpdateHaloSOAhelper {
	UpdateHaloSOAhelper(std::vector<BUFFER*> const & buffers, const hardware::System& system, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE, unsigned reqd_width)
	{
		const auto main_device = buffers[0]->get_device();
		LatticeGrid lG(main_device->getGridSize());
		grid_size = lG.tExtent;
		halo_size = main_device->getHaloExtent();
		volspace = system.getHardwareParameters()->getSpatialLatticeVolume() * ELEMS_PER_SITE;
		if(reqd_width > halo_size) {
			std::ostringstream tmp;
			tmp << "Requested halo width is " << reqd_width << ", but maximum halo width is " << halo_size;
			throw std::invalid_argument(tmp.str());
		}
		if(!reqd_width) {
			reqd_width = halo_size;
		}
		total_halo_elems = halo_size * volspace * CHUNKS_PER_LANE;
		total_halo_chunk_elems = halo_size * volspace;
		reqd_halo_elems = reqd_width * volspace * CHUNKS_PER_LANE;
		reqd_halo_chunk_elems = reqd_width * volspace;
		vol4d_local = get_vol4d(main_device->getLocalLatticeExtents()) * ELEMS_PER_SITE;

		logger.trace() << "GRID_SIZE: " << grid_size;
		logger.trace() << "HALO_SIZE: " << halo_size;
		logger.trace() << "reqd_width: " << reqd_width;
		logger.trace() << "eff. VOLSPACE: " << volspace;
		logger.trace() << "TOTAL HALO ELEMS: " << total_halo_elems;
		logger.trace() << "TOTAL HALO ELEMS per CHUNK: " << total_halo_chunk_elems;
		logger.trace() << "REQD HALO ELEMS: " << reqd_halo_elems;
		logger.trace() << "REQD HALO ELEMS per CHUNK: " << reqd_halo_chunk_elems;
		logger.trace() << "eff. VOL4D_LOCAL: " << vol4d_local;
		logger.trace() << "CHUNKS_PER_LANE: " << CHUNKS_PER_LANE;

	}

	unsigned grid_size;
	unsigned halo_size;
	unsigned volspace;
	unsigned total_halo_elems;
	unsigned total_halo_chunk_elems;
	unsigned reqd_halo_elems;
	unsigned reqd_halo_chunk_elems;
	unsigned vol4d_local;
};

template <typename T, class BUFFER> void hardware::buffers::initialize_update_halo_soa(std::vector<BUFFER*> buffers, const hardware::System& system, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE, unsigned reqd_width)
{
	// TODO all transfer imply that buffers are mapped 1:1 to devices, this must be improved!
	const size_t num_buffers = buffers.size();
	if(num_buffers > 1) {
		UpdateHaloSOAhelper<BUFFER> const helper(buffers, system, ELEMS_PER_SITE, CHUNKS_PER_LANE, reqd_width);

		auto const devices = system.get_devices();

		// copy inside of boundaries to host
		for(size_t i = 0; i < num_buffers; ++i) {
			const auto buffer = buffers[i];
			logger.debug() << "Extracting data from buffer " << i;
			auto const up_transfer = system.get_transfer(i, upper_grid_neighbour(i, helper.grid_size), UP_TRANSFER);
			extract_boundary(up_transfer, buffer, helper.vol4d_local - helper.reqd_halo_chunk_elems, helper.reqd_halo_chunk_elems, ELEMS_PER_SITE, CHUNKS_PER_LANE, SynchronizationEvent());
			auto const down_transfer = system.get_transfer(i, lower_grid_neighbour(i, helper.grid_size), DOWN_TRANSFER);
			extract_boundary(down_transfer, buffer, 0, helper.reqd_halo_chunk_elems, ELEMS_PER_SITE, CHUNKS_PER_LANE, SynchronizationEvent());
		}

		for(auto * device: devices) {
			device->flush();
		}

		// trigger transfers (might be a NOOP)
		for(size_t i = 0; i < num_buffers; ++i) {
			auto const up_transfer = system.get_transfer(i, upper_grid_neighbour(i, helper.grid_size), UP_TRANSFER);
			up_transfer->transfer();
			auto const down_transfer = system.get_transfer(i, lower_grid_neighbour(i, helper.grid_size), DOWN_TRANSFER);
			down_transfer->transfer();
		}
	}
}

template <typename T, class BUFFER> void hardware::buffers::finalize_update_halo_soa(std::vector<BUFFER*> buffers, const hardware::System& system, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE, unsigned reqd_width)
{
	const size_t num_buffers = buffers.size();
	if(num_buffers > 1) {
		UpdateHaloSOAhelper<BUFFER> const helper(buffers, system, ELEMS_PER_SITE, CHUNKS_PER_LANE, reqd_width);

		// copy data from host to halo (of course getting what the neighbour stored
		for(size_t i = 0; i < num_buffers; ++i) {
			const auto buffer = buffers[i];
			// our lower halo is the upper bounary of our lower neighbour
			// its storage location is wrapped around to be the last chunk of data in our buffer, that is after local data and upper halo
			auto const up_transfer = system.get_transfer(lower_grid_neighbour(i, helper.grid_size), i, UP_TRANSFER);
			logger.debug() << "Sending data to buffer " << i;
			send_halo(up_transfer, buffer, helper.vol4d_local + 2 * helper.total_halo_chunk_elems - helper.reqd_halo_chunk_elems, helper.reqd_halo_chunk_elems, ELEMS_PER_SITE, CHUNKS_PER_LANE, SynchronizationEvent());
			// our upper halo is the lower bounary of our upper neighbour, it's stored right after our local data
			auto const down_transfer = system.get_transfer(upper_grid_neighbour(i, helper.grid_size), i, DOWN_TRANSFER);
			send_halo(down_transfer, buffer, helper.vol4d_local, helper.reqd_halo_chunk_elems, ELEMS_PER_SITE, CHUNKS_PER_LANE, SynchronizationEvent());
		}
	}
}

template <typename T, class BUFFER> void hardware::buffers::update_halo_soa(std::vector<BUFFER*> buffers, const hardware::System& system, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE, unsigned reqd_width)
{
	initialize_update_halo_soa<T>(buffers, system, ELEMS_PER_SITE, CHUNKS_PER_LANE, reqd_width);
	finalize_update_halo_soa<T>(buffers, system, ELEMS_PER_SITE, CHUNKS_PER_LANE, reqd_width);
}

#endif /* _HARDWARE_BUFFERS_HALO_UPDATE_ */
