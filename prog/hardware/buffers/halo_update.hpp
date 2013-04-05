#ifndef _HARDWARE_BUFFERS_PLAIN_
#define _HARDWARE_BUFFERS_PLAIN_

#include "plain.hpp"

namespace hardware {

namespace buffers {

template <typename T, class BUFFER> void update_halo(std::vector<BUFFER*> buffers, const meta::Inputparameters& params, const float ELEMS_PER_SITE = 1.);
template <typename T, class BUFFER> void update_halo_soa(std::vector<BUFFER*> buffers, const meta::Inputparameters& params, const float ELEMS_PER_SITE = 1., const unsigned CHUNKS_PER_LANE = 1);

}

}


static unsigned lower_grid_neighbour(const unsigned idx, const unsigned GRID_SIZE);
static unsigned upper_grid_neighbour(const unsigned idx, const unsigned GRID_SIZE);
template<class BUFFER> static void send_halo(const BUFFER * buffer, const char* host, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS, const float ELEMS_PER_SITE = 1., const unsigned CHUNKS_PER_LANE = 1);
template<class BUFFER> static void extract_boundary(char* host, const BUFFER * buffer, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS, const float ELEMS_PER_SITE = 1., const unsigned CHUNKS_PER_LANE = 1);

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

template<class BUFFER> static void extract_boundary(char* host, const BUFFER * buffer, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE)
{
	logger.debug() << "Extracting boundary. Offset: " << in_lane_offset << " - Elements per chunk: " << HALO_CHUNK_ELEMS;
	const unsigned NUM_LANES = buffer->get_lane_count();
	const unsigned STORAGE_TYPE_SIZE = buffer->get_storage_type_size();
	const unsigned CHUNK_STRIDE = get_vol4d(buffer->get_device()->get_mem_lattice_size()) * ELEMS_PER_SITE;

	for(size_t lane = 0; lane < NUM_LANES; ++lane) {
		logger.trace() << "Reading lane " << lane;
		size_t lane_offset = lane * buffer->get_lane_stride();
		for(size_t chunk = 0; chunk < CHUNKS_PER_LANE; ++chunk) {
			size_t host_offset = (lane * CHUNKS_PER_LANE + chunk) * HALO_CHUNK_ELEMS;
			size_t dev_offset = lane_offset + in_lane_offset + chunk * CHUNK_STRIDE;
			logger.trace() << "Chunk " << chunk << " - host offset: " << host_offset << " - device offset: " << dev_offset;
			buffer->dump_raw(&host[host_offset * STORAGE_TYPE_SIZE], HALO_CHUNK_ELEMS * STORAGE_TYPE_SIZE, dev_offset * STORAGE_TYPE_SIZE);
		}
	}
}

template<class BUFFER> static void send_halo(const BUFFER * buffer, const char* host, size_t in_lane_offset, size_t HALO_CHUNK_ELEMS, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE)
{
	logger.debug() << "Sending Halo. Offset: " << in_lane_offset << " - Elements per chunk: " << HALO_CHUNK_ELEMS;
	const unsigned NUM_LANES = buffer->get_lane_count();
	const unsigned STORAGE_TYPE_SIZE = buffer->get_storage_type_size();
	const unsigned CHUNK_STRIDE = get_vol4d(buffer->get_device()->get_mem_lattice_size()) * ELEMS_PER_SITE;

	for(size_t lane = 0; lane < NUM_LANES; ++lane) {
		logger.trace() << "Sending lane " << lane;
		size_t lane_offset = lane * buffer->get_lane_stride();
		for(size_t chunk = 0; chunk < CHUNKS_PER_LANE; ++chunk) {
			size_t host_offset = (lane * CHUNKS_PER_LANE + chunk) * HALO_CHUNK_ELEMS;
			size_t dev_offset = lane_offset + in_lane_offset + chunk * CHUNK_STRIDE;
			logger.trace() << "Chunk " << chunk << " - host offset: " << host_offset << " - device offset: " << dev_offset;
			buffer->load_raw(&host[host_offset * STORAGE_TYPE_SIZE], HALO_CHUNK_ELEMS * STORAGE_TYPE_SIZE, dev_offset * STORAGE_TYPE_SIZE);
		}
	}
}

template <typename T, class BUFFER> void hardware::buffers::update_halo_soa(std::vector<BUFFER*> buffers, const meta::Inputparameters& params, const float ELEMS_PER_SITE, const unsigned CHUNKS_PER_LANE)
{
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
		const unsigned HALO_ELEMS = HALO_SIZE * VOLSPACE * CHUNKS_PER_LANE;
		const unsigned HALO_CHUNK_ELEMS = HALO_SIZE * VOLSPACE;
		const unsigned VOL4D_LOCAL = get_vol4d(main_device->get_local_lattice_size()) * ELEMS_PER_SITE;

		logger.trace() << "GRID_SIZE: " << GRID_SIZE;
		logger.trace() << "HALO_SIZE: " << HALO_SIZE;
		logger.trace() << "eff. VOLSPACE: " << VOLSPACE;
		logger.trace() << "HALO ELEMS: " << HALO_ELEMS;
		logger.trace() << "HALO ELEMS per CHUNK: " << HALO_CHUNK_ELEMS;
		logger.trace() << "eff. VOL4D_LOCAL: " << VOL4D_LOCAL;
		logger.trace() << "CHUNKS_PER_LANE: " << CHUNKS_PER_LANE;

		std::vector<char*> upper_boundaries;
		upper_boundaries.reserve(num_buffers);
		std::vector<char*> lower_boundaries;
		lower_boundaries.reserve(num_buffers);
		for(size_t i = 0; i < num_buffers; ++i) {
			upper_boundaries.push_back(new char[HALO_ELEMS * sizeof(T)]);
			lower_boundaries.push_back(new char[HALO_ELEMS * sizeof(T)]);
		}

		// copy inside of boundaries to host
		for(size_t i = 0; i < num_buffers; ++i) {
			const auto buffer = buffers[i];
			logger.debug() << "Extracting data from buffer " << i;
			extract_boundary(upper_boundaries[i], buffer, VOL4D_LOCAL - HALO_CHUNK_ELEMS, HALO_CHUNK_ELEMS, ELEMS_PER_SITE, CHUNKS_PER_LANE);
			extract_boundary(lower_boundaries[i], buffer, 0, HALO_CHUNK_ELEMS, ELEMS_PER_SITE, CHUNKS_PER_LANE);
		}

		// copy data from host to halo (of coure getting what the neighbour stored
		for(size_t i = 0; i < num_buffers; ++i) {
			const auto buffer = buffers[i];
			logger.debug() << "Sending data to buffer " << i;
			// our lower halo is the upper bounary of our lower neighbour
			// its storage location is wrapped around to be the last chunk of data in our buffer, that is after local data and upper halo
			send_halo(buffer, upper_boundaries[lower_grid_neighbour(i, GRID_SIZE)], VOL4D_LOCAL + HALO_CHUNK_ELEMS, HALO_CHUNK_ELEMS, ELEMS_PER_SITE, CHUNKS_PER_LANE);
			// our upper halo is the lower bounary of our upper neighbour, it's stored right after our local data
			send_halo(buffer, lower_boundaries[upper_grid_neighbour(i, GRID_SIZE)], VOL4D_LOCAL, HALO_CHUNK_ELEMS, ELEMS_PER_SITE, CHUNKS_PER_LANE);
		}

		// clean up host
		for(size_t i = 0; i < num_buffers; ++i) {
			delete[] upper_boundaries[i];
			delete[] lower_boundaries[i];
		}
	}
}

#endif /* _HARDWARE_BUFFERS_PLAIN_ */
