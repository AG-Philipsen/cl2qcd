#ifndef _HARDWARE_BUFFERS_PLAIN_
#define _HARDWARE_BUFFERS_PLAIN_

#include "plain.hpp"

namespace hardware {

namespace buffers {

template <typename T, class BUFFER> void update_halo(std::vector<BUFFER*> buffers, const meta::Inputparameters& params, const float ELEMS_PER_SITE = 1.);

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

#endif /* _HARDWARE_BUFFERS_PLAIN_ */
