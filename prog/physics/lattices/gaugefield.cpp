/** @file
 * Implementation of the physics::lattices::Gaugefield class
 */

#include "gaugefield.hpp"
#include "../../meta/util.hpp"
#include "../../logger.hpp"
#include "../../host_operations_gaugefield.h"

static std::vector<const hardware::buffers::SU3 *> allocate_buffers(hardware::System& system);

static void set_hot(std::vector<const hardware::buffers::SU3 *> buffers, physics::PRNG& prng);
static void set_cold(std::vector<const hardware::buffers::SU3 *> buffers);
static void set_cold(Matrixsu3 * field, size_t elems);
static void set_hot(Matrixsu3 * field, physics::PRNG& prng, size_t elems);

physics::lattices::Gaugefield::Gaugefield(hardware::System& system, physics::PRNG& prng, bool hot)
	: system(system), prng(prng), buffers(allocate_buffers(system))
{
	if(hot) {
		set_hot(buffers, prng);
	} else {
		set_cold(buffers);
	}
}

static  std::vector<const hardware::buffers::SU3 *> allocate_buffers(hardware::System& system)
{
	using hardware::buffers::SU3;

	// only use device 0 for now
	hardware::Device * device = system.get_devices().at(0);
	std::vector<const SU3 *> buffers;
	buffers.push_back(new SU3(meta::get_vol4d(system.get_inputparameters()) * 4, device));
	return buffers;
}

static void set_hot(std::vector<const hardware::buffers::SU3 *> buffers, physics::PRNG& prng)
{
	using hardware::Device;

for(auto buffer: buffers) {
		size_t elems = buffer->get_elements();
		Matrixsu3 * tmp = new Matrixsu3[elems];
		set_hot(tmp, prng, elems);
		Device * device = buffer->get_device();
		device->get_gaugefield_code()->importGaugefield(buffer, tmp);
		device->synchronize();
		delete[] tmp;
	}
}

static void set_cold(std::vector<const hardware::buffers::SU3 *> buffers)
{
	using hardware::Device;

for(auto buffer: buffers) {
		size_t elems = buffer->get_elements();
		Matrixsu3 * tmp = new Matrixsu3[elems];
		set_cold(tmp, elems);
		Device * device = buffer->get_device();
		device->get_gaugefield_code()->importGaugefield(buffer, tmp);
		device->synchronize();
		delete[] tmp;
	}
}

void set_cold(Matrixsu3 * field, size_t elems)
{
	for(size_t i = 0; i < elems; ++i) {
		field[i] = unit_matrixsu3();
	}
}

void set_hot(Matrixsu3 * field, physics::PRNG& prng, size_t elems)
{
	for(size_t i = 0; i < elems; ++i) {
		field[i] = random_matrixsu3();
	}
}
