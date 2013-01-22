/** @file
 * Implementation of the physics::lattices::Gaugemomenta class
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "gaugemomenta.hpp"
#include "../../meta/util.hpp"

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
	buffers.push_back(new Gaugemomentum(NDIM * meta::get_vol4d(system.get_inputparameters()), device));
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
