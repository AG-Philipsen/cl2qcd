/** @file
 * Implementation of the physics::lattices::Spinorfield_eo class
 */

#include "spinorfield_eo.hpp"
#include "../../meta/util.hpp"
#include <cassert>

static std::vector<const hardware::buffers::Spinor *> allocate_buffers(hardware::System& system);

physics::lattices::Spinorfield_eo::Spinorfield_eo(hardware::System& system)
	: system(system), buffers(allocate_buffers(system))
{
}

static  std::vector<const hardware::buffers::Spinor *> allocate_buffers(hardware::System& system)
{
	using hardware::buffers::Spinor;

	// only use device 0 for now
	hardware::Device * device = system.get_devices().at(0);
	std::vector<const Spinor*> buffers;
	buffers.push_back(new Spinor(meta::get_eoprec_spinorfieldsize(system.get_inputparameters()), device));
	return buffers;
}

physics::lattices::Spinorfield_eo::~Spinorfield_eo()
{
for(auto buffer: buffers) {
		delete buffer;
	}
}
