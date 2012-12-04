/** @file
 * Implementation of the physics::lattices::Spinorfield class
 */

#include "spinorfield.hpp"
#include "../../meta/util.hpp"
#include <cassert>

static std::vector<const hardware::buffers::Plain<spinor> *> allocate_buffers(hardware::System& system);

physics::lattices::Spinorfield::Spinorfield(hardware::System& system)
	: system(system), buffers(allocate_buffers(system))
{
}

static  std::vector<const hardware::buffers::Plain<spinor> *> allocate_buffers(hardware::System& system)
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

const std::vector<const physics::lattices::Spinorfield *> physics::lattices::create_spinorfields(hardware::System& system, const size_t n)
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
