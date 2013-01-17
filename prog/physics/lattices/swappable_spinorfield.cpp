/** @file
 * Implementation of the physics::lattices::SwappableSpinorfield class
 */

#include "swappable_spinorfield.hpp"

physics::lattices::SwappableSpinorfield::SwappableSpinorfield(const hardware::System& system, const bool place_on_host) : Spinorfield(system, place_on_host)
{
	// nothing to do
}

physics::lattices::SwappableSpinorfield::~SwappableSpinorfield()
{
	// nothing to do
}

void physics::lattices::SwappableSpinorfield::swap_in()
{
	if(get_buffers().size() != 0) {
		return;
	}

	fill_buffers();

	auto buffers = get_buffers();
	size_t num_bufs = buffers.size();
	if(num_bufs != swap.size()) {
		throw Print_Error_Message("Spinorfield does not have same structure as its own swap, something went really wrong!", __FILE__, __LINE__);
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		buffers[i]->load(swap[i]);
	}
	swap.clear();
}

void physics::lattices::SwappableSpinorfield::swap_out()
{
	auto buffers = get_buffers();
	if(buffers.size() == 0) {
		return;
	}

	for(auto buffer: buffers) {
		spinor* host_mem = new spinor[buffer->get_elements()];
		buffer->dump(host_mem);
		swap.push_back(host_mem);
	}

	clear_buffers();
}
