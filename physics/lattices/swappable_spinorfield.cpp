/** @file
 * Implementation of the physics::lattices::SwappableSpinorfield class
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#include "swappable_spinorfield.hpp"

physics::lattices::SwappableSpinorfield::SwappableSpinorfield(const hardware::System& system, const physics::lattices::SpinorfieldParametersInterface & spinorfieldParametersInterface, const bool place_on_host) : Spinorfield(system, spinorfieldParametersInterface, place_on_host)
{
	// nothing to do
}

physics::lattices::SwappableSpinorfield::~SwappableSpinorfield()
{
	for(auto host_mem: swap) {
		delete[] host_mem;
	}
	swap.clear();
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
		spinor* host_mem = swap[i];
		buffers[i]->load(host_mem);
		delete[] host_mem;
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

std::vector<physics::lattices::Spinorfield *> physics::lattices::create_swappable_spinorfields(const hardware::System& system, const size_t n, physics::InterfacesHandler & interfacesHandler, const bool place_on_host)
{
	std::vector<Spinorfield *> fields;
	fields.reserve(n);

	for(size_t i = 0; i < n; ++i) {
		auto field = new SwappableSpinorfield(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>(), place_on_host);
		field->swap_out();
		fields.push_back(field);
	}

	return fields;
}

void physics::lattices::swap_out(const std::vector<Spinorfield*>& fields)
{
	logger.trace() << "Swapping out spinorfields.";
	for(auto field: fields) {
		auto swappable_field = dynamic_cast<Swappable*>(field);
		if(!swappable_field) {
			throw Print_Error_Message("swap_out was given a field that is not swappable", __FILE__, __LINE__);
		}
		swappable_field->swap_out();
	}
	logger.trace() << "Finished swapping out spinorfields.";
}

void physics::lattices::swap_in(const std::vector<Spinorfield*>& fields)
{
	logger.trace() << "Swapping in spinorfields.";
	for(auto field: fields) {
		auto swappable_field = dynamic_cast<Swappable*>(field);
		if(swappable_field) {
			swappable_field->swap_in();
		}
	}
	logger.trace() << "Finished swapping in spinorfields.";
}
