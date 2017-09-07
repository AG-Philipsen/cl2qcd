/*
 * Copyright 2016 Francesca Cuteri
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

#include "spinorfield.hpp"
#include "../device.hpp"
#include "../buffers/halo_update.hpp"
#include "../code/fermions.hpp"
#include "../code/spinors.hpp"
#include "../../geometry/parallelization.hpp"

hardware::lattices::Spinorfield::Spinorfield(const hardware::System& system, const bool place_on_host)
	: system(system), buffers(allocate_buffers()), place_on_host(place_on_host)
{}

std::vector<const hardware::buffers::Plain<spinor> *> hardware::lattices::Spinorfield::allocate_buffers()
{
	using hardware::buffers::Plain;

	std::vector<const Plain<spinor>*> buffers;
	for(auto device: system.get_devices()) {
		buffers.push_back(new Plain<spinor>(hardware::code::get_spinorfieldsize(device->getLocalLatticeMemoryExtents()), device, place_on_host));
	}
	return buffers;
}

hardware::lattices::Spinorfield::~Spinorfield()
{
    clear_buffers();
}

void hardware::lattices::Spinorfield::clear_buffers()
{
for(auto buffer: buffers) {
		delete buffer;
	}
	buffers.clear();
}

void hardware::lattices::Spinorfield::fill_buffers()
{
	if(buffers.size() != 0) {
		return;
	}

	buffers = allocate_buffers();
}

const std::vector<const hardware::buffers::Plain<spinor> *> hardware::lattices::Spinorfield::get_buffers() const noexcept
{
	return buffers;
}

void hardware::lattices::Spinorfield::update_halo() const
{
	hardware::buffers::update_halo<spinor>(buffers, system);
}

void hardware::lattices::Spinorfield::import(const spinor * const host) const
{
	logger.trace() << "importing spinorfield";
	if(buffers.size() == 1) {
		buffers[0]->load(host);
	} else {
		for(auto const buffer: buffers) {
			auto device = buffer->get_device();

//			//todo: put these calls into own fct.!
			TemporalParallelizationHandlerNonLink tmp2(device->getGridPos(), device->getLocalLatticeExtents(), sizeof(spinor), device->getHaloExtent());
			buffer->load( &host[tmp2.getMainPartIndex_source()] , tmp2.getMainPartSize());
			buffer->load( &host[tmp2.getFirstHaloPartIndex_source()] , tmp2.getHaloPartSize(), tmp2.getMainPartSize());
			buffer->load( &host[tmp2.getSecondHaloPartIndex_source()], tmp2.getHaloPartSize(), tmp2.getMainPartSize() + tmp2.getHaloPartSize());

		}
	}
	logger.trace() << "import complete";
}
