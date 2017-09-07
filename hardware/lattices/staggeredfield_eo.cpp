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

#include "staggeredfield_eo.hpp"
#include "../../hardware/code/spinors.hpp" //For hardware::code::get_eoprec_spinorfieldsize()

hardware::lattices::Staggeredfield_eo::Staggeredfield_eo(const hardware::System& system)
	: system(system), buffers(allocate_buffers())
{}

std::vector<const hardware::buffers::SU3vec *> hardware::lattices::Staggeredfield_eo::allocate_buffers()
{
	using hardware::buffers::SU3vec;

	auto devices = system.get_devices();
	std::vector<const SU3vec*> buffers;
	buffers.reserve(devices.size());
	for(auto device: devices) {
		buffers.push_back(new SU3vec(hardware::code::get_eoprec_spinorfieldsize(device->getLocalLatticeMemoryExtents()), device));
	}
	return buffers;
}

hardware::lattices::Staggeredfield_eo::~Staggeredfield_eo()
{
    for(auto buffer: buffers)
        delete buffer;
}

hardware::lattices::Staggeredfield_eo::Staggeredfield_eo(hardware::lattices::Staggeredfield_eo&& movedFrom)
    : system(std::move(movedFrom.system)), //use move even though a copy is done (it is a const reference!)
      buffers(std::move(movedFrom.buffers))
{
    (const_cast<std::vector<const hardware::buffers::SU3vec *>&>(movedFrom.buffers)).clear();
}

const std::vector<const hardware::buffers::SU3vec *> hardware::lattices::Staggeredfield_eo::get_buffers() const noexcept
{
	return buffers;
}

void hardware::lattices::Staggeredfield_eo::update_halo() const
{
	throw Print_Error_Message("Halo update is not implemented, yet.", __FILE__, __LINE__);
}

void hardware::lattices::Staggeredfield_eo::import(const su3vec * const host) const
{
	logger.trace() << "importing staggeredfield_eo";
	if(buffers.size() == 1) {
		buffers[0]->load(host);
	}
	else
	{
		throw Print_Error_Message("Import not implemented for multi device staggeredfield_eo!", __FILE__, __LINE__);
	}
	logger.trace() << "import complete";
}
