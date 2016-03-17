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

#include "gaugemomenta.hpp"
#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include "../buffers/halo_update.hpp"
#include "../code/gaugemomentum.hpp"
#include "../../geometry/parallelization.hpp"

hardware::lattices::Gaugemomenta::Gaugemomenta(const hardware::System& system):
system(system), buffers(allocate_buffers())
{}

hardware::lattices::Gaugemomenta::~Gaugemomenta()
{
    for(auto buffer: buffers) {
		delete buffer;
	}
}

const std::vector<const hardware::buffers::Gaugemomentum *> hardware::lattices::Gaugemomenta::get_buffers() const noexcept
{
	return buffers;
}

std::vector<const hardware::buffers::Gaugemomentum *> hardware::lattices::Gaugemomenta::allocate_buffers() const
{
	using hardware::buffers::Gaugemomentum;

	// only use device 0 for now
	auto devices = system.get_devices();
	std::vector<const Gaugemomentum*> buffers;
	buffers.reserve(devices.size());
	for(auto device: devices) {
		buffers.push_back(new Gaugemomentum(NDIM * (device->getLocalLatticeMemoryExtents().getLatticeVolume()), device)); //dont do calculations here!
	}
	return buffers;
}

void hardware::lattices::Gaugemomenta::update_halo_aos(const std::vector<const hardware::buffers::Gaugemomentum *> buffers, const hardware::System& system) const
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo<ae>(buffers, system, NDIM);
}

void hardware::lattices::Gaugemomenta::update_halo_soa(const std::vector<const hardware::buffers::Gaugemomentum *> buffers, const hardware::System& system) const
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(!buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo_soa<ae>(buffers, system, .5, 2 * NDIM);
}

void hardware::lattices::Gaugemomenta::update_halo() const
{
	if(buffers.size() > 1) { // for a single device this will be a noop
		// currently either all or none of the buffers must be SOA
		if(buffers[0]->is_soa()) {
			update_halo_soa(buffers, system);
		} else {
			update_halo_aos(buffers, system);
		}
	}
}

void hardware::lattices::Gaugemomenta::import(const ae * const host) const
{
	logger.trace() << "importing gaugemomenta";
	if(buffers.size() == 1) {
		auto device = buffers[0]->get_device();
		device->getGaugemomentumCode()->importGaugemomentumBuffer(buffers[0], host);
	} else {
		for(auto const buffer: buffers) {
			auto device = buffer->get_device();
			ae * mem_host = new ae[buffer->get_elements()];

//			//todo: put these calls into own fct.! With smart pointers?
			TemporalParallelizationHandlerLink tmp2(device->getGridPos(), device->getLocalLatticeExtents(), sizeof(ae), device->getHaloExtent());
			memcpy(&mem_host[tmp2.getMainPartIndex_destination()]  , &host[tmp2.getMainPartIndex_source()]  , tmp2.getMainPartSizeInBytes());
			memcpy(&mem_host[tmp2.getFirstHaloIndex_destination()] , &host[tmp2.getFirstHaloPartIndex_source()] , tmp2.getHaloPartSizeInBytes());
			memcpy(&mem_host[tmp2.getSecondHaloIndex_destination()], &host[tmp2.getSecondHaloPartIndex_source()], tmp2.getHaloPartSizeInBytes());

			device->getGaugemomentumCode()->importGaugemomentumBuffer(buffer, mem_host);

			delete[] mem_host;
		}
	}
	logger.trace() << "import complete";
}