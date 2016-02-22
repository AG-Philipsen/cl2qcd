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

#include "gaugefield.hpp"
#include "../../host_functionality/logger.hpp"
#include "../../hardware/device.hpp"
#include "../../hardware/buffers/halo_update.hpp"
#include "../../hardware/code/gaugefield.hpp"
#include "../../geometry/parallelization.hpp"

hardware::lattices::Gaugefield::Gaugefield(const hardware::System& system):
system(system), buffers(allocate_buffers()),unsmeared_buffers()
{}

hardware::lattices::Gaugefield::~Gaugefield()
{
	release_buffers(&buffers);
	release_buffers(&unsmeared_buffers);
}

std::vector<const hardware::buffers::SU3 *> hardware::lattices::Gaugefield::allocate_buffers() const
{
	using hardware::buffers::SU3;

	std::vector<const SU3 *> buffers;

	auto const devices = system.get_devices();
	for(auto device: devices) 
	{
		buffers.push_back(new SU3(device->getLocalLatticeMemoryExtents().getLatticeVolume() * 4, device)); //todo: do not calculate here!
	}
	return buffers;
}

void hardware::lattices::Gaugefield::release_buffers(std::vector<const hardware::buffers::SU3 *>* buffers) const
{
	for(auto buffer: *buffers) 
	{
		delete buffer;
	}
	buffers->clear();
}

void hardware::lattices::Gaugefield::send_gaugefield_to_buffers(const std::vector<const hardware::buffers::SU3 *> buffers, const Matrixsu3 * const gf_host) const
{
	logger.trace() << "importing gaugefield";
// 	if(buffers.size() == 1) {
// 		auto device = buffers[0]->get_device();
// 		device->getGaugefieldCode()->importGaugefield(buffers[0], gf_host);
// 		device->synchronize();
// 	} else {
		for(auto const buffer: buffers) {
			logger.fatal() << buffers.size();
			auto device = buffer->get_device();
			TemporalParallelizationHandlerLink tmp2(device->getGridPos(), device->getLocalLatticeExtents(), sizeof(Matrixsu3), device->getHaloExtent());

			if(buffers.size() == 1) device->getGaugefieldCode()->importGaugefield(buffer, gf_host);
			else{
				Matrixsu3 * mem_host = new Matrixsu3[buffer->get_elements()];
//				//todo: put these calls into own fct.! With smart pointers?
				memcpy(&mem_host[tmp2.getMainPartIndex_destination()]  , &gf_host[tmp2.getMainPartIndex_source()]  , tmp2.getMainPartSizeInBytes());
				memcpy(&mem_host[tmp2.getFirstHaloIndex_destination()] , &gf_host[tmp2.getFirstHaloPartIndex_source()] , tmp2.getHaloPartSizeInBytes());
				memcpy(&mem_host[tmp2.getSecondHaloIndex_destination()], &gf_host[tmp2.getSecondHaloPartIndex_source()], tmp2.getHaloPartSizeInBytes());

				device->getGaugefieldCode()->importGaugefield(buffer, mem_host);
				delete[] mem_host;
			}
			device->synchronize();
		}
// 	}
	logger.trace() << "import complete";
}

void hardware::lattices::Gaugefield::fetch_gaugefield_from_buffers(const std::vector<const hardware::buffers::SU3 *> buffers, Matrixsu3 * const gf_host) const
{
	logger.trace() << "fetching gaugefield";
	if(buffers.size() == 1) {
		auto device = buffers[0]->get_device();
		device->getGaugefieldCode()->exportGaugefield(gf_host, buffers[0]);
		device->synchronize();
	} else {
		for(auto const buffer: buffers) {
			// fetch local part for each device
			auto device = buffer->get_device();
			Matrixsu3 * mem_host = new Matrixsu3[buffer->get_elements()];

			device->getGaugefieldCode()->exportGaugefield(mem_host, buffer);

			TemporalParallelizationHandlerLink tmp2(device->getGridPos(), device->getLocalLatticeExtents(), sizeof(Matrixsu3), device->getHaloExtent());
			memcpy(&gf_host[tmp2.getMainPartIndex_source()]  , &mem_host[tmp2.getMainPartIndex_destination()]  , tmp2.getMainPartSizeInBytes());

			delete[] mem_host;
		}
	}
}

void hardware::lattices::Gaugefield::update_halo_aos(const std::vector<const hardware::buffers::SU3 *> buffers, const hardware::System& system) const
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo<Matrixsu3>(buffers, system, NDIM);
}

void hardware::lattices::Gaugefield::update_halo_soa(const std::vector<const hardware::buffers::SU3 *> buffers, const hardware::System& system) const
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(!buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo_soa<Matrixsu3>(buffers, system, .5, 2 * NDIM);
}