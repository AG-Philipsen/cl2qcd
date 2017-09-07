/** @file
 * Interface for the DirectGMA based transfer method
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#ifndef _HARDWARE_TRANSFER_DGMA_HPP_
#define _HARDWARE_TRANSFER_DGMA_HPP_

//TODO: This header should not be comipled on Apple since DGMA won't be available.
//      Avoid ifdef in this case.
#ifdef __APPLE__
#include <OpenCL/cl_ext.h>
#else
#include <CL/cl_ext.h>
#endif

#ifdef CL_MEM_BUS_ADDRESSABLE_AMD // make sure definitions for DGMA are available

#include "transfer.hpp"
#include <map>
#include <memory>

namespace hardware {

namespace transfer {

class DGMAGhostBuffer;

/**
 * This class is thrown to indicate that one of the involved devices does not support DirectGMA.
 */
struct DGMAUnsupported { };

/**
 * A transfer method using normal opencl buffer copies but transferring independent of the base command streams.
 */
class DirectGMA : public Transfer {

public:
	/**
	 * Create a DirectGMA transfer link.
	 *
	 * \throw DGMAUnsupported if DirectGMA is not supported by all involved devices.
	 */
	DirectGMA(hardware::Device * from, hardware::Device * to, hardware::System const & context);
	virtual ~DirectGMA();

	SynchronizationEvent load(const hardware::buffers::Buffer* orig, const size_t *src_origin, const size_t *region, size_t src_row_pitch, size_t src_slice_pitch, const hardware::SynchronizationEvent& event) override;
	SynchronizationEvent transfer() override;
	SynchronizationEvent dump(const hardware::buffers::Buffer* dest, const size_t *dest_origin, const size_t *region, size_t dest_row_pitch, size_t dest_slice_pitch, const hardware::SynchronizationEvent& event) override;

private:
	std::unique_ptr<hardware::buffers::Buffer> src_cache;
	std::unique_ptr<hardware::transfer::DGMAGhostBuffer> ghost;
	std::unique_ptr<hardware::buffers::Buffer> dest_cache;
	hardware::SynchronizationEvent load_event;
	hardware::SynchronizationEvent transfer_event;
	hardware::SynchronizationEvent dump_event;
	cl_command_queue transfer_queue;
	size_t active_size;
	hardware::System const & system;

	/**
	 * Make sure that the transfer buffers have the proper size.
	 */
	void ensure_buffers(size_t const size);
};

}

}

#endif

#endif
