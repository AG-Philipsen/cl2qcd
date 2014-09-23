/** @file
 * Interface for the simple opencl transfer method
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

#ifndef _HARDWARE_TRANSFER_OCL_COPY_HPP_
#define _HARDWARE_TRANSFER_OCL_COPY_HPP_

#include "transfer.hpp"
#include <map>
#include <memory>

namespace hardware {

namespace transfer {

/**
 * A transfer method using normal opencl buffer copies.
 */
class OclCopy : public Transfer {

public:
	OclCopy(hardware::Device * from, hardware::Device * to);
	virtual ~OclCopy();

	SynchronizationEvent load(const hardware::buffers::Buffer* orig, const size_t *src_origin, const size_t *region, size_t src_row_pitch, size_t src_slice_pitch, const hardware::SynchronizationEvent& event) override;
	SynchronizationEvent transfer() override;
	SynchronizationEvent dump(const hardware::buffers::Buffer* dest, const size_t *dest_origin, const size_t *region, size_t dest_row_pitch, size_t dest_slice_pitch, const hardware::SynchronizationEvent& event) override;

private:
	std::map<size_t, std::unique_ptr<hardware::buffers::Buffer>> transfer_buffers;
	hardware::SynchronizationEvent load_event; // TODO these events could be stored with each buffer
	hardware::SynchronizationEvent dump_event; // TODO these events could be stored with each buffer

	hardware::buffers::Buffer * get_transfer_buffer(const size_t * region);
};

}

}

#endif
