/** @file
 * Interface for buffer transfer methods
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

#ifndef _HARDWARE_TRANSFER_TRANSFER_HPP_
#define _HARDWARE_TRANSFER_TRANSFER_HPP_

#include "../buffers/buffer.hpp"
#include <memory>

namespace hardware {

class System;

/**
 * A transfer interface.
 *
 * Note that these classes are not stateless -> make sure you don't start with a second transfer before the first ist finished.
 * @todo get rid of this restriction.
 *
 * You might wonder why the transfer is split into multiple phases.
 * The major reason is that loading and dumping (aka. writing to the src and destination buffers) should usually happen in the same
 * command queue that is used by the kernels. Therefore you might first one all devices to dump their local data into transfers and
 * afterwards load the transfered data on each device to avoid serialization of the transfers.
 */
class Transfer {
	protected:
		Transfer(Device * src, Device * dest)
			: src_device(src), dest_device(dest) { };

	public:
		virtual ~Transfer() { };
		Transfer(Transfer const&) = delete;
		Transfer & operator=(Transfer const&) = delete;

		/**
		 * Load data into the transfer.
		 *
		 * Once the returned event is complete the src buffer may be modified without interfering with the data transfer.
		 */
		virtual SynchronizationEvent load(const hardware::buffers::Buffer* orig, const size_t *src_origin, const size_t *region, size_t src_row_pitch, size_t src_slice_pitch, const hardware::SynchronizationEvent& event) = 0;
		/**
		 * Send data to the destination device.
		 *
		 * The transfer might also haben in the load or dump phase, in this case this could be a noop returning a dummy, always complete event.
		 */
		virtual SynchronizationEvent transfer() = 0;
		/**
		 * Dump transferred data into the 
		 *
		 * Once the returned event is complete the src buffer may be modified without interfering with the data transfer.
		 */
		virtual SynchronizationEvent dump(const hardware::buffers::Buffer* dest, const size_t *dest_origin, const size_t *region, size_t dest_row_pitch, size_t dest_slice_pitch, const hardware::SynchronizationEvent& event) = 0;

		Device * get_src_device() const { return src_device; }
		Device * get_dest_device() const { return dest_device; }

	private:
		Device * const src_device;
		Device * const dest_device;
};

/**
 * Set up a transfer link from the source to the destination device.
 */
std::unique_ptr<Transfer> create_transfer(Device * src, Device * dest, hardware::System const & system);

}

#endif
