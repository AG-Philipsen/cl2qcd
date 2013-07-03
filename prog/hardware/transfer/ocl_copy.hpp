/** @file
 * Interface for the simple opencl transfer method
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _HARDWARE_TRANSFER_OCL_COPY_HPP_
#define _HARDWARE_TRANSFER_OCL_COPY_HPP_

#include "transfer.hpp"

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
		hardware::buffers::Buffer * transfer_buffer;
		hardware::SynchronizationEvent load_event;
		hardware::SynchronizationEvent dump_event;

		void ensure_proper_transfer_buffer_size(const size_t * region);
};

}

}

#endif
