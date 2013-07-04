/** @file
 * Interface for the simple opencl transfer method
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _HARDWARE_TRANSFER_ASYNC_OCL_COPY_HPP_
#define _HARDWARE_TRANSFER_ASYNC_OCL_COPY_HPP_

#include "transfer.hpp"
#include <map>
#include <memory>

namespace hardware {

namespace transfer {

/**
 * A transfer method using normal opencl buffer copies but transferring independent of the base command streams.
 */
class AsyncOclCopy : public Transfer {

	public:
		AsyncOclCopy(hardware::Device * from, hardware::Device * to, hardware::System const & context);
		virtual ~AsyncOclCopy();

		SynchronizationEvent load(const hardware::buffers::Buffer* orig, const size_t *src_origin, const size_t *region, size_t src_row_pitch, size_t src_slice_pitch, const hardware::SynchronizationEvent& event) override;
		SynchronizationEvent transfer() override;
		SynchronizationEvent dump(const hardware::buffers::Buffer* dest, const size_t *dest_origin, const size_t *region, size_t dest_row_pitch, size_t dest_slice_pitch, const hardware::SynchronizationEvent& event) override;

	private:
		std::map<size_t, std::unique_ptr<hardware::buffers::Buffer>> src_cache;
		std::map<size_t, std::unique_ptr<hardware::buffers::Buffer>> dest_cache;
		hardware::SynchronizationEvent load_event;
		hardware::SynchronizationEvent transfer_event; 
		hardware::SynchronizationEvent dump_event; 
		size_t active_size;
		cl_command_queue transfer_queue;

		hardware::buffers::Buffer * get_src_cache(size_t const size);
		hardware::buffers::Buffer * get_dest_cache(size_t const size);
};

}

}

#endif
