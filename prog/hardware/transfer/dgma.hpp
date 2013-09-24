/** @file
 * Interface for the DirectGMA based transfer method
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _HARDWARE_TRANSFER_DGMA_HPP_
#define _HARDWARE_TRANSFER_DGMA_HPP_

#include <CL/cl_ext.h>

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
