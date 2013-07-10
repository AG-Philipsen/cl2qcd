/** @file
 * Implementation of the transfer method factory
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "transfer.hpp"
#include "ocl_copy.hpp"
#include "async_ocl_copy.hpp"
#include "dgma.hpp"
#include "../../logger.hpp"
#include <array>
#include "../../klepsydra/klepsydra.hpp"
#include "../device.hpp"

namespace {
/**
 * Benchmark the given Transfer link for a buffer of 4 MiB.
 *
 * @return Performance in GB/s
 */
float benchmark(hardware::Transfer *);
}

std::unique_ptr<hardware::Transfer> hardware::create_transfer(hardware::Device * const src, hardware::Device * const dest, hardware::System const & system)
{
#ifdef ASYNC_HALO_UPDATES
	std::unique_ptr<hardware::Transfer> standard_transfer(new transfer::AsyncOclCopy(src, dest, system));
#else
	std::unique_ptr<hardware::Transfer> standard_transfer(new transfer::OclCopy(src, dest));
#endif
// required headers for DirectGMA might be missing
// in that case it is unavailable
#ifndef CL_MEM_BUS_ADDRESSABLE_AMD
	return std::move(standard_transfer);
#else
	try {
		std::unique_ptr<hardware::Transfer> dgma_transfer(new transfer::DirectGMA(src, dest, system));
		// DirectGMA is only fast for two GPUs on the same PCIe
		// compare to standard method and use faster one
		auto const standard_performance = benchmark(standard_transfer.get());
		auto const dgma_performance = benchmark(dgma_transfer.get());

		logger.debug() << "Standard Transfer Performance: " << standard_performance << " GB/s - DirectGMA Performance: " << dgma_performance << " GB/s";

		if(standard_performance > dgma_performance) {
			logger.warn() << "DirectGMA performs badly from " << src->get_grid_pos() << " to " << dest->get_grid_pos() << ". Falling back to standard transfer methods.";
			return std::move(standard_transfer);
		} else {
			logger.debug() << "Using DirectGMA from " << src->get_grid_pos() << " to " << dest->get_grid_pos() << ".";
			return std::move(dgma_transfer);
		}

	} catch(hardware::transfer::DGMAUnsupported) {
		logger.warn() << "DirectGMA is unavailable from " << src->get_grid_pos() << " to " << dest->get_grid_pos() << ". Falling back to standard transfer methods.";
		return std::move(standard_transfer);
	}
#endif
}

namespace {

float benchmark(hardware::Transfer * const transfer)
{
	using hardware::buffers::Buffer;
	using hardware::SynchronizationEvent;

	auto const size = 4u * 1024u *1024u;
	auto const iterations = 25;

	Buffer src(size, transfer->get_src_device());
	Buffer dest(size, transfer->get_dest_device());

	// warmup
	size_t const offset[3] = {0, 0, 0};
	size_t const layout[3] = {size, 1, 1};
	transfer->load(&src, offset, layout, 0, 0, SynchronizationEvent());
	transfer->transfer();
	auto event = transfer->dump(&dest, offset, layout, 0, 0, SynchronizationEvent());

	event.wait();
	klepsydra::Monotonic timer;

	for(auto i = 0; i < iterations - 1; ++i) {
		transfer->load(&src, offset, layout, 0, 0, event);
		transfer->transfer();
		event = transfer->dump(&dest, offset, layout, 0, 0, SynchronizationEvent());
	}

	event.wait();
	auto const elapsed_mus = timer.getTime();

	return (size * iterations / elapsed_mus / 1e3f);
}

}
