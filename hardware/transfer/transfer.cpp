/** @file
 * Implementation of the transfer method factory
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

#include "transfer.hpp"
#include "ocl_copy.hpp"
#include "async_ocl_copy.hpp"
#include "dgma.hpp"
#include "../../host_functionality/logger.hpp"
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
			logger.warn() << "DirectGMA performs badly from " << src->getGridPos() << " to " << dest->getGridPos() << ". Falling back to standard transfer methods.";
			return std::move(standard_transfer);
		} else {
			logger.debug() << "Using DirectGMA from " << src->getGridPos() << " to " << dest->getGridPos() << ".";
			return std::move(dgma_transfer);
		}

	} catch(hardware::transfer::DGMAUnsupported) {
		logger.warn() << "DirectGMA is unavailable from " << src->getGridPos() << " to " << dest->getGridPos() << ". Falling back to standard transfer methods.";
		return std::move(standard_transfer);
	}
#endif
}

namespace {

float benchmark(hardware::Transfer * const transfer)
{
	using hardware::buffers::Buffer;
	using hardware::SynchronizationEvent;

	auto const size = 4u * 1024u * 1024u;
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
