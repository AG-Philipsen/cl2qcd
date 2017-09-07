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

#ifndef _HARDWARE_LATTICES_SCALAR_
#define _HARDWARE_LATTICES_SCALAR_

#include <vector>
#include <functional>
#include "../system.hpp"
#include "../buffers/plain.hpp"
#include "../../executables/exceptions.h"
#include "../../meta/type_ops.hpp"

namespace hardware {

namespace lattices {

	template<typename SCALAR> static std::vector<const hardware::buffers::Plain<SCALAR>* > create_scalar_buffers(const hardware::System& system);

	template<typename SCALAR> class Scalar {
public:
	Scalar(const hardware::System& system) : system(system), buffers(create_scalar_buffers<SCALAR>(system)) { };

	Scalar& operator=(const Scalar&) = delete;
	Scalar(const Scalar&) = delete;
	Scalar() = delete;

	~Scalar();

	SCALAR get() const;

	void sum() const;

	SCALAR get_sum() const;
	
	void store(const SCALAR& val) const;
	
	const std::vector<const hardware::buffers::Plain<SCALAR> *> get_buffers() const noexcept;

private:
	const hardware::System& system;
	const std::vector<const hardware::buffers::Plain<SCALAR>* > buffers;
	};
}

/*
 * TEMPLATE IMPLEMENTATION
 */

template<typename SCALAR> static std::vector<const hardware::buffers::Plain<SCALAR>* > hardware::lattices::create_scalar_buffers(const hardware::System& system)
{
	using hardware::buffers::Plain;

	std::vector<const Plain<SCALAR> *> buffers;


	auto const devices = system.get_devices();
	for(auto device: devices) {
		cl_mem_flags buffer_flags = 0;
#ifdef CL_MEM_USE_PERSISTENT_MEM_AMD
		if(device->check_extension("cl_amd_device_persistent_memory")) {
			buffer_flags |= CL_MEM_USE_PERSISTENT_MEM_AMD;
		}
#endif

		buffers.push_back(new Plain<SCALAR>(1, device, false, buffer_flags));
	}

	return buffers;
}

template<typename SCALAR> hardware::lattices::Scalar<SCALAR>::~Scalar()
{
for(auto buffer: buffers) {
		delete buffer;
	}
}

template<typename SCALAR> SCALAR hardware::lattices::Scalar<SCALAR>::get() const
{
	// if this is a scalar we can read from any buffer
	auto buffer = buffers[0];
	SCALAR host_val;
	buffer->dump(&host_val);
	return host_val;
}

template<typename SCALAR> void hardware::lattices::Scalar<SCALAR>::sum() const
{
	// TODO make this run async
	size_t num_buffers = buffers.size();
	if(num_buffers > 1) {
		std::vector<std::unique_ptr<hardware::buffers::MappedBufferHandle>> handles(num_buffers);
		for(size_t i = 0; i < num_buffers; ++i) {
			handles[i] = buffers[i]->map();
		}

		std::vector<SCALAR*> mapped_ptrs(num_buffers);
		std::vector<hardware::SynchronizationEvent> events(num_buffers);
		for(size_t i = 0; i < num_buffers; ++i) {
			mapped_ptrs[i] = static_cast<SCALAR*>(handles[i]->get_mapped_ptr());
			events[i] = handles[i]->get_map_event();
		}
		hardware::wait(events);

		SCALAR res = mapped_ptrs[0][0] + mapped_ptrs[1][0];
		for(size_t i = 2; i < num_buffers; ++i) {
			res += mapped_ptrs[i][0];
		}
		for(size_t i = 0; i < num_buffers; ++i) {
			mapped_ptrs[i][0] = res;
		}
		logger.trace() << "Summed scalar: " << std::setprecision(16) << res;

		// buffers are automatically unmapped by leaving scope
	}
}

template<typename SCALAR> SCALAR hardware::lattices::Scalar<SCALAR>::get_sum() const
{
	size_t num_buffers = buffers.size();
	if(num_buffers > 1) {
		std::vector<std::unique_ptr<hardware::buffers::MappedBufferHandle>> handles(num_buffers);
		for(size_t i = 0; i < num_buffers; ++i) {
			handles[i] = buffers[i]->map(CL_MAP_READ);
		}

		std::vector<SCALAR*> mapped_ptrs(num_buffers);
		std::vector<hardware::SynchronizationEvent> events(num_buffers);
		for(size_t i = 0; i < num_buffers; ++i) {
			mapped_ptrs[i] = static_cast<SCALAR*>(handles[i]->get_mapped_ptr());
			events[i] = handles[i]->get_map_event();
		}
		hardware::wait(events);

		SCALAR res = mapped_ptrs[0][0] + mapped_ptrs[1][0];
		for(size_t i = 2; i < num_buffers; ++i) {
			res += mapped_ptrs[i][0];
		}
		logger.trace() << "Summed scalar: " << std::setprecision(16) << res;

		// buffers are automatically unmapped by leaving scope
		return res;
	} else {
		return get();
	}
}

template<typename SCALAR> void hardware::lattices::Scalar<SCALAR>::store(const SCALAR& val) const
{
	size_t num_buffers = buffers.size();
	if(num_buffers > 1) {
		std::vector<hardware::SynchronizationEvent> events(num_buffers);
		for(size_t i = 0; i < num_buffers; ++i) {
			events[i] = buffers[i]->load_async(&val);
		}
		hardware::wait(events);
	} else {
		buffers[0]->load(&val);
	}
}

template<typename SCALAR> const std::vector<const hardware::buffers::Plain<SCALAR> *> hardware::lattices::Scalar<SCALAR>::get_buffers() const noexcept
{
	return buffers;
}

}
#endif /* _HARDWARE_LATTICES_SCALAR_ */
