/** @file
 * Declaration and implementation of the physics::lattices::Scalar template
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _PHYSICS_LATTICES_SCALAR_
#define _PHYSICS_LATTICES_SCALAR_

#include <vector>
#include <functional>
#include "../../hardware/system.hpp"
#include "../../hardware/buffers/plain.hpp"
#include "../../exceptions.h"

namespace physics {
	namespace lattices {

		/**
		 * Utility method to create the buffers for a Scalar
		 */
		template<typename SCALAR> static std::vector<const hardware::buffers::Plain<SCALAR>* > create_scalar_buffers(const hardware::System& system);

		/**
		 * Allows usage of scalar datatypes.
		 */
		template<typename SCALAR> class Scalar {

		public:
			/**
			 * Create a new scalar
			 */
			Scalar(const hardware::System& system) : system(system), buffers(create_scalar_buffers<SCALAR>(system)) { };

			/**
			 * Cleanup
			 */
			~Scalar();

			/*
			 * Don't allow copies
			 */
			Scalar& operator=(const Scalar&) = delete;
			Scalar(const Scalar&) = delete;
			Scalar() = delete;

			/**
			 * Retrieve the value
			 */
			SCALAR get() const;

			/**
			 * Sum up the data of all buffers, getting the scalar into a consistent state again.
			 *
			 * On single-device systems this is a noop.
			 */
			void sum() const;

			/**
			 * Store a value
			 */
			void store(const SCALAR& val) const;

			/**
			 * Get the buffers containing the scalar on the devices.
			 *
			 * @todo with multi-device we might have to give up on the noexcept part
			 */
			const std::vector<const hardware::buffers::Plain<SCALAR> *> get_buffers() const noexcept;

		private:
			/**
			 * The system the scalar operates on
			 */
			const hardware::System& system;

			/**
			 * The buffers of the scalar
			 */
			const std::vector<const hardware::buffers::Plain<SCALAR>* > buffers;
		};
	}
}

/*
 * TEMPLATE IMPLEMENTATION
 */

template<typename SCALAR> static std::vector<const hardware::buffers::Plain<SCALAR>* > physics::lattices::create_scalar_buffers(const hardware::System& system)
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

template<typename SCALAR> physics::lattices::Scalar<SCALAR>::~Scalar()
{
for(auto buffer: buffers) {
		delete buffer;
	}
}

template<typename SCALAR> SCALAR physics::lattices::Scalar<SCALAR>::get() const
{
	// if this is a scalar we can read from any buffer
	auto buffer = buffers[0];
	SCALAR host_val;
	buffer->dump(&host_val);
	return host_val;
}

template<typename SCALAR> void physics::lattices::Scalar<SCALAR>::sum() const
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

template<typename SCALAR> void physics::lattices::Scalar<SCALAR>::store(const SCALAR& val) const
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

template<typename SCALAR> const std::vector<const hardware::buffers::Plain<SCALAR> *> physics::lattices::Scalar<SCALAR>::get_buffers() const noexcept
{
	return buffers;
}

#endif /* _PHYSICS_LATTICES_SCALAR_ */
