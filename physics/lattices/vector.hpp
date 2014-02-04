/** @file
 * Declaration and implementation of the physics::lattices::Vector template
 *
 * Copyright (c) 2014 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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

#ifndef _PHYSICS_LATTICES_VECTOR_
#define _PHYSICS_LATTICES_VECTOR_

#include <vector>
#include <functional>
#include "../../hardware/system.hpp"
#include "../../hardware/buffers/plain.hpp"
#include "../../executables/exceptions.h"

namespace physics {
	namespace lattices {

		/**
		 * Utility method to create the buffers for a Vector
		 */
		template<typename SCALAR> static std::vector<const hardware::buffers::Plain<SCALAR>* > create_vector_buffers(const size_t N, const hardware::System& system);

		/**
		 * Allows usage of scalar datatypes.
		 */
		/* Note: here the size of the vector could be either specified as non type template
		 *       parameter or as constructor standard parameter. The difference is that in
		 *       the former case, at compile-time, a different instatiation of the class
		 *       is created for each size used, while in the latter case the compiler behaves
		 *       as usual. If we needed to create an array or a valarray with the size,
		 *       the first approach would be mandatory, then. Here we will use the second one.
		 */
		template<typename SCALAR> class Vector {

		public:
			/**
			 * Create a new Vector
			 */
			Vector(const size_t N, const hardware::System& system) : N(N), system(system), buffers(create_vector_buffers<SCALAR>(N, system)) { };

			/**
			 * Cleanup
			 */
			~Vector();

			/*
			 * Don't allow copies
			 */
			Vector& operator=(const Vector&) = delete;
			Vector(const Vector&) = delete;
			Vector() = delete;

			/**
			 * Retrieve the values
			 */
			std::vector<SCALAR> get() const;

			/**
			 * Store values
			 */
			void store(const std::vector<SCALAR> val) const;
			
			/**
			 * Retrieve the number of elements of each buffer
			 */
			size_t get_vector_size() const noexcept;

			/**
			 * Get the buffers containing the vectors on the devices.
			 *
			 * @todo with multi-device we might have to give up on the noexcept part
			 */
			const std::vector<const hardware::buffers::Plain<SCALAR> *> get_buffers() const noexcept;

		private:
			/**
			 * The number of elements of the Vector
			 */
			const size_t N;
			
			/**
			 * The system the Vector operates on
			 */
			const hardware::System& system;

			/**
			 * The buffers of the Vector
			 */
			const std::vector<const hardware::buffers::Plain<SCALAR>* > buffers;
		};
	}
}

/*
 * TEMPLATE IMPLEMENTATION
 */

template<typename SCALAR> static std::vector<const hardware::buffers::Plain<SCALAR>* > physics::lattices::create_vector_buffers(const size_t N, const hardware::System& system)
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

		buffers.push_back(new Plain<SCALAR>(N, device, false, buffer_flags));
	}

	return buffers;
}

template<typename SCALAR> physics::lattices::Vector<SCALAR>::~Vector()
{
for(auto buffer: buffers) {
		delete buffer;
	}
}

template<typename SCALAR> std::vector<SCALAR> physics::lattices::Vector<SCALAR>::get() const
{
	// We can read from any buffer
	auto buffer = buffers[0];
	std::vector<SCALAR> host_vec(N);
	buffer->dump(&host_vec[0]);
	return host_vec;
}

template<typename SCALAR> void physics::lattices::Vector<SCALAR>::store(const std::vector<SCALAR> vec) const
{
	if(vec.size() != N)
	  throw Invalid_Parameters("The vector given to Vector::store has the wrong size!", N, vec.size());
	size_t num_buffers = buffers.size();
	if(num_buffers > 1) {
		std::vector<hardware::SynchronizationEvent> events(num_buffers);
		for(size_t i = 0; i < num_buffers; ++i) {
			events[i] = buffers[i]->load_async(&vec[0]);
		}
		hardware::wait(events);
	} else {
		buffers[0]->load(&vec[0]);
	}
}

template<typename SCALAR> const std::vector<const hardware::buffers::Plain<SCALAR> *> physics::lattices::Vector<SCALAR>::get_buffers() const noexcept
{
	return buffers;
}

template<typename SCALAR> size_t physics::lattices::Vector<SCALAR>::get_vector_size() const noexcept
{
	return N;
}

#endif /* _PHYSICS_LATTICES_VECTOR_ */
