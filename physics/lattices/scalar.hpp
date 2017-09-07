/** @file
 * Declaration and implementation of the physics::lattices::Scalar template
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

#ifndef _PHYSICS_LATTICES_SCALAR_
#define _PHYSICS_LATTICES_SCALAR_

#include "../../hardware/lattices/scalar.hpp"

namespace physics {
	namespace lattices {
		/**
		 * Allows usage of scalar datatypes.
		 */
		template<typename SCALAR> class Scalar {

		public:
			/**
			 * Create a new scalar
			 */
			Scalar(const hardware::System& system) : system(system), scalar(system) { };

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
			 * Return the sum of all scalars.
			 *
			 * @warning Does not put the scalar into a consistent state. Manually store afterwards of this is desired.
			 */
			SCALAR get_sum() const;

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
			const hardware::System& system;
			hardware::lattices::Scalar<SCALAR> scalar;
		};
	}
}

template<typename SCALAR> physics::lattices::Scalar<SCALAR>::~Scalar()
{}

template<typename SCALAR> SCALAR physics::lattices::Scalar<SCALAR>::get() const
{
	return scalar.get();
}

template<typename SCALAR> void physics::lattices::Scalar<SCALAR>::sum() const
{
	scalar.sum();
}

template<typename SCALAR> SCALAR physics::lattices::Scalar<SCALAR>::get_sum() const
{
	return scalar.get_sum();
}

template<typename SCALAR> void physics::lattices::Scalar<SCALAR>::store(const SCALAR& val) const
{
	scalar.store(val);
}

template<typename SCALAR> const std::vector<const hardware::buffers::Plain<SCALAR> *> physics::lattices::Scalar<SCALAR>::get_buffers() const noexcept
{
	return scalar.get_buffers();
}

#endif /* _PHYSICS_LATTICES_SCALAR_ */
