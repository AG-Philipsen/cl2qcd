/** @file
 * PRNG PRNG unit declaration
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#ifndef _PHYSICS_PRNG_
#define _PHYSICS_PRNG_

#include "../hardware/system.hpp"
#include "../hardware/buffers/prng_buffer.hpp"

/**
 * This package contains the actual "business" logic of the library.
 */
namespace physics {

	/**
	 * The PRNG PRNG for host and device.
	 *
	 * WARNING: Must only be used as a singleton!
	 */
	class PRNG {
		public:
			/**
			 * Initialize the ranlux instance
			 */
			PRNG(const hardware::System& system);

			~PRNG();

			/*
			 * non-copyable
			 */
			PRNG& operator=(const PRNG&) = delete;
			PRNG(const PRNG&) = delete;
			PRNG() = delete;

			/**
			 * Get random numbers
			 */
			double get_double() noexcept;

			/**
			 * Get the buffers containing the random state on the device.
			 */
			const std::vector<const hardware::buffers::PRNGBuffer*> get_buffers() const noexcept;

			/**
			 * Store the current PRNG state to a file of the given name.
			 *
			 * \param filename The name of the file to store the random state in
			 */
			void store(const std::string filename) const;

		private:
			/**
			 * Reference to the PRNG Buffers used on each device
			 */
			std::vector<const hardware::buffers::PRNGBuffer*> buffers;

			/**
			 * Reference to the system this PRNG is for.
			 */
			const hardware::System& system;
	};

	void gaussianComplexVector(hmc_complex * vector, int length, hmc_float sigma, physics::PRNG& prng);
	void gaussianNormalPair(hmc_float * z1, hmc_float * z2, physics::PRNG& prng);
	Matrixsu3 random_matrixsu3(physics::PRNG& prng);
}

#endif /* _PHYSICS_PRNG_ */
