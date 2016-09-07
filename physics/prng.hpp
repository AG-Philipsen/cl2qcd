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
#include "../hardware/lattices/prng.hpp"
#include "prngInterface.hpp"

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
			PRNG(const hardware::System& system, const physics::PrngParametersInterface * parametersIn);

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
			double get_double() const noexcept;

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

			/**
			 * Save current PRNG state to a file with name prng_prefix + 'save' + prng_postfix
			 * @param[in] number The trajectory number to be stored in the file
			 */
			void save();
			/**
			 * Save current PRNG state to a file with specific name based on
			 * current trajectory number.
			 * @param[in] number The trajectory number to be stored in the file
			 */
			void saveToSpecificFile(int number);

			void verifyWritingWasSuccessful(const std::string filename) const;

			bool operator == (const physics::PRNG & prng) const;
			bool operator != (const physics::PRNG & prng) const;

			std::string getName(int = -1) const noexcept;
		private:

			/**
			 * Reference to the system this PRNG is for.
			 */
			const hardware::System& system;
			const physics::PrngParametersInterface * parameters;
			hardware::lattices::PRNG const hPrng;
	};
}

#endif /* _PHYSICS_PRNG_ */
