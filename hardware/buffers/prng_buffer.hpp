/** @file
 * Declaration of the hardware::buffers::Buffer class
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#ifndef _HARDWARE_BUFFERS_PRN_GBUFFER_
#define _HARDWARE_BUFFERS_PRN_GBUFFER_

#include "buffer.hpp"
#include "../device.hpp"

namespace hardware {
namespace buffers {

/**
 * Get the default size PRNG buffer on the given device.
 *
 * \param device The device the buffer is to be used with
 * \return The size of the buffer in element
 * @Todo: move elsewhere
 */
size_t get_prng_buffer_size(const Device * device, const bool useSameRandomNumbers);

/**
 * A PRNG OpenCL buffer.
 */
class PRNGBuffer : public Buffer {

public:
#ifdef USE_PRNG_RANLUX
	typedef cl_float4 prng_state_t[7];
	static_assert(sizeof(prng_state_t) == 7 * sizeof(cl_float4), "PRNG state type mockup is of wrong size.");
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX

	/**
	 * Allocate a buffer with the default number of
	 * elements for this device.
	 *
	 * \param device The device to locate the buffer on
	 * \param useSameRandomNumbers This influences the size of the buffer
	 */
	PRNGBuffer(Device * device, const bool useSameRandomNumbers);

	/*
	 * PRNGBuffers cannot be copied
	 */
	PRNGBuffer& operator=(const PRNGBuffer&) = delete;
	PRNGBuffer(const PRNGBuffer&) = delete;
	PRNGBuffer() = delete;

	/**
	 * Get the buffer size in bytes
	 */
	size_t get_elements() const noexcept;

	void load(const prng_state_t *) const;

	void dump(prng_state_t *) const;

private:
	/**
	 * Allocate a buffer of the given size in elements.
	 *
	 * \param elems The size of the buffer in elements
	 * \param device The device to locate the buffer on
	 */
	PRNGBuffer(size_t elems, Device * device);

	/**
	 * The size of the buffer in elements.
	 */
	const size_t elems;
};
}
}

#endif /* _HARDWARE_BUFFERS_PRNG_BUFFER_ */

