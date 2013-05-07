/** @file
 * Implementation of the hardware::buffers::Buffer class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "prng_buffer.hpp"

#include "../device.hpp"
#include "../../meta/util.hpp"

hardware::buffers::PRNGBuffer::PRNGBuffer(size_t elems, Device * device)
	: Buffer(elems * sizeof(prng_state_t), device), elems(elems)
{
	// already inited
}

hardware::buffers::PRNGBuffer::PRNGBuffer(Device * device, const meta::Inputparameters& params)
	: PRNGBuffer(get_prng_buffer_size(device, params), device)
{
	// already inited
}

size_t hardware::buffers::PRNGBuffer::get_elements() const noexcept
{
	return elems;
}

size_t hardware::buffers::get_prng_buffer_size(const hardware::Device * device, const meta::Inputparameters& params)
{
	if(params.get_use_same_rnd_numbers()) {
		    return 1.;
	} else {
#ifdef USE_PRNG_NR3
		// Prepare random number arrays, for each task and device separately
		if(device->get_device_type() == CL_DEVICE_TYPE_GPU) {
			return 5120;
		} else {
			return 64;
		}
#elif defined(USE_PRNG_RANLUX)
		// make num of random states equal to default num of global threads
		// TODO make this somewhat more automatic (avoid code duplication)
		if(device->get_device_type() == CL_DEVICE_TYPE_GPU) {
			return 4 * device->get_preferred_local_thread_num() * device->get_num_compute_units();
		} else {
			return device->get_preferred_local_thread_num() * device->get_num_compute_units();
		}
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
	}
}

void hardware::buffers::PRNGBuffer::load(const prng_state_t * array) const
{
	Buffer::load(array);
}

void hardware::buffers::PRNGBuffer::dump(prng_state_t * array) const
{
	Buffer::dump(array);
}
