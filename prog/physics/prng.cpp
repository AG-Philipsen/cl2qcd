/** @file
 * Ranlux PRNG implementation
 */

#include "prng.hpp"

#include "../host_random.h"
#include "../hardware/buffers/prng_buffer.hpp"
//#include "../hardware/code/prng.hpp"

physics::PRNG::~PRNG()
{
for(const hardware::buffers::PRNGBuffer * buffer : buffers) {
		delete buffer;
	}
}

physics::PRNG::PRNG(const hardware::System& system) :
	system(system)
{
	using hardware::buffers::PRNGBuffer;

	// initialize host prng
	uint32_t seed = system.get_inputparameters().get_host_seed();
	prng_init(seed);

	// initialize devices
for(hardware::Device * device : system.get_devices()) {
		// create a buffer for each device
		const PRNGBuffer * buffer = new PRNGBuffer(device);
		auto code = device->get_prng_code();
		code->initialize(buffer, ++seed);
		buffers.push_back(buffer);
	}
}
