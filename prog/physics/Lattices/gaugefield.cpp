/** @file
 * Implementation of the physics::lattices::Gaugefield class
 */

#include "gaugefield.hpp"
#include "../../meta/util.hpp"
#include "../../logger.hpp"
//#include "../../hardware/code/gaugefield.hpp"
//#include "../../hardware/code/prng.hpp"

static std::vector<const hardware::buffers::SU3 *> allocate_buffers(hardware::System& system);

static void set_hot(std::vector<const hardware::buffers::SU3 *> buffers, std::vector<const hardware::buffers::PRNGBuffer *> prngs);
static void set_cold(std::vector<const hardware::buffers::SU3 *> buffers);

physics::lattices::Gaugefield::Gaugefield(hardware::System& system, physics::PRNG& prng, bool hot)
	: system(system), prng(prng), buffers(allocate_buffers(system))
{
	// TODO initialize as hot or cold
	if(hot) {
		set_hot(buffers, prng.buffers);
	} else {
		set_cold(buffers);
	}
}

static  std::vector<const hardware::buffers::SU3 *> allocate_buffers(hardware::System& system)
{
//	using hardware::buffers::SU3;
//
//	// only use device 0 for now
//	hardware::Device * device = system.get_devices().at(0);
//	std::vector<const SU3 *> buffers;
//	buffers.push_back(new SU3(meta::get_vol4d(system.get_inputparameters()) * 4, device));
//	return buffers;
}

static void set_hot(std::vector<const hardware::buffers::SU3 *> buffers, std::vector<const hardware::buffers::PRNGBuffer *> prngs)
{
//	using hardware::Device;
//	using namespace hardware::buffers;
//
//	for(size_t i = 0; i < buffers.size(); ++i) {
//		buffers[i]->get_device()->get_gaugefield().set_hot(buffers[i], prngs.at(i));
//	}
}

static void set_cold(std::vector<const hardware::buffers::SU3 *> buffers)
{
//	using hardware::Device;
//	using hardware::buffers::SU3;
//
//for(const SU3 * buffer: buffers) {
//		buffer->get_device()->get_gaugefield().set_cold(buffer);
//	}
}
