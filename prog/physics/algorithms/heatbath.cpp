/** @file
 * Implementation of the heatbath algorithm
 */

#include "heatbath.hpp"

#include <cassert>
#include "../hardware/code/heatbath.hpp"

void physics::algorithms::heatbath(physics::lattices::Gaugefield& gf, physics::PRNG& prng, int overrelax)
{
	assert(overrelax >= 0);
	assert(gf.get_buffers().size() == 1);

	// run heatbath
	auto gf_dev = gf.get_buffers()[0];
	auto prng_dev = prng.get_buffers()[0];
	auto code = gf_dev->get_device()->get_heatbath_code();

	code->run_heatbath(gf_dev, prng_dev);

	// add overrelaxation
	if(overrelax > 0) {
		physics::algorithms::overrelax(gf, prng, overrelax);
	}
}

void physics::algorithms::overrelax(physics::lattices::Gaugefield& gf, physics::PRNG& prng, int steps)
{
	assert(steps > 0);
	assert(gf.get_buffers().size() == 1);

	auto gf_dev = gf.get_buffers()[0];
	auto prng_dev = prng.get_buffers()[0];
	auto code = gf_dev->get_device()->get_heatbath_code();

	for(int i = 0; i < steps; ++i)
		code->run_overrelax(gf_dev, prng_dev);
}
