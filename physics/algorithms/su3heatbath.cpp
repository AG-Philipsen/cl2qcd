/** @file
 * Implementation of the su3heatbath algorithm
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

#include "su3heatbath.hpp"

#include <cassert>
#include "../../hardware/device.hpp"
#include "../../hardware/code/heatbath.hpp"

void physics::algorithms::su3heatbath(physics::lattices::Gaugefield& gf, physics::PRNG& prng, int overrelax)
{
	assert(overrelax >= 0);
	assert(gf.get_buffers().size() == 1);

	// run su3heatbath
	auto gf_dev = gf.get_buffers()[0];
	auto prng_dev = prng.get_buffers()[0];
	auto code = gf_dev->get_device()->getHeatbathCode();

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
	auto code = gf_dev->get_device()->getHeatbathCode();

	for(int i = 0; i < steps; ++i)
		code->run_overrelax(gf_dev, prng_dev);
}
