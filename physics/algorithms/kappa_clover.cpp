/** @file
 * Implementation of the heatbath algorithm
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

#include "kappa_clover.hpp"

#include <cassert>
#include "../../hardware/device.hpp"
#include "../../hardware/code/kappa.hpp"

hmc_float physics::algorithms::kappa_clover(physics::lattices::Gaugefield& gf, hmc_float beta)
{
	assert(gf.get_buffers().size() == 1);

	auto gf_dev = gf.get_buffers()[0];

	auto device = gf_dev->get_device();
	hardware::buffers::Plain<hmc_float> kappa_clover_dev(1, device);
	gf_dev->get_device()->get_kappa_code()->run_kappa_clover(&kappa_clover_dev, gf_dev, beta);

	hmc_float kappa_clover_host;
	kappa_clover_dev.dump(&kappa_clover_host);
	return kappa_clover_host;
}
