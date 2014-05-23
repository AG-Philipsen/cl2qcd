/** @file
 * physics::observables::wilson::TwoFlavourCorrelators class
 *
 * Copyright 2014 Christopher Pinke
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

#include "wilsonTwoFlavourCorrelators.hpp"

#include "../physics/sources.hpp"
#include "../physics/algorithms/flavour_doublet.hpp"
#include "../physics/algorithms/inversion.hpp"

void physics::observables::wilson::measureTwoFlavourDoubletCorrelatorsOnGaugefield(const physics::lattices::Gaugefield * gaugefield, std::string currentConfigurationName) 
{
	auto parameters = gaugefield->getParameters();
	auto system = gaugefield->getSystem();
	auto prng = gaugefield->getPrng();
	
	std::string filenameForCorrelatorData = meta::get_ferm_obs_corr_file_name(*parameters, currentConfigurationName);
	// for the correlator calculation, all sources are needed on the device
	const std::vector<physics::lattices::Spinorfield*> sources = physics::create_swappable_sources(*system, *prng, parameters->get_num_sources());
	const std::vector<physics::lattices::Spinorfield*> result = physics::lattices::create_swappable_spinorfields(*system, sources.size(), parameters->get_place_sources_on_host());
	swap_out(sources);
	swap_out(result);
	physics::algorithms::perform_inversion(&result, gaugefield, sources, *system);
	logger.info() << "Finished inversion. Starting measurements.";
	swap_in(sources);
	swap_in(result);
	physics::algorithms::flavour_doublet_correlators(result, sources, filenameForCorrelatorData, *system);
	release_spinorfields(result);
	release_spinorfields(sources);
}