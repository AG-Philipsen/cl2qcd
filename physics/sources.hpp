/** @file
 * Functions for working with sources
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

#ifndef _PHYSICS_SOURCES_
#define _PHYSICS_SOURCES_

#include "../hardware/system.hpp"
#include "prng.hpp"
#include "lattices/spinorfield.hpp"
#include "lattices/swappable_spinorfield.hpp"
#include "lattices/staggeredfield_eo.hpp"

namespace physics {

	/**
	 * Create sources as specified by the input parameters of the system.
	 * TODO: Turn these two function into a template and figure out how to treat the
	 *       place_on_host parameter that the Spinorfield has and that Staggeredfield_eo
	 *       has not for the moment!
	 */
	std::vector<lattices::Spinorfield *> create_sources(const hardware::System& system, const PRNG& prng, const size_t n_sources,
	                                                    physics::InterfacesHandler & interfacesHandler);
	std::vector<physics::lattices::Staggeredfield_eo *> create_staggered_sources(const hardware::System& system, const PRNG& prng,
	                                                                             const size_t numberOfSources, physics::InterfacesHandler & interfacesHandler);

	/**
	 * Create a set of spinorfields that can be swapped.
	 * Return normal spinorfield pointers for compatibility reasons.
	 */
	std::vector<lattices::Spinorfield *> create_swappable_sources(const hardware::System& system, const PRNG& prng, const size_t n_sources,
	                                                              physics::InterfacesHandler & interfacesHandler);

	void set_point_source(const physics::lattices::Spinorfield *, int k, const physics::SourcesParametersInterface& params);
	void set_volume_source(const physics::lattices::Spinorfield *, const PRNG& prng);
	void set_timeslice_source(const physics::lattices::Spinorfield *, const PRNG& prng, int t);
	void set_zslice_source(const physics::lattices::Spinorfield *, const PRNG& prng, int z);
	//Staggered sources
	void set_volume_source(const physics::lattices::Staggeredfield_eo *, const PRNG& prng);
	void set_point_source(const physics::lattices::Staggeredfield_eo *, int k, const physics::SourcesParametersInterface& params);
}

#endif /* _PHYSICS_SOURCES_ */
