/** @file
 * Implementation of functions working with sources.
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

#include "sources.hpp"
#include <cassert>
#include <stdexcept>
#include "../host_functionality/host_geometry.h"
#include "../hardware/code/correlator.hpp"
#include "../hardware/code/correlator_staggered.hpp"

void physics::set_point_source(const physics::lattices::Spinorfield * spinorfield, int k, const meta::Inputparameters& params)
{
	if(k >= 12 || k < 0) {
		throw std::invalid_argument("k must be within 0..11");
	}

	spinorfield->zero();

	auto buffers = spinorfield->get_buffers();

	// only execute on the buffer were the given position results.
	int t_pos = params.get_source_t();
	unsigned local_lattice_size = buffers[0]->get_device()->get_local_lattice_size().t;
	auto buffer = buffers[t_pos / local_lattice_size];
	auto device = buffer->get_device();
	int local_t = t_pos % local_lattice_size;

	device->get_correlator_code()->create_point_source_device(buffer, k, get_source_pos_spatial(params), local_t);

	spinorfield->update_halo();
}

void physics::set_volume_source(const physics::lattices::Spinorfield * spinorfield, PRNG& prng)
{
	spinorfield->zero();

	auto buffers = spinorfield->get_buffers();

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto buffer = buffers[i];
		auto prng_buffer = prng.get_buffers().at(i);

		buffer->get_device()->get_correlator_code()->create_volume_source_device(buffer, prng_buffer);
	}

	spinorfield->update_halo();
}

void physics::set_timeslice_source(const physics::lattices::Spinorfield * spinorfield, PRNG& prng, int t_pos)
{
	spinorfield->zero();

	auto buffers = spinorfield->get_buffers();

	unsigned local_lattice_size = buffers[0]->get_device()->get_local_lattice_size().t;
	unsigned t_buf = t_pos / local_lattice_size;
	auto buffer = buffers[t_buf];
	auto device = buffer->get_device();
	int local_t = t_pos % local_lattice_size;
	auto prng_buffer = prng.get_buffers().at(t_buf);

	device->get_correlator_code()->create_timeslice_source_device(buffer, prng_buffer, local_t);

	spinorfield->update_halo();
}

void physics::set_zslice_source(const physics::lattices::Spinorfield * spinorfield, PRNG& prng, int z)
{
	spinorfield->zero();

	auto buffers = spinorfield->get_buffers();

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto buffer = buffers[i];
		auto prng_buffer = prng.get_buffers().at(i);

		buffer->get_device()->get_correlator_code()->create_zslice_source_device(buffer, prng_buffer, z);
	}

	spinorfield->update_halo();
}

//Steggered source
void physics::set_volume_source(const physics::lattices::Staggeredfield_eo * inout, PRNG& prng)
{
	auto buffers = inout->get_buffers();

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto buffer = buffers[i];
		auto prng_buffer = prng.get_buffers().at(i);

		buffer->get_device()->get_correlator_staggered_code()->create_volume_source_stagg_eoprec_device(buffer, prng_buffer);
	}

	if(buffers.size()!=1)
	  inout->update_halo();
}

static void fill_sources(const std::vector<physics::lattices::Spinorfield *>& sources, physics::PRNG& prng, const meta::Inputparameters& params);

std::vector<physics::lattices::Spinorfield *> physics::create_sources(hardware::System& system, PRNG& prng, const size_t n_sources)
{
	auto params = system.get_inputparameters();

	auto sources = lattices::create_spinorfields(system, n_sources, params.get_place_sources_on_host());
	fill_sources(sources, prng, params);
	return sources;
}

static void fill_sources(const std::vector<physics::lattices::Spinorfield *>& sources, physics::PRNG& prng, const meta::Inputparameters& params)
{
	using namespace physics;

	for(size_t k = 0; k < sources.size(); k++) {
		auto source = sources[k];
		try_swap_in(source);

		switch(params.get_sourcetype()) {
			case meta::Inputparameters::point:
				logger.debug() << "start creating point-source...";
				//CP: k has to be between 0..12 (corresponds to spin-color index)
				set_point_source(source, (k % 12), params);
				break;
			case meta::Inputparameters::sourcetypes::volume:
				logger.debug() << "start creating volume-source...";
				set_volume_source(source, prng);
				break;
			case meta::Inputparameters::sourcetypes::timeslice:
				logger.debug() << "start creating timeslice-source...";
				set_timeslice_source(source, prng, params.get_source_t());
				break;
			case meta::Inputparameters::sourcetypes::zslice:
				logger.debug() << "start creating zslice-source...";
				set_zslice_source(source, prng, params.get_source_z());
				break;
			default:
				throw std::domain_error("no such sourcetype");
		}

		try_swap_out(source);
	}
}

std::vector<physics::lattices::Spinorfield *> physics::create_swappable_sources(hardware::System& system, PRNG& prng, const size_t n_sources)
{
	auto params = system.get_inputparameters();

	auto sources = lattices::create_swappable_spinorfields(system, n_sources, params.get_place_sources_on_host());
	fill_sources(sources, prng, params);
	return sources;
}
