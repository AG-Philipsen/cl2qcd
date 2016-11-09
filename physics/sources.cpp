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
#include "../hardware/code/correlator.hpp"
#include "../hardware/code/correlator_staggered.hpp"
#include "../geometry/index.hpp"

void physics::set_point_source(const physics::lattices::Spinorfield * spinorfield, int k, const physics::SourcesParametersInterface& params)
{
	if(k >= 12 || k < 0) {
		throw std::invalid_argument("k must be within 0..11");
	}

	spinorfield->zero();

	auto buffers = spinorfield->get_buffers();

	// only execute on the buffer where the given position results.
	int t_pos = params.getSourceT();
	unsigned local_lattice_size = buffers[0]->get_device()->getLocalLatticeExtents().tExtent;
	auto buffer = buffers[t_pos / local_lattice_size];
	auto device = buffer->get_device();
	int local_t = t_pos % local_lattice_size;

	device->getCorrelatorCode()->create_point_source_device(buffer, k, Index(params.getSourceX(), params.getSourceY(), params.getSourceZ(), params.getSourceT(), LatticeExtents(params.getNs(), params.getNt())).spatialIndex, local_t);

	spinorfield->update_halo();
}

void physics::set_volume_source(const physics::lattices::Spinorfield * spinorfield, const PRNG& prng)
{
	spinorfield->zero();

	auto buffers = spinorfield->get_buffers();

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto buffer = buffers[i];
		auto prng_buffer = prng.get_buffers().at(i);

		buffer->get_device()->getCorrelatorCode()->create_volume_source_device(buffer, prng_buffer);
	}

	spinorfield->update_halo();
}

void physics::set_timeslice_source(const physics::lattices::Spinorfield * spinorfield, const PRNG& prng, int t_pos)
{
	spinorfield->zero();

	auto buffers = spinorfield->get_buffers();

	unsigned local_lattice_size = buffers[0]->get_device()->getLocalLatticeExtents().tExtent;
	unsigned t_buf = t_pos / local_lattice_size;
	auto buffer = buffers[t_buf];
	auto device = buffer->get_device();
	int local_t = t_pos % local_lattice_size;
	auto prng_buffer = prng.get_buffers().at(t_buf);

	device->getCorrelatorCode()->create_timeslice_source_device(buffer, prng_buffer, local_t);

	spinorfield->update_halo();
}

void physics::set_zslice_source(const physics::lattices::Spinorfield * spinorfield, const PRNG& prng, int z)
{
	spinorfield->zero();

	auto buffers = spinorfield->get_buffers();

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto buffer = buffers[i];
		auto prng_buffer = prng.get_buffers().at(i);

		buffer->get_device()->getCorrelatorCode()->create_zslice_source_device(buffer, prng_buffer, z);
	}

	spinorfield->update_halo();
}


static void fill_sources(const std::vector<physics::lattices::Spinorfield *>& sources, const physics::PRNG& prng, const physics::SourcesParametersInterface& params);

std::vector<physics::lattices::Spinorfield *> physics::create_sources(const hardware::System& system, const PRNG& prng, const size_t n_sources, physics::InterfacesHandler & interfacesHandler)
{
    const physics::SourcesParametersInterface & sourcesParameters = interfacesHandler.getSourcesParametersInterface();

    auto sources = lattices::create_spinorfields(system, n_sources, interfacesHandler, sourcesParameters.placeSourcesOnHost());
    fill_sources(sources, prng, sourcesParameters);
    return sources;
}

static void fill_sources(const std::vector<physics::lattices::Spinorfield *>& sources, const physics::PRNG& prng, const physics::SourcesParametersInterface& params)
{
    using namespace physics;

    for(size_t k = 0; k < sources.size(); k++) {
        auto source = sources[k];
        try_swap_in(source);

        switch(params.getSourceType()) {
            case common::point:
                logger.debug() << "start creating point-source...";
                //CP: k has to be between 0..12 (corresponds to spin-color index)
                set_point_source(source, (k % 12), params);
                break;
            case common::sourcetypes::volume:
                logger.debug() << "start creating volume-source...";
                set_volume_source(source, prng);
                break;
            case common::sourcetypes::timeslice:
                logger.debug() << "start creating timeslice-source...";
                set_timeslice_source(source, prng, params.getSourceT());
                break;
            case common::sourcetypes::zslice:
                logger.debug() << "start creating zslice-source...";
                set_zslice_source(source, prng, params.getSourceZ());
                break;
            default:
                throw std::domain_error("no such sourcetype");
        }

        try_swap_out(source);
    }
}

std::vector<physics::lattices::Spinorfield *> physics::create_swappable_sources(const hardware::System& system, const PRNG& prng, const size_t n_sources, physics::InterfacesHandler & interfacesHandler)
{
    const physics::SourcesParametersInterface & sourcesParameters = interfacesHandler.getSourcesParametersInterface();

    auto sources = lattices::create_swappable_spinorfields(system, n_sources, interfacesHandler, sourcesParameters.placeSourcesOnHost());
    fill_sources(sources, prng, sourcesParameters);
    return sources;
}

/*******************************************************************************************************************************/

void physics::set_volume_source(const physics::lattices::Staggeredfield_eo * inout, const PRNG& prng)
{
	auto buffers = inout->get_buffers();

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto buffer = buffers[i];
		auto prng_buffer = prng.get_buffers().at(i);

		buffer->get_device()->getCorrelatorStaggeredCode()->create_volume_source_stagg_eoprec_device(buffer, prng_buffer);
	}

	if(buffers.size()!=1)
	  inout->update_halo();
}

void physics::set_point_source(const physics::lattices::Staggeredfield_eo * source, int k, const physics::SourcesParametersInterface& params)
{
	if(k > 2 || k < 0) {
		throw std::invalid_argument("k must be within 0..2");
	}

	auto buffers = source->get_buffers();

	// only execute on the buffer where the given position results (with single device it is trivial but it works!).
	int t_pos = params.getSourceT();
	unsigned local_lattice_size = buffers[0]->get_device()->getLocalLatticeExtents().tExtent;
	auto buffer = buffers[t_pos / local_lattice_size];
	auto device = buffer->get_device();
	int local_t = t_pos % local_lattice_size;

	device->getCorrelatorStaggeredCode()->create_point_source_stagg_eoprec_device(buffer, k, Index(params.getSourceX(), params.getSourceY(), params.getSourceZ(), params.getSourceT(), LatticeExtents(params.getNs(), params.getNt())).spatialIndex, local_t);

	if(buffers.size() != 1){
	    throw Print_Error_Message("Point source for staggered fermions not implemented for multiple devices!", __FILE__, __LINE__);
	    source->update_halo();
	}

}


static void fill_sources(const std::vector<physics::lattices::Staggeredfield_eo *>& sources, const physics::PRNG& prng, const physics::SourcesParametersInterface& params)
{
    using namespace physics;

    for(size_t k = 0; k < sources.size(); k++) {
        auto source = sources[k];

        switch(params.getSourceType()) {
            case common::point:
                logger.debug() << "Creating point-source...";
                //k has to be between 0..2 (corresponds to color index, no spin in staggered)
                set_point_source(source, (k % 3), params);
                break;
            case common::sourcetypes::volume:
                logger.debug() << "Creating volume-source...";
                set_volume_source(source, prng);
                break;
            default:
                throw std::domain_error("no such sourcetype");
        }

        logger.debug() << "...source created!";
    }
}

std::vector<physics::lattices::Staggeredfield_eo *> physics::create_staggered_sources(const hardware::System& system, const PRNG& prng,
                                                                                      const size_t numberOfSources, physics::InterfacesHandler & interfacesHandler)
{
    const physics::SourcesParametersInterface & sourcesParameters = interfacesHandler.getSourcesParametersInterface();

    std::vector<physics::lattices::Staggeredfield_eo *> sources;
    sources.reserve(numberOfSources);

    for(size_t i = 0; i < numberOfSources; ++i) {
        sources.push_back(new physics::lattices::Staggeredfield_eo(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>()));
    }

    fill_sources(sources, prng, sourcesParameters);

    return sources;
}

