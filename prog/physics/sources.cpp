/** @file
 * Implementation of functions working with sources.
 */

#include "sources.hpp"
#include <cassert>
#include <stdexcept>
#include "../host_geometry.h"

void physics::set_point_source(const physics::lattices::Spinorfield * spinorfield, int k, const meta::Inputparameters& params)
{
	if(k >= 12 || k < 0) {
		throw std::invalid_argument("k must be within 0..11");
	}

	auto buffers = spinorfield->get_buffers();

	// for now require a single device
	assert(buffers.size() == 1);

	auto buffer = buffers[0];
	auto device = buffer->get_device();

	device->get_correlator_code()->create_point_source_device(buffer, k, get_source_pos_spatial(params), params.get_source_t());
}

void physics::set_volume_source(const physics::lattices::Spinorfield * spinorfield, PRNG& prng)
{
	auto buffers = spinorfield->get_buffers();

	// for now require a single device
	assert(buffers.size() == 1);

	auto buffer = buffers[0];
	auto prng_buffer = prng.get_buffers().at(0);

	buffer->get_device()->get_correlator_code()->create_volume_source_device(buffer, prng_buffer);
}

void physics::set_timeslice_source(const physics::lattices::Spinorfield * spinorfield, PRNG& prng, int t)
{
	auto buffers = spinorfield->get_buffers();

	// for now require a single device
	assert(buffers.size() == 1);

	auto buffer = buffers[0];
	auto prng_buffer = prng.get_buffers().at(0);

	buffer->get_device()->get_correlator_code()->create_timeslice_source_device(buffer, prng_buffer, t);
}

void physics::set_zslice_source(const physics::lattices::Spinorfield * spinorfield, PRNG& prng, int z)
{
	auto buffers = spinorfield->get_buffers();

	// for now require a single device
	assert(buffers.size() == 1);

	auto buffer = buffers[0];
	auto prng_buffer = prng.get_buffers().at(0);

	buffer->get_device()->get_correlator_code()->create_zslice_source_device(buffer, prng_buffer, z);
}

const std::vector<const physics::lattices::Spinorfield *> physics::create_sources(hardware::System& system, PRNG& prng, const size_t n_sources)
{
	auto params = system.get_inputparameters();
	//	const size_t n_sources = params.get_num_sources();

	const std::vector<const physics::lattices::Spinorfield *> sources = lattices::create_spinorfields(system, n_sources);

	for(int k = 0; k < n_sources; k++) {
		auto source = sources[k];

		switch(system.get_inputparameters().get_sourcetype()) {
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
	}
	return sources;
}
