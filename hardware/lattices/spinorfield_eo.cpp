/*
 * Copyright 2016 Francesca Cuteri
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

#include "spinorfield_eo.hpp"
#include "../code/spinors.hpp"
#include "../code/fermions.hpp"
#include "../buffers/halo_update.hpp"
#include <algorithm>

hardware::lattices::Spinorfield_eo::Spinorfield_eo(const hardware::System& system)
	: system(system), buffers(allocate_buffers())
#ifdef LAZY_HALO_UPDATES
	  , valid_halo_width(0)
#endif
{
}

std::vector<const hardware::buffers::Spinor *> hardware::lattices::Spinorfield_eo::allocate_buffers()
{
	using hardware::buffers::Spinor;

	auto devices = system.get_devices();
	std::vector<const Spinor*> buffers;
	buffers.reserve(devices.size());
	for(auto device: devices) {
		buffers.push_back(new Spinor(hardware::code::get_eoprec_spinorfieldsize(device->getLocalLatticeMemoryExtents()), device));
	}
	return buffers;
}

hardware::lattices::Spinorfield_eo::~Spinorfield_eo()
{
for(auto buffer: buffers) {
		delete buffer;
	}
}

const std::vector<const hardware::buffers::Spinor *> hardware::lattices::Spinorfield_eo::get_buffers() const noexcept
{
	return buffers;
}

void hardware::lattices::Spinorfield_eo::mark_halo_dirty() const
{
#ifdef LAZY_HALO_UPDATES
	logger.trace() << "Halo of Spinorfield_eo " << this << " marked as dirty.";
	valid_halo_width = 0;
#else
	update_halo();
#endif
}

void hardware::lattices::Spinorfield_eo::require_halo(unsigned reqd_width) const
{
#ifdef LAZY_HALO_UPDATES
	if(!reqd_width) {
		reqd_width = buffers[0]->get_device()->getHaloExtent();
	}
	logger.debug() << "Halo of Spinorfield_eo " << this << " required with width " << reqd_width << ". Width of valid halo: " << valid_halo_width;
	if(valid_halo_width < reqd_width) {
		update_halo(reqd_width);
		valid_halo_width = reqd_width;
	}
#endif
}

hardware::lattices::Spinorfield_eoHaloUpdate hardware::lattices::Spinorfield_eo::require_halo_async(unsigned reqd_width) const
{
#ifdef LAZY_HALO_UPDATES
	if(!reqd_width) {
		reqd_width = buffers[0]->get_device()->getHaloExtent();
	}
	logger.debug() << "Async Halo of Spinorfield_eo " << this << " required with width " << reqd_width << ". Width of valid halo: " << valid_halo_width;
	if(valid_halo_width < reqd_width) {
		return update_halo_async(reqd_width);
	} else {
		return Spinorfield_eoHaloUpdate(*this);
	}
#endif
}

void hardware::lattices::Spinorfield_eo::mark_halo_clean(unsigned width) const
{
#ifdef LAZY_HALO_UPDATES
	valid_halo_width = width ? width : buffers[0]->get_device()->getHaloExtent();
	logger.trace() << "Halo of Spinorfield_eo " << this << " marked as clean (width " << valid_halo_width << ").";
#endif
}

void hardware::lattices::Spinorfield_eo::update_halo(unsigned width) const
{
	logger.trace() << "Updating halo of Spinorfield_eo " << this;
	if(buffers.size() > 1) { // for a single device this will be a noop
		// currently either all or none of the buffers must be SOA
		if(buffers[0]->is_soa()) {
			update_halo_soa(width);
		} else {
			update_halo_aos();
		}
	}
}

hardware::lattices::Spinorfield_eoHaloUpdate hardware::lattices::Spinorfield_eo::update_halo_async(unsigned width) const
{
	logger.debug() << "Starting async update of halo of Spinorfield_eo " << this;
	if(buffers.size() > 1) { // for a single device this will be a noop
		// currently either all or none of the buffers must be SOA
		if(buffers[0]->is_soa()) {
			update_halo_soa_async(width);
		} else {
			update_halo_aos();
		}
	}
	return Spinorfield_eoHaloUpdate(*this, width);
}

void hardware::lattices::Spinorfield_eo::update_halo_soa(const unsigned width) const
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(!buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo_soa<spinor>(buffers, system, .5 /* only even or odd sites */, 1, width );
}

void hardware::lattices::Spinorfield_eo::update_halo_soa_async(const unsigned width) const
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(!buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::initialize_update_halo_soa<spinor>(buffers, system, .5 /* only even or odd sites */, 1, width );
}

void hardware::lattices::Spinorfield_eo::update_halo_aos() const
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo<spinor>(buffers, system, .5 /* only even or odd sites */ );
}

void hardware::lattices::Spinorfield_eoHaloUpdate::finalize()
{
	logger.trace() << "Finalizing update on Spinorfield_eo " << &target;
#ifdef LAZY_HALO_UPDATES
	// 0 is used to indicate that the update is already complete (or was not required)
	if(reqd_halo_width) {
		target.update_halo_finalize(reqd_halo_width);
		target.valid_halo_width = reqd_halo_width;
		reqd_halo_width = 0; // mark self as done to avoid duplicate call
	}
#endif
}

void hardware::lattices::Spinorfield_eo::update_halo_finalize(unsigned width) const
{
	logger.debug() << "Finalizing async update of halo of Spinorfield_eo " << this;
	if(buffers.size() > 1) { // for a single device this will be a noop
		// currently either all or none of the buffers must be SOA
		if(buffers[0]->is_soa()) {
			logger.trace() << "This is SoA";
			update_halo_soa_finalize(width);
		}
		// NOOP for AoS
	}
}

void hardware::lattices::Spinorfield_eo::update_halo_soa_finalize(const unsigned width) const
{
	hardware::buffers::finalize_update_halo_soa<spinor>(buffers, system, .5 /* only even or odd sites */, 1, width );
}
