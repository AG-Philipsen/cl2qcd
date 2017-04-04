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

#ifndef _HARDWARE_LATTICES_SPINORFIELD_EO_
#define _HARDWARE_LATTICES_SPINORFIELD_EO_

#include "../system.hpp"
#include "../buffers/spinor.hpp"
#include "spinorfield.hpp"
#include "scalar.hpp"
#include "../../common_header_files/types_fermions.h"


namespace hardware {

namespace lattices {

class Spinorfield_eoHaloUpdate;

class Spinorfield_eo
{
	friend Spinorfield_eoHaloUpdate;

public:

	Spinorfield_eo(const hardware::System&);

	Spinorfield_eo& operator=(const Spinorfield_eo&) = delete;
	Spinorfield_eo(const Spinorfield_eo&) = delete;
	Spinorfield_eo() = delete;

	virtual ~Spinorfield_eo();

	std::vector<const hardware::buffers::Spinor *> allocate_buffers();
	const std::vector<const hardware::buffers::Spinor *> get_buffers() const noexcept;

	void mark_halo_dirty() const;

	void require_halo(unsigned width = 0) const;

	Spinorfield_eoHaloUpdate require_halo_async(unsigned width = 0) const;

	void mark_halo_clean(unsigned width = 0) const;

private:
	hardware::System const& system;
	const std::vector<const hardware::buffers::Spinor *> buffers;
	/**
	 * Unconditionally update the halo.
	 *
	 * \param widh Up to which thickness to update the halo. Use 0 to indicate the full halo shall be updated.
	 */
	void update_halo(unsigned width = 0) const;
	Spinorfield_eoHaloUpdate update_halo_async(unsigned width = 0) const;
	void update_halo_finalize(unsigned width = 0) const;
	void update_halo_soa(const unsigned width) const;
	void update_halo_soa_async(const unsigned width) const;
	void update_halo_soa_finalize(const unsigned width) const;
	void update_halo_aos() const;
#ifdef LAZY_HALO_UPDATES
	mutable unsigned valid_halo_width;
#endif
};

class Spinorfield_eoHaloUpdate {
		friend Spinorfield_eo;
	public:
		/**
		 * Complete the halo update.
		 *
		 * Access the the halo data happens synchroneous to the default queue of each buffer.
		 * Therefore commands using the default queue of this buffer can safely operate the halo elements if queued after this call.
		 *
		 * Note that depending on the exact transfer methods this might cause data transfer between devices.
		 */
		void finalize();
	private:
		/**
		 * Construct the handler object. Only done by Spinorfield_eo.
		 *
		 * \param target The Spinorfield_eo on which the update is performed.
		 * \param reqd_halo_width Width of the updated halo segment.
		 *        0 indicates that update is finished / has alredy been completed and makes finish() a NOOP.
		 */
		Spinorfield_eoHaloUpdate(Spinorfield_eo const & target, unsigned const & reqd_halo_width = 0)
		 : target(target), reqd_halo_width(reqd_halo_width) {  };

		Spinorfield_eo const & target;
		unsigned reqd_halo_width;
};

}

}
#endif /* _HARDWARE_LATTICES_SPINORFIELD_EO_ */
