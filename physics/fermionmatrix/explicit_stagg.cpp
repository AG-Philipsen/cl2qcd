/** @file
 * Implementation of explicit fermionamtrix operations
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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

#include "fermionmatrix_stagg.hpp"


void physics::fermionmatrix::DKS_eo(const physics::lattices::Staggeredfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Staggeredfield_eo& in, int evenodd)
{
	auto out_bufs = out->get_buffers();
	auto gf_bufs = gf.get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != gf_bufs.size() || num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->getFermionStaggeredCode();
		fermion_code->D_KS_eo_device(in_bufs[i], out_bufs[i], gf_bufs[i], evenodd);
	}
	if(num_bufs!=1)
	  out->update_halo();
}

