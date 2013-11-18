/*
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

__kernel void convertGaugefieldToSOA(__global Matrixsu3StorageType * const restrict out, __global const Matrixsu3 * const restrict in)
{
	// we need to take care of index converion. in the AOS storage the dimension is the continious index
	// in the soa storage we want the space indices to be continuous and have the dimension as outermost.
	for(uint d = 0; d < NDIM; ++d) {
		for(site_idx s = get_global_id(0); s < VOL4D_MEM; s += get_global_size(0)) {
			const st_idx site = get_st_idx_from_site_idx(s);
			Matrixsu3 tmp = in[get_link_idx_AOS(d, site)];

			putSU3(out, get_link_idx_SOA(d, site), tmp);
		}
	}
}
__kernel void convertGaugefieldFromSOA(__global Matrixsu3 * const restrict out, __global const Matrixsu3StorageType * const restrict in)
{
	// we need to take care of index converion. in the AOS storage the dimension is the continious index
	// in the soa storage we want the space indices to be continuous and have the dimension as outermost.
	for(uint d = 0; d < NDIM; ++d) {
		for(site_idx s = get_global_id(0); s < VOL4D_MEM; s += get_global_size(0)) {
			const st_idx site = get_st_idx_from_site_idx(s);
			Matrixsu3 tmp = getSU3(in, get_link_idx_SOA(d, site));

			out[get_link_idx_AOS(d, site)] = tmp;
		}
	}
}
