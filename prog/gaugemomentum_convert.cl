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

__kernel void gaugemomentum_convert_to_soa(__global aeStorageType * const restrict dest, __global const ae * const restrict src)
{
	PARALLEL_FOR(s, VOL4D_MEM) {
		for(uint d = 0; d < 4; ++d) {
			const st_idx site = get_st_idx_from_site_idx(s);
			const link_idx aos_i = get_link_idx_AOS(d, site);
			const link_idx soa_i = get_link_idx(d, site);
			putAe(dest, soa_i, src[aos_i]);
		}
	}
}

__kernel void gaugemomentum_convert_from_soa(__global ae * const restrict dest, __global const aeStorageType * const restrict src)
{
	PARALLEL_FOR(s, VOL4D_MEM) {
		for(uint d = 0; d < 4; ++d) {
			const st_idx site = get_st_idx_from_site_idx(s);
			const link_idx aos_i = get_link_idx_AOS(d, site);
			const link_idx soa_i = get_link_idx(d, site);
			dest[aos_i] = getAe(src, soa_i);
		}
	}
}
