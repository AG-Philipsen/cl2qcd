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

/** @todo add args for reduction... */
/// NOTE: The reduction used in this kernel is only safe with ls being a power of 2 and bigger than 8!
__kernel void gaugemomentum_squarenorm(__global const aeStorageType * const restrict in, __global hmc_float * const restrict out, __local hmc_float * const restrict result_local)
{
	hmc_float sum = 0.f;
	PARALLEL_FOR(id_local, VOL4D_LOCAL) {
		site_idx site = get_site_idx((id_local % 2 == 0) ? get_even_st_idx_local(id_local / 2) : get_odd_st_idx_local(id_local / 2));
		for(uint d = 0; d < 4; ++d) {
			const link_idx id_mem = get_link_idx(d, get_st_idx_from_site_idx(site));
			sum += ae_squarenorm(getAe(in, id_mem));
		}
	}

	const size_t group_id = get_group_id(0);
	const size_t idx = get_local_id(0);
	const size_t local_size = get_local_size(0);

	if(local_size == 1) {
		out[ group_id ] = sum;
	} else {
		result_local[idx] = sum;
		// sync threads
		barrier(CLK_LOCAL_MEM_FENCE);

		//reduction until threads 0-7 hold all partial sums
		int cut1;
		int cut2 = local_size;
		for(cut1 = local_size / 2; cut1 > 4; cut1 /= 2) {
		  for(int i = idx + cut1; i < cut2; i += cut1) {
		    result_local[idx] +=  result_local[i];
		  }
		  barrier(CLK_LOCAL_MEM_FENCE);
		  cut2 = cut1;
		}
		//thread 0 sums up the last 8 results and stores them in the global buffer
		if (idx == 0) {
		  out[ group_id ] =  result_local[0] + result_local[1] + result_local[2] + result_local[3] +
		    result_local[4] + result_local[5] + result_local[6] + result_local[7];
		}
	}
}
