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

// real part of complex (!!!) scalarproduct, return in result
// --> use 2 kernels: 1 for the summation in one block and 1 for summation over blockresults
/// NOTE: The reduction used in this kernel is only safe with ls being a power of 2 and bigger than 8!

// Description of variables of scalar_product_staggered_eoprec kernel:
//  - x: The first staggered field (an su3vec per each site => vector of VOL4D/2
//       components that are su3vec varibles)
//  - y: The second staggered field (an su3vec per each site => vector of VOL4D/2
//       components that are su3vec varibles)
//  - result: Vector of hmc_float that will contain the sums of the components
//            of the result_local vectors. In other words, each component of 
//            this vector will contain the sum of all scalarproduct that have
//            been mapped to the threads within a group. Therefore result is a
//            vector with num_groups components.
//  - result_local: Vector with local_size components. At the end of the scalar_prduct
//                  kernel, its first component will be the sum of all its
//                  components (and will be put in result[group_id]).
//                  Observe that some components of result_local can include
//                  the sum of squarenorms of several fields (if EOPREC_SPINORFIELDSIZE_LOCAL>global_size).

__kernel void scalar_product_real_part_staggered_eoprec( __global const staggeredStorageType  * const x, __global const staggeredStorageType * const y, __global hmc_float * const result, __local hmc_float * const result_local )
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int group_id = get_group_id (0);
	int idx = get_local_id(0);

	hmc_float sum = 0.0;
	
	for(int id_local = id; id_local < EOPREC_SPINORFIELDSIZE_LOCAL; id_local += global_size) {
		site_idx id_mem = get_eo_site_idx_from_st_idx(get_even_st_idx_local(id_local));
		su3vec x_tmp = get_su3vec_from_field_eo(x, id_mem);
		su3vec y_tmp = get_su3vec_from_field_eo(y, id_mem);
		hmc_float tmp = su3vec_scalarproduct_real_part(x_tmp, y_tmp);
		sum += tmp;
	}

	if(local_size == 1) {
		result[ group_id ] = sum;
	} else {
	  	result_local[idx] = sum;
		// sync threads
		barrier(CLK_LOCAL_MEM_FENCE);

		//reduction until threads 0-7 hold all partial sums
		int cut1;
		int cut2 = local_size;
		for(cut1 = local_size / 2; cut1 > 4; cut1 /= 2) {
		  for(int i = idx + cut1; i < cut2; i += cut1) {
		    result_local[idx] += result_local[i];
		  }
		  barrier(CLK_LOCAL_MEM_FENCE);
		  cut2 = cut1;
		}
		//thread 0 sums up the last 8 results and stores them in the global buffer
		if (idx == 0) {
		  result[ group_id ] =  result_local[0] + result_local[1] + result_local[2] + result_local[3] +
		                        result_local[4] + result_local[5] + result_local[6] + result_local[7];
		}
	}
}
