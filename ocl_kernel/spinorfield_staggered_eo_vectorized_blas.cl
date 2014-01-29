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

// This kernels performs the sax operation with a set of alpha constants and then
// it evaluates the squarenorm of each resulting field.

// Description of variables:
//  - x: The input staggered field (an su3vec per each site => vector of VOL4D/2
//       components that are su3vec varibles)
//  - alpha: vector of real numbers
//  - num_fields: number of constants alpha, i.e. number of output fields
//  - result: Vector of hmc_float that will contain the sums of the components of the result_local vectors.
//            Note that it has num_groups*num_fields components.
//  - result_local: Vector with local_size components. At the end of the local reduction, its first num_fields
//                  components will be the sum of all its components.

__kernel void sax_vectorized_and_squarenorm_eoprec(__global const staggeredStorageType * const x, __global const hmc_float * alpha, const int num_fields, __global hmc_float * const restrict result, __local hmc_float * const restrict result_local)
{
	const int id = get_global_id(0);
	const int local_size = get_local_size(0);
	const int global_size = get_global_size(0);
	const int group_id = get_group_id(0);
	const int idx = get_local_id(0);

// 	if(id == 0){
// 	  for(uint i=0; i<num_fields; i++)
// 	    sqnorms[i]=0.0;
// 	}
// 	//sync threads
// 	barrier(CLK_LOCAL_MEM_FENCE);
	
	/* Here we have to give up having a completely general kernel, since we cannot allocate an array in
	 * the private memory with num_fields components (num_fields is not known at compilation time,
	 * but only at run time). Then we give this number as option to the kernel, but then we must
	 * chose a meaningful number. Since this function will be used to check the residuum in the CGM
	 * algorithm, then here we put the maximum number of equations that can be in the multimass inverter,
	 * namely the highest rational approximation order used in the RHMC.
	 * >> TODO Maybe it is worth renaming this kernel that is not so general any more.
	 */
	hmc_float sum[RA_MAX_ORDER]; 
	for(uint i=0; i<num_fields; i++)
	  sum[i]=0.0;
	
	/* Here I play with indeces: let alpha_idx and site_idx be the indices to cover the set of
	 * staggered fields. Then we have a superindex s = site_idx + EOPREC_SPINORFIELDSIZE_LOCAL * alpha_idx.
	 * It ranges between 0 and EOPREC_SPINORFIELDSIZE_LOCAL*num_fields. Then, we can use this superindex
	 * to parallelize the calculation on the GPU. To get the two indeces from the superindex, we use
	 *  - alpha_idx = s / EOPREC_SPINORFIELDSIZE_LOCAL
	 *  -  site_idx = s - (s / EOPREC_SPINORFIELDSIZE_LOCAL) * EOPREC_SPINORFIELDSIZE_LOCAL
	 * Here below, the superindex is called id_mem.
	 */
	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_LOCAL * num_fields; id_mem += global_size) {
		int alpha_idx = id_mem / EOPREC_SPINORFIELDSIZE_LOCAL;
		int  site_idx = id_mem - alpha_idx * EOPREC_SPINORFIELDSIZE_LOCAL;
		su3vec x_tmp = get_su3vec_from_field_eo(x, site_idx);
		x_tmp = su3vec_times_real(x_tmp, alpha[alpha_idx]);
		//sqnorms[alpha_idx] += su3vec_squarenorm(x_tmp);
		sum[alpha_idx] += su3vec_squarenorm(x_tmp);
	}
	
	/* At this point each thread has the sum vector with inside some numbers. I must now
	 * collect these numbers in the right way. The procedure is quite similar to the reduction
	 * used in the squarenorm kernel, but here result has num_fields components per each working
	 * group. It has then num_fields*num_groups components and we have to use again a superindex.
	 * It is s = alpha_idx + num_fields * group_idx and we have
	 *  - group_idx = s / num_fields
	 *  - alpha_idx = s - (s / num_fields) * num_fields
	 * 
	 * The same is valid for result_local that will be a vector with local_size*num_fields components.
	 * Its superindex will be s = alpha_idx + num_fields * local_idx and we have
	 *  - local_idx = s / num_fields
	 *  - alpha_idx = s - (s / num_fields) * num_fields
	 * 
	 * Here below, the group_idx is group_id and the local_idx is idx.
	 */
	if(local_size == 1) {
		for(uint i=0; i<num_fields; i++)
		  result[ i + num_fields * group_id ] = sum[i];
	} else {
		for(uint i=0; i<num_fields; i++)
		  result_local[i + num_fields * idx] = sum[i];
		// sync threads
		barrier(CLK_LOCAL_MEM_FENCE);

		//reduction until threads 0-7 hold all partial sums
		int cut1;
		int cut2 = local_size;
		for(cut1 = local_size / 2; cut1 > 4; cut1 /= 2) {
		  for(uint k=0; k<num_fields; k++){
		    for(int i = idx + cut1; i < cut2; i += cut1) {
		      result_local[k + num_fields * idx] +=  result_local[k + num_fields* i];
		    }
		  }
		  barrier(CLK_LOCAL_MEM_FENCE);
		  cut2 = cut1;
		}
		//thread 0 sums up the last 8 results and stores them in the global buffer
		if (idx == 0) {
		  for(uint i=0; i<num_fields; i++){
		    result[ i + num_fields * group_id ] =  result_local[i + num_fields * 0] + 
		                                           result_local[i + num_fields * 1] +
		                                           result_local[i + num_fields * 2] +
		                                           result_local[i + num_fields * 3] +
		                                           result_local[i + num_fields * 4] +
		                                           result_local[i + num_fields * 5] +
		                                           result_local[i + num_fields * 6] +
		                                           result_local[i + num_fields * 7];
		  }
		}
	}
}


// Description of variables of sax_vectorized_and_squarenorm_reduction kernel:
//  - out: The final result, namely a vector with num_fields components, each with a total squarenorm
//  - result_tmp: This is the vector result filled by the kernel sax_vectorized_and_squarenorm_eoprec
//  - elems: It is the number of components of result_tmp, namely the variable "num_groups" of the kernel
//           sax_vectorized_and_squarenorm_eoprec times num_fields
//  - num_fields: number of output fields

__kernel void sax_vectorized_and_squarenorm_reduction(__global hmc_float* out, __global hmc_float* result_tmp, const uint elems, const uint num_fields)
{
	const uint id = get_global_id(0);
	for(uint i=0; i<num_fields; i++)
		out[i]=0.0;
	  
	if(id == 0) {
		for (uint s = 0; s < elems; s++) {
		  /* We have s = field_idx + num_fields * group_idx and we have
		   *  - group_idx = s / num_fields
		   *  - field_idx = s - (s / num_fields) * num_fields
		   */
		  out[s - (s / num_fields) * num_fields] += result_tmp[s];
		}
	}
}







