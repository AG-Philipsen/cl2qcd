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

// This is the merge of the saxpy kernel and the squarenorm kernel
// saxpy:
// out = -alpha*x + y
//CP: defined with a minus!!!
// the squarenorm result |out|^2 is stored in "result"
__kernel void saxpy_AND_squarenorm_eo(__global const spinorStorageType * const restrict x, __global const spinorStorageType * const restrict y, __global const hmc_complex * const restrict alpha, __global spinorStorageType * const restrict out, __global hmc_float * const restrict result, __local hmc_float *\
 const restrict result_local )
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int group_id = get_group_id (0);
	int idx = get_local_id(0);

	hmc_float sum;
	sum = 0.;

	const hmc_complex alpha_tmp = complexLoadHack(alpha);
	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		spinor x_tmp = getSpinor_eo(x, id_mem);
		spinor y_tmp = getSpinor_eo(y, id_mem);
		x_tmp = spinor_times_complex(x_tmp, alpha_tmp);
		x_tmp = spinor_dim(y_tmp, x_tmp);
		//calc squarenorm of resulting spinor
//		hmc_float tmp = spinor_squarenorm(x_tmp);
//		sum += tmp;
sum +=		spinor_squarenorm(x_tmp);
		putSpinor_eo(out, id_mem, x_tmp);
	}

	//perform local reduction
	if(local_size == 1) {
		result[ group_id ] = sum;
	} else {
		// sync threads
		barrier(CLK_LOCAL_MEM_FENCE);
		//reduction
		(result_local[idx]) = sum;
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 64)
			result_local[idx % 64] += result_local[idx];
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 32)
			result_local[idx - 32] += result_local[idx];
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 16)
			result_local[idx - 16] += result_local[idx];
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 8)
			result_local[idx - 8] += result_local[idx];
		barrier(CLK_LOCAL_MEM_FENCE);
		//thread 0 sums up the result_local and stores it in array result
		if (idx == 0) {
			if(local_size >= 8) {
				result[ group_id ] =  result_local[0] + result_local[1] + result_local[2] + result_local[3] +
				                      result_local[4] + result_local[5] + result_local[6] + result_local[7];
			} else {
				for(int i = 0; i < local_size; i++)
					result[group_id] += result_local[i];
			}
		}
	}
	return;


}
