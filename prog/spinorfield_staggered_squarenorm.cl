// hmc_float squarenorm, return in result
// --> use 2 kernels: 1 for the summation in one block and 1 for summation over blockresults

// Description of variables of global_squarenorm_staggered kernel:
//  - x: The staggered field (an su3vec per each site => vector of VOL4D components that are su3vec varibles)
//  - result: Vector of hmc_float that will contain the sums of the components of the result_local vectors.
//            In other words, each component of this vector will contain the sum of all squarenorms that have
//            been mapped to the threads within a group. Therefore result is a vector with num_groups components.
//  - result_local: Vector with local_size components. At the end of the local reduction, the sum of its first
//                  components will be the sum of all its components (and will be put in result[group_id]).
//                  Observe that some components of result_local can include the sum of squarenorms of
//                  several fields (if SPINORFIELDSIZE>global_size).

__kernel void global_squarenorm_staggered(__global const su3vec * const restrict x, __global hmc_float * const restrict result, __local hmc_float * const restrict result_local)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int idx = get_local_id(0);

	hmc_float sum;
	sum = 0.;

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {
		su3vec x_tmp = x[id_tmp];
		hmc_float tmp = su3vec_squarenorm(x_tmp);
		sum += tmp;
	}

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

// Description of variables of global_squarenorm_reduction kernel:
//  - dest: The final result, namely the sum of the squarenorms of all staggered fields.
//  - result_tmp: This is the vector result filled by the kernel global_squarenorm_staggered.
//  - elems: It is the number of components of the vector result_tmp, namely the variable "num_groups"
//            of the kernel global_squarenorm_staggered.

__kernel void global_squarenorm_reduction(__global hmc_float* dest, __global hmc_float* result_tmp, const uint elems)
{
	uint id = get_global_id(0);
	hmc_float tmp = 0;
	if(id == 0) {
		for (uint i = 0; i < elems; i++) {
			tmp += result_tmp[i];
		}
		*dest = tmp;
	}
}
