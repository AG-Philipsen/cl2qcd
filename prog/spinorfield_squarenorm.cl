// hmc_float squarenorm, return in result
// --> use 2 kernels: 1 for the summation in one block and 1 for summation over blockresults
__kernel void global_squarenorm(__global const spinor * const restrict x, __global hmc_float * const restrict result, __local hmc_float * const restrict result_local)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int idx = get_local_id(0);

	hmc_float sum;
	sum = 0.;

	for(int id_local = id; id_local < SPINORFIELDSIZE_LOCAL; id_local += global_size) {
		site_idx id_mem = get_site_idx((id_local % 2 == 0) ? get_even_st_idx_local(id_local /   2) : get_odd_st_idx_local(id_local / 2));
		spinor x_tmp = x[id_mem];
		hmc_float tmp = spinor_squarenorm(x_tmp);
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
