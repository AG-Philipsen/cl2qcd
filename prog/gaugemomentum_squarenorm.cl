/** @todo add args for reduction... */
__kernel void gaugemomentum_squarenorm(__global const aeStorageType * const restrict in, __global hmc_float * const restrict out, __local hmc_float * const restrict result_local)
{
	hmc_float sum = 0.f;
	PARALLEL_FOR(i, GAUGEMOMENTASIZE) {
		sum += ae_squarenorm(getAe(in, i));
	}

	const size_t group_id = get_group_id(0);
	const size_t idx = get_local_id(0);
	const size_t local_size = get_local_size(0);

	if(local_size == 1) {
		out[ group_id ] = sum;
	} else {
		// FIXME ensure local size is 128!!!
		// sync threads
		barrier(CLK_LOCAL_MEM_FENCE);
		//reduction
		result_local[idx] = sum;
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx < 64)
			result_local[idx] += result_local[idx + 64];
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx < 32)
			result_local[idx] += result_local[idx + 32];
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx < 16)
			result_local[idx] += result_local[idx + 16];
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx < 8)
			result_local[idx] += result_local[idx + 8];
		barrier(CLK_LOCAL_MEM_FENCE);
		//thread 0 sums up the result_local and stores it in array result
		if (idx == 0) {
			if(local_size >= 8) {
				out[ group_id ] =  result_local[0] + result_local[1] + result_local[2] + result_local[3] +
				                   result_local[4] + result_local[5] + result_local[6] + result_local[7];
			} else {
				for(size_t i = 0; i < local_size; i++)
					out[group_id] += result_local[i];
			}
		}
	}
}
