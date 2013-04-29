/** @file
 * Rectangles calculation kernels.
 * This is very similar to the Plaquette kernel.
 * It is:
 *  Rectangles = \sum_{sites} \sum_{mu \neq nu} local_rectangles (i, mu, nu)
 * NOTE: The reduction used in this kernel is only safe with ls being a power of 2 and bigger than 8!
 */
__kernel void rectangles(__global Matrixsu3StorageType * field, __global hmc_float * rect_out, __local hmc_float * rect_loc)
{
	int id_local;
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id_tmp = get_global_id(0);
	int idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	if(idx == 0) {
		(rect_out)[group_id] = 0.0f;
	}

	hmc_float rect = 0;
	hmc_float tmpfloat = 0;

	Matrixsu3 prod;

	for(id_local = id_tmp; id_local < VOL4D_LOCAL; id_local += global_size) {
		st_index pos = (id_local < VOL4D_LOCAL / 2) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local - (VOL4D_LOCAL / 2));

		for(int mu = 0; mu < NDIM; mu++) {
			for(int nu = 0; nu < NDIM; nu++) {
				if(nu == mu) continue;
				prod = local_rectangles(field, pos.space, pos.time, mu, nu );
				tmpfloat = trace_matrixsu3(prod).re;
				rect += tmpfloat / NC;
			}
		}
	}

	if(local_size == 1) {
		rect_out[ group_id ]  += rect;
	} else {
	  	rect_loc[idx] = rect;
		// sync threads
		barrier(CLK_LOCAL_MEM_FENCE);

		//reduction until threads 0-7 hold all partial sums
		int cut1;
		int cut2 = local_size;
		for(cut1 = local_size / 2; cut1 > 4; cut1 /= 2) {
		  for(int i = idx + cut1; i < cut2; i += cut1) {
		    rect_loc[idx] +=  rect_loc[i];
		  }
		  barrier(CLK_LOCAL_MEM_FENCE);
		  cut2 = cut1;
		}
		//thread 0 sums up the last 8 results and stores them in the global buffer
		if (idx == 0) {
		  rect_out[ group_id ] =  rect_loc[0] + rect_loc[1] + rect_loc[2] + rect_loc[3] +
		    rect_loc[4] + rect_loc[5] + rect_loc[6] + rect_loc[7];
		}
	}
	return;
}

__kernel void rectangles_reduction(__global hmc_float* rect_buf, __global hmc_float* rect, const uint bufElems)
{
	int id = get_global_id(0);
	if(id == 0) {
		for (uint i = 1; i < bufElems; i++) {
			rect_buf[0]  += rect_buf[i];
		}
		(*rect)  = rect_buf[0];
	}

	return;
}
