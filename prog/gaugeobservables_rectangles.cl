/** @file
 * Plaquette calculation kernels
 */

__kernel void rectangles(__global ocl_s_gaugefield * field, __global hmc_float * rect_out, __local hmc_float * rect_loc)
{

	int id;
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

	for(id = id_tmp; id < VOLSPACE * NTIME; id += global_size) {
		st_index pos = (id < VOLSPACE * NTIME / 2) ? get_even_site(id) : get_odd_site(id - (VOLSPACE * NTIME / 2));

		for(int mu = 1; mu < NDIM; mu++) {
			for(int nu = 0; nu < mu; nu++) {
				prod = local_rectangles(field, pos.space, pos.time, mu, nu );
				tmpfloat = trace_matrixsu3(prod).re;
				rect += tmpfloat;
			}
		}
	}

	if(local_size == 1) {
		rect_out[ group_id ]  += rect;
	} else {
		barrier(CLK_LOCAL_MEM_FENCE);
		//reduction
		rect_loc[ idx ]  = rect;

		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 64) {
			rect_loc[ idx%64 ]  += rect_loc[ idx ];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 32) {
			rect_loc[ idx-32 ]  += rect_loc[ idx ];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 16) {
			rect_loc[ idx-16 ]  += rect_loc[ idx ];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 8) {
			rect_loc[ idx-8 ]  += rect_loc[ idx ];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//thread 0 sums up the result_local and stores it in array result
		if (idx == 0) {
			if(local_size >= 8) {
				for (int i = 0; i < 8; i++) {
					rect_out[ group_id ]  += rect_loc[ i ];
				}
			} else {
				for (int i = 0; i < local_size; i++) {
					rect_out[ group_id ]  += rect_loc[ i ];
				}
			}
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
