/** @file
 * Plaquette calculation kernels
 */

__kernel void plaquette(__global ocl_s_gaugefield * field, __global hmc_float * plaq_out, __global hmc_float* tplaq_out, __global hmc_float* splaq_out, __local hmc_float * plaq_loc, __local hmc_float* tplaq_loc, __local hmc_float* splaq_loc)
{

	int t, pos, id;
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id_tmp = get_global_id(0);
	int idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	if(idx == 0) {
		(plaq_out)[group_id] = 0.0f;
		(splaq_out)[group_id] = 0.0f;
		(tplaq_out)[group_id] = 0.0f;
	}

	hmc_float plaq = 0;
	hmc_float splaq = 0;
	hmc_float tplaq = 0;
	hmc_float tmpfloat = 0;

	Matrixsu3 prod;

	for(id = id_tmp; id < VOLSPACE * NTIME; id += global_size) {
		if(id < VOLSPACE * NTIME / 2)
			get_even_site(id, &pos, &t);
		else
			get_odd_site(id - (VOLSPACE * NTIME / 2), &pos, &t);

		for(int mu = 0; mu < NDIM; mu++) {
			for(int nu = 0; nu < mu; nu++) {
				prod = local_plaquette(field, pos, t, mu, nu );
				tmpfloat = trace_matrixsu3(prod).re;
				plaq += tmpfloat;
				if(mu == 0 || nu == 0) {
					tplaq += tmpfloat;
				} else {
					splaq += tmpfloat;
				}
			}
		}
	}

	if(local_size == 1) {
		plaq_out[ group_id ]  += plaq;
		tplaq_out[ group_id ] += tplaq;
		splaq_out[ group_id ] += splaq;
	} else {
		barrier(CLK_LOCAL_MEM_FENCE);
		//reduction
		plaq_loc[ idx ]  = plaq;
		tplaq_loc[ idx ] = tplaq;
		splaq_loc[ idx ] = splaq;

		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 64) {
			plaq_loc[ idx%64 ]  += plaq_loc[ idx ];
			tplaq_loc[ idx%64 ] += tplaq_loc[ idx ];
			splaq_loc[ idx%64 ] += splaq_loc[ idx ];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 32) {
			plaq_loc[ idx-32 ]  += plaq_loc[ idx ];
			tplaq_loc[ idx-32 ] += tplaq_loc[ idx ];
			splaq_loc[ idx-32 ] += splaq_loc[ idx ];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 16) {
			plaq_loc[ idx-16 ]  += plaq_loc[ idx ];
			tplaq_loc[ idx-16 ] += tplaq_loc[ idx ];
			splaq_loc[ idx-16 ] += splaq_loc[ idx ];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 8) {
			plaq_loc[ idx-8 ]  += plaq_loc[ idx ];
			tplaq_loc[ idx-8 ] += tplaq_loc[ idx ];
			splaq_loc[ idx-8 ] += splaq_loc[ idx ];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//thread 0 sums up the result_local and stores it in array result
		if (idx == 0) {
			if(local_size >= 8) {
				for (int i = 0; i < 8; i++) {
					plaq_out[ group_id ]  += plaq_loc[ i ];
					tplaq_out[ group_id ] += tplaq_loc[ i ];
					splaq_out[ group_id ] += splaq_loc[ i ];
				}
			} else {
				for (int i = 0; i < local_size; i++) {
					plaq_out[ group_id ]  += plaq_loc[ i ];
					tplaq_out[ group_id ] += tplaq_loc[ i ];
					splaq_out[ group_id ] += splaq_loc[ i ];
				}
			}
		}
	}
	return;
}

__kernel void plaquette_reduction(__global hmc_float* plaq_buf, __global hmc_float* tplaq_buf, __global hmc_float* splaq_buf, __global hmc_float* plaq, __global hmc_float* tplaq, __global hmc_float* splaq, const uint bufElems)
{
	int id = get_global_id(0);
	if(id == 0) {
		for (uint i = 1; i < bufElems; i++) {
			plaq_buf[0]  += plaq_buf[i];
			tplaq_buf[0] += tplaq_buf[i];
			splaq_buf[0] += splaq_buf[i];
		}
		(*plaq)  = plaq_buf[0];
		(*tplaq) = tplaq_buf[0];
		(*splaq) = splaq_buf[0];
	}

	return;
}
