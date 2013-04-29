/** @file
 * Plaquette calculation kernels
 * It is:
 *  Plaquette = \sum_{sites} \sum_{mu > nu} local_plaquette (i, mu, nu)
 * NOTE: The reduction used in this kernel is only safe with ls being a power of 2 and bigger than 8!
 */
__kernel void plaquette(__global Matrixsu3StorageType * field, __global hmc_float * plaq_out, __global hmc_float* tplaq_out, __global hmc_float* splaq_out, __local hmc_float * plaq_loc, __local hmc_float* tplaq_loc, __local hmc_float* splaq_loc)
{
	int id_local;
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

	for(id_local = id_tmp; id_local < VOL4D_LOCAL; id_local += global_size) {
		st_index pos = (id_local < VOL4D_LOCAL / 2) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local - (VOL4D_LOCAL / 2));

		for(int mu = 1; mu < NDIM; mu++) {
			for(int nu = 0; nu < mu; nu++) {
				prod = local_plaquette(field, pos.space, pos.time, mu, nu );

				tmpfloat = trace_matrixsu3(prod).re;
				tmpfloat /= NC;
				plaq += tmpfloat;

				if(nu == 0) {
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
	  	plaq_loc[idx] = plaq;
	  	tplaq_loc[idx] = tplaq;
	  	splaq_loc[idx] = splaq;
		// sync threads
		barrier(CLK_LOCAL_MEM_FENCE);

		//reduction until threads 0-7 hold all partial sums
		int cut1;
		int cut2 = local_size;
		for(cut1 = local_size / 2; cut1 > 4; cut1 /= 2) {
		  for(int i = idx + cut1; i < cut2; i += cut1) {
		     plaq_loc[idx] +=  plaq_loc[i];
		    tplaq_loc[idx] += tplaq_loc[i];
		    splaq_loc[idx] += splaq_loc[i];
		  }
		  barrier(CLK_LOCAL_MEM_FENCE);
		  cut2 = cut1;
		}
		//thread 0 sums up the last 8 results and stores them in the global buffer
		if (idx == 0) {
		  plaq_out[ group_id ] =  plaq_loc[0] + plaq_loc[1] + plaq_loc[2] + plaq_loc[3] +
		    plaq_loc[4] + plaq_loc[5] + plaq_loc[6] + plaq_loc[7];
		  tplaq_out[ group_id ] =  tplaq_loc[0] + tplaq_loc[1] + tplaq_loc[2] + tplaq_loc[3] +
		    tplaq_loc[4] + tplaq_loc[5] + tplaq_loc[6] + tplaq_loc[7];
		  splaq_out[ group_id ] =  splaq_loc[0] + splaq_loc[1] + splaq_loc[2] + splaq_loc[3] +
		    splaq_loc[4] + splaq_loc[5] + splaq_loc[6] + splaq_loc[7];
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
