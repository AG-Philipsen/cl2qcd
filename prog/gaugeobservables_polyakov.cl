/** @file
 * Polyakov calculation kernels
 */

__kernel void polyakov_reduction(__global hmc_complex* poly_buf,  __global hmc_complex* poly, const uint bufElems)
{
	int id = get_global_id(0);
	if(id == 0) {
		hmc_complex sum = hmc_complex_zero;
		for (int i = 0; i < bufElems; i++) {
			sum.re += poly_buf[i].re;
			sum.im += poly_buf[i].im;
		}
		*poly = sum;
	}

	return;
}

__kernel void polyakov(__global Matrixsu3StorageType * field, __global hmc_complex * out, __local hmc_complex * out_loc)
{

	int id;
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id_tmp = get_global_id(0);
	int idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	hmc_complex tmp_pol;
	hmc_complex tmpcomplex;
	tmp_pol.re = 0.;
	tmp_pol.im = 0.;

	if(idx == 0) {
		out[group_id].re = 0.0f;
		out[group_id].im = 0.0f;
	}

	for(id = id_tmp; id < VOLSPACE; id += global_size) {
		Matrixsu3 prod;
		prod = local_polyakov(field, id);
		tmpcomplex = trace_matrixsu3(prod);
		(tmp_pol).re += tmpcomplex.re / NC;
		(tmp_pol).im += tmpcomplex.im / NC;
	}

	//reduction
	if(local_size == 1) {
		((out))[group_id].re += tmp_pol.re;
		((out))[group_id].im += tmp_pol.im;
	} else {
		//wait for all threads to end calculations
		barrier(CLK_LOCAL_MEM_FENCE);

		//!!CP: this should be checked by someone else than me
		//perform reduction
		out_loc[idx].re = tmp_pol.re;
		out_loc[idx].im = tmp_pol.im;
		int cut1;
		int cut2 = local_size;
		for(cut1 = local_size / 2; cut1 > 0; cut1 /= 2) {
			for(int i = idx + cut1; i < cut2; i += cut1) {
				((out_loc)[idx]).re +=  ((out_loc)[i]).re;
				((out_loc)[idx]).im +=  ((out_loc)[i]).im;
			}
			//!!CP: is this dangerous inside a for-loop?
			barrier(CLK_LOCAL_MEM_FENCE);
			cut2 = cut1;
		}
		if(idx == 0) {
			out[group_id].re = out_loc[0].re;
			out[group_id].im = out_loc[0].im;
		}
	}
	return;
}
