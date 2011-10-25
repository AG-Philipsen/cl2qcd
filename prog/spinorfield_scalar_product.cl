// complex (!!!) scalarproduct, return in result
// --> use 2 kernels: 1 for the summation in one block and 1 for summation over blockresults
__kernel void scalar_product( __global spinorfield *x, __global spinorfield *y, __global hmc_complex* result, __local hmc_complex* result_local )
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int idx = get_local_id(0);

	hmc_complex sum;
	sum.re = 0.;
	sum.im = 0.;

	//!! CP: perhaps here one can first copy a whole spinor and then do the dot-prod
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {
		spinor x_tmp = x[id_tmp];
		spinor y_tmp = y[id_tmp];
		hmc_complex tmp = spinor_scalarproduct(x_tmp, y_tmp);
		sum.re += tmp.re;
		sum.im += tmp.im;
	}


	if(local_size == 1) {
		result[ group_id ].re = sum.re;
		result[ group_id ].im = sum.im;
	} else {
		// sync threads
		barrier(CLK_LOCAL_MEM_FENCE);
		//reduction
		(result_local[idx]).re = sum.re;
		(result_local[idx]).im = sum.im;

		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 64) {
			result_local[idx%64].re += result_local[idx].re;
			result_local[idx%64].im += result_local[idx].im;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 32) {
			result_local[idx-32].re += result_local[idx].re;
			result_local[idx-32].im += result_local[idx].im;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 16) {
			result_local[idx-16].re += result_local[idx].re;
			result_local[idx-16].im += result_local[idx].im;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 8) {
			result_local[idx-8].re += result_local[idx].re;
			result_local[idx-8].im += result_local[idx].im;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//thread 0 sums up the result_local and stores it in array result
		if (idx == 0) {
			if(local_size >= 8) {
				for (int i = 1; i < 8; i++) {
					result_local[idx].re += result_local[i].re;
					result_local[idx].im += result_local[i].im;
				}
				result[ group_id ].re = result_local[idx].re;
				result[ group_id ].im = result_local[idx].im;
			} else {
				for (int i = 1; i < local_size; i++) {
					result_local[idx].re += result_local[i].re;
					result_local[idx].im += result_local[i].im;
				}
				result[ group_id ].re = result_local[idx].re;
				result[ group_id ].im = result_local[idx].im;
			}
		}
	}
	return;
}


__kernel void scalar_product_reduction(__global hmc_complex* result_tmp, __global hmc_complex* result)
{
	//!!CP: complex_acc cannot handle __global
	int id = get_global_id(0);
	if(id == 0) {
		hmc_complex tmp1;
		hmc_complex tmp2;
		tmp2 = complexLoadHack(&result_tmp[0]);
		for (int i = 1; i < get_num_groups(0); i++) {
			tmp1 = complexLoadHack(&result_tmp[i]);
			tmp2 = complexadd(tmp2, tmp1);
		}
		(*result) = tmp2;
	}
	return;
}
