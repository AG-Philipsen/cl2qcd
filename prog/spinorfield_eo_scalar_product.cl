// complex (!!!) scalarproduct, return in result
// --> use 2 kernels: 1 for the summation in one block and 1 for summation over blockresults
__kernel void scalar_product_eoprec( __global const spinorStorageType  * const x, __global const spinorStorageType * const y, __global hmc_complex * const result, __local hmc_complex * const result_local )
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int group_id = get_group_id (0);
	int idx = get_local_id(0);

	hmc_complex sum;
	sum.re = 0.;
	sum.im = 0.;

	for(int id_local = id; id_local < EOPREC_SPINORFIELDSIZE_LOCAL; id_local += global_size) {
		site_idx id_mem = get_eo_site_idx_from_st_idx(get_even_st_idx_local(id_local));
		spinor x_tmp = getSpinor_eo(x, id_mem);
		spinor y_tmp = getSpinor_eo(y, id_mem);
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
			result_local[idx % 64].re += result_local[idx].re;
			result_local[idx % 64].im += result_local[idx].im;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 32) {
			result_local[idx - 32].re += result_local[idx].re;
			result_local[idx - 32].im += result_local[idx].im;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 16) {
			result_local[idx - 16].re += result_local[idx].re;
			result_local[idx - 16].im += result_local[idx].im;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx >= 8) {
			result_local[idx - 8].re += result_local[idx].re;
			result_local[idx - 8].im += result_local[idx].im;
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
