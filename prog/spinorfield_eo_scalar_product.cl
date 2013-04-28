// complex (!!!) scalarproduct, return in result
// --> use 2 kernels: 1 for the summation in one block and 1 for summation over blockresults
/// NOTE: The reduction used in this kernel is only safe with ls being a power of 2 and bigger than 8!
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
	  	(result_local[idx]).re = sum.re;
		(result_local[idx]).im = sum.im;
		// sync threads
		barrier(CLK_LOCAL_MEM_FENCE);

		//reduction until threads 0-7 hold all partial sums
		int cut1;
		int cut2 = local_size;
		for(cut1 = local_size / 2; cut1 > 4; cut1 /= 2) {
		  for(int i = idx + cut1; i < cut2; i += cut1) {
		    result_local[idx].re += result_local[i].re;
		    result_local[idx].im += result_local[i].im;
		  }
		  barrier(CLK_LOCAL_MEM_FENCE);
		  cut2 = cut1;
		}
		//thread 0 sums up the last 8 results and stores them in the global buffer
		if (idx == 0) {
		  result[ group_id ].re =  result_local[0].re + result_local[1].re + result_local[2].re + result_local[3].re +
		    result_local[4].re + result_local[5].re + result_local[6].re + result_local[7].re;
		  result[ group_id ].im =  result_local[0].im + result_local[1].im + result_local[2].im + result_local[3].im +
		    result_local[4].im + result_local[5].im + result_local[6].im + result_local[7].im;
		}
	}
	return;
}
