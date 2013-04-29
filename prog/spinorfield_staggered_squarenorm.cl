// hmc_float squarenorm, return in result
// --> use 2 kernels: 1 for the summation in one block and 1 for summation over blockresults

// Description of variables of global_squarenorm_staggered kernel:
//  - x: The staggered field (an su3vec per each site => vector of VOL4D components that are su3vec varibles)
//  - result: Vector of hmc_float that will contain the sums of the components of the result_local vectors.
//            In other words, each component of this vector will contain the sum of all squarenorms that have
//            been mapped to the threads within a group. Therefore result is a vector with num_groups components.
//  - result_local: Vector with local_size components. At the end of the local reduction, its first
//                  component will be the sum of all its components (and will be put in result[group_id]).
//                  Observe that some components of result_local can include the sum of squarenorms of
//                  several fields (if SPINORFIELDSIZE>global_size).

__kernel void global_squarenorm_staggered(__global const su3vec * const restrict x, __global hmc_float * const restrict result, __local hmc_float * const restrict result_local)
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
		su3vec x_tmp = x[id_mem];
		hmc_float tmp = su3vec_squarenorm(x_tmp);
		sum += tmp;
	}

	if(local_size == 1) {
	      result[ group_id ] = sum;
	} else {
	       result_local[idx] = sum;
	       // sync threads
	       barrier(CLK_LOCAL_MEM_FENCE);
	       //reduction until threads 0-7 hold all partial sums
	       int cut1;
	       int cut2 = local_size;
	       for(cut1 = local_size / 2; cut1 > 4; cut1 /= 2) {
	         for(int i = idx + cut1; i < cut2; i += cut1) {
		    result_local[idx] += result_local[i];
		 }
		 barrier(CLK_LOCAL_MEM_FENCE);
		 cut2 = cut1;
	       }
	       //thread 0 sums up the last 8 results and stores them in the global buffer
	       if (idx == 0) {
	         result[ group_id ] = result_local[0] + result_local[1] +
		 	 	      result_local[2] + result_local[3] +
				      result_local[4] + result_local[5] +
				      result_local[6] + result_local[7];
	       }
	}
	return;
}

// Description of variables of global_squarenorm_reduction kernel:
//  - dest: The final result, namely the sum of the squarenorms of all staggered fields.
//  - result_tmp: This is the vector result filled by the kernel global_squarenorm_staggered.
//  - elems: It is the number of components of the vector result_tmp, namely the variable "num_groups"
//            of the kernel global_squarenorm_staggered.

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
