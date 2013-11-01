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