// The variable num_values could not be passed to this function since it is
// equal to get_num_groups(0). Nevertheless, since in get_work_size ls, gs
// and num_groups are checked to be corrected and since it is useful
// to have num_groups variable to allocate buffers, we decided to use here
// the num_values variable.
__kernel void scalar_product_reduction(__global hmc_complex* result_tmp, __global hmc_complex* result, const uint num_values)
{
	//!!CP: complex_acc cannot handle __global
	int id = get_global_id(0);
	if(id == 0) {
		hmc_complex tmp1;
		hmc_complex tmp2;
		tmp2 = complexLoadHack(&result_tmp[0]);
		for (int i = 1; i < num_values; i++) {
			tmp1 = complexLoadHack(&result_tmp[i]);
			tmp2 = complexadd(tmp2, tmp1);
		}
		(*result) = tmp2;
	}
	return;
}