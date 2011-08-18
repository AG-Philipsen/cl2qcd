/**
 * @file stout-smearing of the gaugefield
 */

__kernel void stout_smear(int iterations, hmc_float rho, __global ocl_s_gaugefield * u_inout)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	
	//DUMMY!!!!
	int n, t, index;
	Matrixsu3 tmp;
	Matrixsu3 tmp2;
	
	for(int id_tmp = id; id_tmp < VOL4D; id_tmp += global_size) {
		/** @todo this must be done more efficient */
		if(id_tmp < VOLSPACE * NTIME / 2)
			get_even_site(id_tmp, &n, &t);
		else
			get_odd_site(id_tmp-(VOLSPACE*NTIME/2), &n, &t);
		
		for(int mu = 0; mu < NDIM; mu++) {
			index = get_global_link_pos(mu, n, t);
			// an su3 algebra element has NC*NC-1 = 8 hmc_float entries
			// &(p_in[index*8]) should point to the right position for the pos-th element of the long gaugemomentum vector p_in
// 			tmp2 = build_su3matrix_by_exponentiation((p_in[index]), eps);
// 
// 			tmp = get_matrixsu3(u_inout, n, t, mu);
// 			tmp2 = multiply_matrixsu3( tmp2, tmp);
			put_matrixsu3(u_inout, tmp2, n, t, mu);
		}
	}
}
