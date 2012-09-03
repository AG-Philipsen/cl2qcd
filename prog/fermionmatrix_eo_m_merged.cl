__kernel void M_tm_sitediagonal_AND_gamma5(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, hmc_float mubar_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	spinor out_tmp;
	spinor plus;

	hmc_complex twistfactor = {1., mubar_in};
	hmc_complex twistfactor_minus = {1., -1.*mubar_in};

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = getSpinor_eo(in, id_tmp);
		//Diagonalpart:
		out_tmp = M_diag_tm_local(plus, twistfactor, twistfactor_minus);
		
		//gamma5 part
		out_tmp = gamma5_local(out_tmp);

		putSpinor_eo(out, id_tmp, out_tmp);
	}
}

__kernel void M_tm_sitediagonal_minus_AND_gamma5(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, hmc_float mubar_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	spinor out_tmp;
	spinor plus;

	hmc_complex twistfactor = {1., mubar_in};
	hmc_complex twistfactor_minus = {1., -1.*mubar_in};

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = getSpinor_eo(in, id_tmp);
		//Diagonalpart:
		out_tmp = M_diag_tm_local(plus, twistfactor_minus, twistfactor);

		//gamma5 part
		out_tmp = gamma5_local(out_tmp);

		putSpinor_eo(out, id_tmp, out_tmp);
	}
}
