__kernel void M_tm_sitediagonal_AND_gamma5_eo(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, hmc_float mubar_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	spinor out_tmp;
	spinor plus;

	//there is a sign change in the factors compared because of the gamma5
	hmc_complex twistfactor = {1., mubar_in};
	hmc_complex twistfactor_minus = {-1., 1.*mubar_in};

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = getSpinor_eo(in, id_mem);
		out_tmp = M_diag_tm_local(plus, twistfactor, twistfactor_minus);

		putSpinor_eo(out, id_mem, out_tmp);
	}
}

__kernel void M_tm_sitediagonal_minus_AND_gamma5_eo(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, hmc_float mubar_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	spinor out_tmp;
	spinor plus;

	//there is a sign change in the factors compared because of the gamma5
	hmc_complex twistfactor = {1., -mubar_in};
	hmc_complex twistfactor_minus = {-1., -1.*mubar_in};

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = getSpinor_eo(in, id_mem);
		out_tmp = M_diag_tm_local(plus, twistfactor, twistfactor_minus);

		putSpinor_eo(out, id_mem, out_tmp);
	}
}
