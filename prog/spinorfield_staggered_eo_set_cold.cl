__kernel void set_cold_spinorfield_stagg_eoprec(__global staggeredStorageType * const restrict out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	//NOTE: the normalization is the same as in the above case, putting |phi|^2 of the WHOLE field to 1!!
	hmc_float norm = 1. / sqrt((hmc_float)(3. * VOL4D_GLOBAL));
	for(int n = id; n < EOPREC_SPINORFIELDSIZE_MEM; n += global_size) {
		put_su3vec_to_field_eo(out, n, su3vec_times_real(set_su3vec_cold(), norm));
	}
}
