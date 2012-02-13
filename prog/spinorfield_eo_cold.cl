__kernel void set_eoprec_spinorfield_cold(__global hmc_complex * const restrict out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	//NOTE: the normalization is the same as in the above case, putting |phi|^2 of the WHOLE field to 1!!
	hmc_float norm = 1. / sqrt((hmc_float)(12. * VOLSPACE * NTIME));
	for(int n = id; n < EOPREC_SPINORFIELDSIZE; n += global_size) {
		putSpinorSOA_eo(out, n, real_multiply_spinor(set_spinor_cold(), norm));
	}
}
