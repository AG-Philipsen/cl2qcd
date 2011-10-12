__kernel void set_eoprec_spinorfield_cold(__global spinorfield* in)
{
	int id = get_global_id(0);
	if (id == 0) {
		//NOTE: the normalization is the same as in the above case, putting |phi|^2 of the WHOLE field to 1!!
		hmc_float norm = 1. / sqrt((hmc_float)(12. * VOLSPACE * NTIME));
		for(int n = 0; n < EOPREC_SPINORFIELDSIZE; n++) {
			in[n] = set_spinor_cold();
			in[n] = real_multiply_spinor(in[n], norm);
		}
	}
}
