// alpha*x
__kernel void sax_eoprec(__global const spinorStorageType * const x, __global const hmc_complex * const alpha, __global spinorStorageType * const out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	hmc_complex alpha_tmp = complexLoadHack(alpha);
	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		spinor x_tmp = getSpinor_eo(x, id_mem);
		x_tmp = spinor_times_complex(x_tmp, alpha_tmp);
		putSpinor_eo(out, id_mem, x_tmp);
	}
}
