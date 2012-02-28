//Here, two cases are possible:
//evenodd = ODD or EVEN
//	ODD corresponds to the D_oe case: Dslash acts on even indices (the "x+mu" in the formulae) and the function
//	saves the outcoming spinor with an odd index.
//	EVEN is then D_eo.
__kernel void dslash_eoprec(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, __global const hmc_complex * const restrict field, const int evenodd)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {
		st_idx pos = (evenodd == ODD) ? get_even_st_idx(id_tmp) : get_odd_st_idx(id_tmp);

		spinor out_tmp = set_spinor_zero();
		spinor out_tmp2;

		//calc dslash (this includes mutliplication with kappa)

		out_tmp2 = dslash_eoprec_unified_local(in, field, pos, TDIR);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_unified_local(in, field, pos, XDIR);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_unified_local(in, field, pos, YDIR);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_unified_local(in, field, pos, ZDIR);
		out_tmp = spinor_dim(out_tmp, out_tmp2);

		putSpinorSOA_eo(out, id_tmp, out_tmp);
	}
}

__kernel void convertGaugefieldToSOA(__global hmc_complex * const restrict out, __global const Matrixsu3 * const restrict in)
{
	// we need to take care of index converion. in the AOS storage the dimension is the continious index
	// in the soa storage we want the space indices to be continuous and have the dimension as outermost.
	for(uint d = 0; d < NDIM; ++d) {
		for(site_idx s = get_global_id(0); s < VOL4D; s += get_global_size(0)) {
			Matrixsu3 tmp = in[d + NDIM * s];

			const st_idx site = get_st_idx_from_site_idx(s);
			putSU3SOA(out, get_link_idx_SOA(d, site), tmp);
		}
	}
}
__kernel void convertGaugefieldFromSOA(__global Matrixsu3 * const restrict out, __global const hmc_complex * const restrict in)
{
	// we need to take care of index converion. in the AOS storage the dimension is the continious index
	// in the soa storage we want the space indices to be continuous and have the dimension as outermost.
	for(uint d = 0; d < NDIM; ++d) {
		for(site_idx s = get_global_id(0); s < VOL4D; s += get_global_size(0)) {
			const st_idx site = get_st_idx_from_site_idx(s);
			Matrixsu3 tmp = getSU3SOA(in, get_link_idx_SOA(d, site));

			out[d + NDIM * s] = tmp;
		}
	}
}
