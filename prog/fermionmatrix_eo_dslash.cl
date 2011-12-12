//Here, two cases are possible:
//evenodd = ODD or EVEN
//	ODD corresponds to the D_oe case: Dslash acts on even indices (the "x+mu" in the formulae) and the function
//	saves the outcoming spinor with an odd index.
//	EVEN is then D_eo.
__kernel void dslash_eoprec(__global const hmc_complex * const restrict in, __global hmc_complex * const restrict out, __global const hmc_complex * const restrict field, const int evenodd)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {
		st_index pos = (evenodd == ODD) ? get_even_site(id_tmp) : get_odd_site(id_tmp);

		spinor out_tmp = set_spinor_zero();
		spinor out_tmp2;

		//calc dslash (this includes mutliplication with kappa)

		out_tmp2 = dslash_eoprec_local_0(in, field, pos.space, pos.time);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_local_1(in, field, pos.space, pos.time);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_local_2(in, field, pos.space, pos.time);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_local_3(in, field, pos.space, pos.time);
		out_tmp = spinor_dim(out_tmp, out_tmp2);

		putSpinorSOA_eo(out, id_tmp, out_tmp);
	}
}

__kernel void convertSpinorfieldToSOA_eo(__global hmc_complex * const restrict out, __global const spinor * const restrict in)
{
	for(uint i = get_global_id(0); i < EOPREC_SPINORFIELDSIZE; i += get_global_size(0)) {
		putSpinorSOA_eo(out, i, in[i]);
	}
}
__kernel void convertSpinorfieldFromSOA_eo(__global spinor * const restrict out, __global const hmc_complex * const restrict in)
{
	for(uint i = get_global_id(0); i < EOPREC_SPINORFIELDSIZE; i += get_global_size(0)) {
		out[i] = getSpinorSOA_eo(in, i);
	}
}

__kernel void convertGaugefieldToSOA(__global hmc_complex * const restrict out, __global const Matrixsu3 * const restrict in)
{
	// we need to take care of index converion. in the AOS storage the dimension is the continious index
	// in the soa storage we want the space indices to be continuous and have the dimension as outermost.
	for(uint d = 0; d < NDIM; ++d) {
		for(uint s = get_global_id(0); s < VOL4D; s += get_global_size(0)) {
			Matrixsu3 tmp = in[d + NDIM * s];

			const st_index site = get_site(s);
			putSU3SOA(out, get_global_link_pos_SOA(d, site.space, site.time), tmp);
		}
	}
}
__kernel void convertGaugefieldFromSOA(__global Matrixsu3 * const restrict out, __global const hmc_complex * const restrict in)
{
	// we need to take care of index converion. in the AOS storage the dimension is the continious index
	// in the soa storage we want the space indices to be continuous and have the dimension as outermost.
	for(uint d = 0; d < NDIM; ++d) {
		for(uint s = get_global_id(0); s < VOL4D; s += get_global_size(0)) {
			const st_index site = get_site(s);
			Matrixsu3 tmp = getSU3SOA(in, get_global_link_pos_SOA(d, site.space, site.time));

			out[d + NDIM * s] = tmp;
		}
	}
}
