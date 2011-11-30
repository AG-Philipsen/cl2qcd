//Here, two cases are possible:
//evenodd = ODD or EVEN
//	ODD corresponds to the D_oe case: Dslash acts on even indices (the "x+mu" in the formulae) and the function
//	saves the outcoming spinor with an odd index.
//	EVEN is then D_eo.
#ifdef _USEGPU_
__attribute__((reqd_work_group_size(128, 1, 1)))
#endif
__kernel void dslash_eoprec(__global const spinorfield_eoprec * const restrict in, __global const spinorfield_eoprec * const restrict out, __global const ocl_s_gaugefield * const restrict field, const int evenodd)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {
		st_idx pos = (evenodd == ODD) ? get_even_st_idx(id_tmp) : get_odd_st_idx(id_tmp);

		spinor out_tmp = set_spinor_zero();
		spinor out_tmp2;

		//calc dslash (this includes mutliplication with kappa)

		out_tmp2 = dslash_eoprec_local_0(in, field, pos);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_local_1(in, field, pos);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_local_2(in, field, pos);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_local_3(in, field, pos);
		out_tmp = spinor_dim(out_tmp, out_tmp2);

		put_spinor_to_eoprec_field(out_tmp, out, id_tmp);
	}
}
