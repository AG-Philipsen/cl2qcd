
//Here, two cases are possible:
//evenodd = ODD or EVEN
//	ODD corresponds to the D_oe case: Dslash acts on even indices (the "x+mu" in the formulae) and the function
//	saves the outcoming spinor with an odd index.
//	EVEN is then D_eo.
__kernel void dslash_eoprec(__global spinorfield_eoprec* in, __global spinorfield_eoprec* out, __global ocl_s_gaugefield* field, int evenodd)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int n, t;
#ifndef _USEGPU_
	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp +=	global_size) {
#else
	int id_tmp = id;
	if (id >= EOPREC_SPINORFIELDSIZE) return;
#endif
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
		put_spinor_to_eoprec_field(out_tmp, out, id_tmp);
#ifndef _USEGPU_
	}
#endif
}
