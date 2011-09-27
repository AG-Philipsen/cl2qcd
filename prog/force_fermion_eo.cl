
__kernel void fermion_force_eoprec(__global ocl_s_gaugefield * field, __global  spinorfield * Y, __global  spinorfield * X, __global  ae * out, int evenodd)
{
	int id = get_global_id(0);
//CP: this is just a dummy kernel at the moment
//	basically it needs the same code as the non-eoprec version,
//	but modified with respect to the even-odd distinction (like in the dslash-case)
	return;
}
