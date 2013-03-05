//eoprec operations
__kernel void convert_to_eoprec(__global spinor * const restrict even, __global spinor * const restrict odd, __global spinor const * const restrict in)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	for(int n = id; n < EOPREC_SPINORFIELDSIZE; n += global_size) {
		st_index pos = get_even_site(n);
		even[n] = in[get_site_idx(pos)];
		pos = get_odd_site(n);
		odd[n] = in[get_site_idx(pos)];
	}
	return;
}

__kernel void convert_from_eoprec(__global spinor const * const restrict even, __global spinor const * const restrict odd, __global spinor * const restrict out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	for(int n = id; n < EOPREC_SPINORFIELDSIZE; n += global_size) {
		st_index pos = get_even_site(n);
		out[get_site_idx(pos)] = even[n];
		pos = get_odd_site(n);
		out[get_site_idx(pos)] = odd[n];
	}
	return;
}

__kernel void convertSpinorfieldToSOA_eo(__global spinorStorageType * const restrict out, __global const spinor * const restrict in)
{
	for(uint i = get_global_id(0); i < EOPREC_SPINORFIELDSIZE; i += get_global_size(0)) {
		putSpinor_eo(out, i, in[i]);
	}
}

__kernel void convertSpinorfieldFromSOA_eo(__global spinor * const restrict out, __global const spinorStorageType * const restrict in)
{
	for(uint i = get_global_id(0); i < EOPREC_SPINORFIELDSIZE; i += get_global_size(0)) {
		out[i] = getSpinor_eo(in, i);
	}
}
