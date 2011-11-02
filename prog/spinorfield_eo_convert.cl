//eoprec operations
__kernel void convert_to_eoprec(__global spinorfield_eoprec* even, __global spinorfield_eoprec* odd, __global spinorfield* in)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	for(int n = id; n < EOPREC_SPINORFIELDSIZE; n+=global_size) {
		st_index pos = get_even_site(n);
		even[n] = in[get_global_pos(pos.space, pos.time)];
		pos = get_odd_site(n);
		odd[n] = in[get_global_pos(pos.space, pos.time)];
	}
	return;
}

__kernel void convert_from_eoprec(__global spinorfield_eoprec* even, __global spinorfield_eoprec* odd, __global spinorfield* out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	for(int n = id; n < EOPREC_SPINORFIELDSIZE; n+=global_size) {
			st_index pos = get_even_site(n);
			out[get_global_pos(pos.space, pos.time)] = even[n];
			pos = get_odd_site(n);
			out[get_global_pos(pos.space, pos.time)] = odd[n];
	}
	return;
}
