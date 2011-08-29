__kernel void extendKernel(__global Matrixsu3 * out, __global Matrixsu2 * in, __global int * rand, const ulong num_elems) {
	const size_t id = get_global_id(0);
	if( id >= num_elems )
		return;

	out[id] = extend(rand[id], in[id]);
}
