/** @file
 * Kernel for bandwidth microbenchmark
 */

__kernel void copyFloat( __global hmc_float * const restrict out, __global const hmc_float * const restrict in, const ulong elems, const ulong threads_per_group )
{
	size_t local_id = get_local_id(0);

	if( local_id >= threads_per_group )
		return;

	size_t id = get_group_id(0) * threads_per_group + local_id;

	for( size_t i = id; i < elems;  i += threads_per_group * get_num_groups(0) ) {
		hmc_float tmp = in[i];
		out[i] = tmp;
	}
}

__kernel void copySU3( __global Matrixsu3 * const restrict out, __global const Matrixsu3 * const restrict in, const ulong elems, const ulong threads_per_group )
{
	size_t local_id = get_local_id(0);

	if( local_id >= threads_per_group )
		return;

	size_t id = get_group_id(0) * threads_per_group + local_id;

	for( size_t i = id; i < elems;  i += threads_per_group * get_num_groups(0) ) {
		Matrixsu3 tmp = in[i];
		out[i] = tmp;
	}
}

