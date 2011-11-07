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

__kernel void copySpinor( __global spinor * const restrict out, __global const spinor * const restrict in, const ulong elems, const ulong threads_per_group )
{
	size_t local_id = get_local_id(0);

	if( local_id >= threads_per_group )
		return;

	size_t id = get_group_id(0) * threads_per_group + local_id;

	for( size_t i = id; i < elems;  i += threads_per_group * get_num_groups(0) ) {
		spinor tmp = in[i];
		out[i] = tmp;
	}
}

// FIXME for now simply assume we have 128 threads
__attribute__((reqd_work_group_size(128, 1, 1)))
__kernel void copySpinorLocal( __global spinor * const restrict out, __global const spinor * const restrict in, const ulong elems, const ulong threads_per_group )
{
	__local spinor scratch[128];
	__local double2 * d_scratch = (__local double2 *) scratch;

	__global const double2 * d_in = (__global double2 *) in;
	__global const double2 * d_out = (__global double2 *) out;

	size_t local_id = get_local_id(0);
	size_t id = get_group_id(0) * threads_per_group + local_id;

	for(size_t group_offset = get_group_id(0) * 128; group_offset < elems; group_offset += get_global_size(0)) {
		// make sure you don't run off the end
		size_t elem_count = 128 * sizeof(spinor) / sizeof(double2);
		event_t copy_in_event = async_work_group_copy(d_scratch, &d_in[group_offset], elem_count, 0);
		event_t copy_out_event = async_work_group_copy(&d_out[group_offset], d_scratch, elem_count, copy_in_event);
		wait_group_events(1, copy_out_event);
	}
}
