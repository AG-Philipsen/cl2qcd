/** @file
 * Some utility kernels for buffer managment
 */

__kernel void copy_16_bytes(__global float4 * const restrict to, __global const float4 * const restrict from)
{
	// one float4 is exactly 16 bytes, so there is only work for one thread.

	if(get_global_id(0) == 0) {
		to[0] = from[0];
	}
}
