__kernel void dummyKernel(__global float * const restrict out, __global const float * const restrict in)
{
	size_t id = get_global_id(0);
	out[id] = dummyFunction(in[id]);
}
