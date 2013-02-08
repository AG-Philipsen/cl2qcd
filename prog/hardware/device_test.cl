__kernel void foo(__global float * a, __global float * b, __global float * c)
{
	if(get_global_id(0) == 0) {
		*c = *a + *b;
	}
}
