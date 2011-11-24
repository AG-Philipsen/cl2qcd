__kernel void fillComplex(__global hmc_complex * const restrict out, const hmc_complex value)
{
	for(size_t i = get_global_id(0); i < VOL4D; i += get_global_size(0)) {
		out[i] = value;
	}
}

__kernel void readComplex(__global float2 * const restrict out, __global const hmc_complex * const restrict in)
{
	for(size_t i = get_global_id(0); i < VOL4D; i += get_global_size(0)) {
		hmc_complex tmp = in[i];
		out[i] = (float2) (tmp.re, tmp.im);
	}
}
