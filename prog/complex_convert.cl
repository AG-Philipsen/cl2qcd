/** @file
 * Conversion of complex numbers
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

/**
 * Converts a float to a complex number
 */
__kernel void convert_float_to_complex(__global hmc_float const * const restrict in, __global hmc_complex * const restrict out)
{
	if(get_global_id(0) == 0) {
		hmc_complex tmp;
		tmp.re = *in;
		tmp.im = 0.f;
		*out = tmp;
	}
}

