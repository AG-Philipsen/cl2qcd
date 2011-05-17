/** @file
 * Device code implementing complex numbers
 */

//opencl_operations_complex.cl

hmc_complex complexconj(hmc_complex in)
{
	in.im = -(in.im);
	return in;
}

hmc_complex complexmult(hmc_complex a, hmc_complex b)
{
	hmc_complex res;
	res.re = a.re * b.re - a.im * b.im;
	res.im = a.im * b.re + a.re * b.im;
	return res;
}

hmc_complex complexadd(hmc_complex a, hmc_complex b)
{
	hmc_complex res;
	res.re = a.re + b.re;
	res.im = a.im + b.im;
	return res;
}

hmc_complex complexsubtract(hmc_complex a, hmc_complex b)
{
	hmc_complex res;
	res.re = a.re - b.re;
	res.im = a.im - b.im;
	return res;
}

hmc_complex complexdivide(hmc_complex numerator, hmc_complex denominator)
{
	hmc_float norm = denominator.re * denominator.re + denominator.im * denominator.im;
	hmc_complex res;
	res.re = (numerator.re * denominator.re + numerator.im * denominator.im ) / norm;
	res.im = (numerator.im * denominator.re - numerator.re * denominator.im ) / norm;
	return res;
}

