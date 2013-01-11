inline hmc_complex complexconj(hmc_complex in)
{
	in.im = -(in.im);
	return in;
}

inline hmc_complex complexmult(const hmc_complex a, const hmc_complex b)
{
	hmc_complex res;
	res.re = a.re * b.re - a.im * b.im;
	res.im = a.im * b.re + a.re * b.im;
	return res;
}

inline hmc_complex complexadd(const hmc_complex a, const hmc_complex b)
{
	hmc_complex res;
	res.re = a.re + b.re;
	res.im = a.im + b.im;
	return res;
}

inline hmc_complex complexsubtract(const hmc_complex a, const hmc_complex b)
{
	hmc_complex res;
	res.re = a.re - b.re;
	res.im = a.im - b.im;
	return res;
}

inline hmc_complex complexdivide(const hmc_complex numerator, const hmc_complex denominator)
{
	hmc_float norm = denominator.re * denominator.re + denominator.im * denominator.im;
	hmc_complex res;
	res.re = (numerator.re * denominator.re + numerator.im * denominator.im ) / norm;
	res.im = (numerator.im * denominator.re - numerator.re * denominator.im ) / norm;
	return res;
}
