/** @file
 * Device code implementing complex numbers
 */

//opencl_operations_complex.cl

hmc_complex complexconj(hmc_complex in)
{
	in.im = -(in.im);
	return in;
}

hmc_complex complexmult(const hmc_complex a, const hmc_complex b)
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


__kernel void ratio(__global hmc_complex * a, __global hmc_complex * b, __global hmc_complex * out)
{
    //!!CP: complexdivide cannot handle __global
    hmc_complex tmp1 = (*a);
    hmc_complex tmp2 = (*b);
    (*out) =  complexdivide(tmp1, tmp2);
    return;
}


__kernel void product(__global hmc_complex * a, __global hmc_complex * b, __global hmc_complex * out)
{
    //!!CP: complexdivide cannot handle __global
    hmc_complex tmp1 = (*a);
    hmc_complex tmp2 = (*b);
    (*out) =  complexmult(tmp1, tmp2);
    return;
}

//this is deprecated...
// void gaussianComplexVector(__global hmc_complex * vector, int length, hmc_float sigma, __global hmc_ocl_ran * rnd){
// 	// SL: this fills real and imaginary part of a vector of "length" complex numbers
// 	//     with components drawn with a Gaussian distribution and variance sigma
// 	for(int idx=0;idx<length;idx++){
// 		gaussianNormalPair(&vector[idx].re,&vector[idx].im);
// 		vector[idx].re*=sigma;
// 		vector[idx].im*=sigma;
// 	}
// 	return HMC_SUCCESS;
// 	// SL: not yet tested
// }

