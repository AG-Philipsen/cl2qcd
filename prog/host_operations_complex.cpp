#include "host_operations_complex.h"
#include "host_random.h"

hmc_complex complexconj(hmc_complex *in){
  hmc_complex z = *in;
  z.im = -z.im;
  return z;
}

hmc_complex complexmult(const hmc_complex *a, const hmc_complex *b){
  hmc_complex res;
  res.re = (*a).re*(*b).re - (*a).im*(*b).im;
  res.im = (*a).im*(*b).re + (*a).re*(*b).im;
  return res;
}

hmc_complex complexadd(hmc_complex *a, hmc_complex *b){
  hmc_complex res;
  res.re = (*a).re + (*b).re;
  res.im = (*a).im + (*b).im;
  return res;
}

hmc_complex complexsubtract(hmc_complex *a, hmc_complex *b){
  hmc_complex res;
  res.re = (*a).re - (*b).re;
  res.im = (*a).im - (*b).im;
  return res;
}

void complexaccumulate(hmc_complex *inout, hmc_complex *incr){
  (*inout).re += (*incr).re;
  (*inout).im += (*incr).im;
  return;
}

hmc_complex complexdivide(hmc_complex* numerator, hmc_complex* denominator){
  hmc_float norm = (*denominator).re*(*denominator).re + (*denominator).im*(*denominator).im;
  hmc_complex res;
  res.re = ((*numerator).re*(*denominator).re+(*numerator).im*(*denominator).im)/norm;
  res.im = ((*numerator).im*(*denominator).re-(*numerator).re*(*denominator).im)/norm;
  return res;
}

/** @todo bitwise copy (memcpy) should work and be faster */
void complexcopy(hmc_complex* source, hmc_complex* dest, int length){
	// copies ``length'' complex numbers from source array to dest array, within cpu memory
	for(int i=0;i<length;i++){
		dest[i] = source[i];
	}
	return;  // SL: function not tested
}

/** @todo bitwise copy (memcpy) should work and be faster */
void hmc_floatcopy(hmc_float* source, hmc_float* dest, int length){
	// copies ``length'' complex numbers from source array to dest array, within cpu memory
	for(int i=0;i<length;i++){
		dest[i] = source[i];
	}
	return;  // SL: function not tested
}

//multiply complex number with real factor
void complexmult_real(hmc_complex *a, hmc_float *b){
  (*a).re *= (*b);
  (*a).im *= (*b);
  return;
}

void gaussianComplexVector(hmc_complex * vector, int length, hmc_float sigma){
	// SL: this fills real and imaginary part of a vector of "length" complex numbers
	//     with components drawn with a Gaussian distribution and variance sigma
	for(int idx=0;idx<length;idx++){
		gaussianNormalPair(&vector[idx].re,&vector[idx].im);
		vector[idx].re*=sigma;
		vector[idx].im*=sigma;
	}
	return;
	// SL: not yet tested
}

//CP: should this go into host_random.h ??
void gaussianNormalPair(hmc_float * z1, hmc_float * z2){
	// Box-Muller method, cartesian form, for extracting two independent normal standard real numbers
	hmc_float u1 = 1.0 - rnd.doub();
	hmc_float u2 = 1.0 - rnd.doub();
	hmc_float p  = sqrt(-2*log(u1));
	*z1 = p * cos(2*PI*u2);
	*z2 = p * sin(2*PI*u2);
	return;
	// SL: not yet tested
}
