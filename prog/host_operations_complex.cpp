#include "host_operations_complex.h"

hmc_complex complexconj(hmc_complex *in){
  hmc_complex z = *in;
  z.im = -z.im;
  return z;
}

hmc_complex complexmult(hmc_complex *a, hmc_complex *b){
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

hmc_error complexaccumulate(hmc_complex *inout, hmc_complex *incr){
  (*inout).re += (*incr).re;
  (*inout).im += (*incr).im;
  return HMC_SUCCESS;
}

hmc_complex complexdivide(hmc_complex* numerator, hmc_complex* denominator){
  hmc_float norm = (*denominator).re*(*denominator).re + (*denominator).im*(*denominator).im;
  hmc_complex res;
  res.re = ((*numerator).re*(*denominator).re+(*numerator).im*(*denominator).im)/norm;
  res.im = ((*numerator).im*(*denominator).re-(*numerator).re*(*denominator).im)/norm;
  return res;
}

hmc_error complexcopy(hmc_complex* source, hmc_complex* dest, int length){
	// copies ``length'' complex numbers from source array to dest array, within cpu memory
	for(int i=0;i<length;i++){
		dest[i] = source[i];
	}
	return HMC_SUCCESS;  // SL: function not tested
}