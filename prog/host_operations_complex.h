#ifndef _OPERATIONS_COMPLEXH_
#define _OPERATIONS_COMPLEXH_
#include <iostream>
#include <math.h>
#include "host_random.h"
#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"

hmc_complex complexconj(hmc_complex *in); 
hmc_complex complexmult(hmc_complex *a, hmc_complex *b); 
hmc_complex complexadd(hmc_complex *a, hmc_complex *b);
hmc_complex complexsubtract(hmc_complex *a, hmc_complex *b); 
hmc_error complexaccumulate(hmc_complex *inout, hmc_complex *incr); 
hmc_complex complexdivide(hmc_complex* numerator, hmc_complex* denominator);
hmc_error complexcopy(hmc_complex* source, hmc_complex* dest, int length);
hmc_error complexmult_real(hmc_complex *a, hmc_float *b); 
hmc_error gaussianComplexVector(hmc_complex * vector, int length, hmc_float sigma);
void gaussianNormalPair(hmc_float * z1, hmc_float * z2);
#endif