/**
 @file types needed in the HMC-algorithm
*/

//CP: this can be done this way only if types.h has been cleaned up!

typedef struct {
  hmc_float e0;
 	hmc_float e1;
  hmc_float e2;
  hmc_float e3;
  hmc_float e4;
  hmc_float e5;
  hmc_float e6;
  hmc_float e7;
} ae;

// #define GAUGEMOMENTASIZE NDIM*VOL4D

// Definition of numeric constants for the symmetric structure constants d_ijk of su(3) suited for OpenCL
/** 1/2 */
#define F_1_2  0.5
/** 1/(2*sqrt(3)) */
#define F_1_2S3 (0.288675134594813)
/** 1/sqrt(3) */
#define F_1_S3  (0.577350269189626)