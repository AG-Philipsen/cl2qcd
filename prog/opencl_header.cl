#ifdef _USEDOUBLEPREC_
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
//#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

#include "globaldefs.h" //NDIM, NSPIN, NC
#include "types.h"

//for hmc_ocl_su3matrix
#ifdef _RECONSTRUCT_TWELVE_
#define SU3SIZE NC*(NC-1)
#define STAPLEMATRIXSIZE NC*NC
#else
#define SU3SIZE NC*NC
#define STAPLEMATRIXSIZE NC*NC
#endif


#ifdef _USEDOUBLEPREC_
__constant hmc_float const projectioneps = 10.e-12;
#else
__constant hmc_float const projectioneps = 10.e-6;
#endif


