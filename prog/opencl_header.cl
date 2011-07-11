/** @file
 * Inclusion and definition of types and definitions required in the Device code.
 */

//opencl_header.cl

#pragma OPENCL EXTENSION cl_amd_printf : enable

#ifdef _USEDOUBLEPREC_
#ifdef _DEVICE_DOUBLE_EXTENSION_AMD_
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif
#ifdef _DEVICE_DOUBLE_EXTENSION_KHR_
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif
#endif


#include "globaldefs.h" //NDIM, NSPIN, NC
#include "types.h"

//!!CP: why is this here?
// #define VOLSPACE NSPACE*NSPACE*NSPACE
#define VOL4D VOLSPACE*NTIME

//for hmc_ocl_su3matrix
#ifdef _RECONSTRUCT_TWELVE_
#define SU3SIZE NC*(NC-1)
#define STAPLEMATRIXSIZE NC*NC
#else
#define SU3SIZE NC*NC
#define STAPLEMATRIXSIZE NC*NC
#endif
