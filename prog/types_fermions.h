/** @file
 * Common fermion types used by HMC, both on host and OpenCL device.
 */

#ifndef _TYPES_FERMIONSH_
#define _TYPES_FERMIONSH_

#include "types.h"

//CP: new defs for spinors
typedef struct {
	hmc_complex e0;
	hmc_complex e1;
	hmc_complex e2;
} su3vec
#ifdef _USEDOUBLEPREC_
__attribute__((aligned(16)));
#else
__attribute__((aligned(8)));
#endif

typedef struct {
	su3vec e0;
	su3vec e1;
	su3vec e2;
	su3vec e3;
} spinor
#ifdef _USEDOUBLEPREC_
__attribute__((aligned(32)));
#else
__attribute__((aligned(16)));
#endif

typedef struct {
	su3vec e0;
	su3vec e1;
} halfspinor;

/**
 * The type used for storing spinors on the device.
 */
#ifdef _USE_SOA_
typedef hmc_complex spinorStorageType;
typedef hmc_complex staggeredStorageType;
#else
typedef spinor spinorStorageType;
typedef su3vec staggeredStorageType;
#endif

#endif /* _TYPES_FERMIONSH_ */

