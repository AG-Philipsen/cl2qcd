/** @file
 * Common fermion types used by HMC, both on host and OpenCL device.
 */

#ifndef _TYPES_FERMIONSH_
#define _TYPES_FERMIONSH_

#include "types.h"

#ifdef _INKERNEL_
//CP: new defs for spinors
typedef struct {
	hmc_complex e0;
	hmc_complex e1;
	hmc_complex e2;
	
} su3vec;

typedef struct {
	su3vec e0;
	su3vec e1;
	su3vec e2;
	su3vec e3;
} spinor;

typedef struct {
	su3vec e0;
	su3vec e1;
} halfspinor;

typedef spinor spinorfield;
typedef spinor spinorfield_eoprec;
#endif /* _INKERNEL_ */

#endif /* _TYPES_FERMIONSH_ */

