/** @file
 * Common types used by all program parts, both on host and OpenCL device.
 */

#ifndef _TYPESH_
#define _TYPESH_

#include "globaldefs.h"
#ifndef _INKERNEL_
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif
#endif

#ifdef _INKERNEL_
#define CONST __constant
#else /* _INKERNEL_ */
#define CONST const
#endif /* _INKERNEL_ */

/** The floating precision type used by hmc, can be 32 or 64 bit. */
#ifdef _USEDOUBLEPREC_
typedef double hmc_float __attribute__((aligned(8)));
#else
typedef float hmc_float __attribute__((aligned(4)));
#endif

/** Complex number type, precision is the same as for hmc_float */
#ifdef _INKERNEL_
typedef struct {
	hmc_float re;
	hmc_float im;
} hmc_complex
#ifdef _USEDOUBLEPREC_
__attribute__((aligned(16)));
#else
__attribute__((aligned(8)));
#endif
typedef struct {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e10;
	hmc_complex e11;
} Matrixsu2;
typedef struct {
	hmc_float e00;
	hmc_float e01;
	hmc_float e10;
	hmc_float e11;
} Matrixsu2_pauli;
#else
struct hmc_complex {
	hmc_float re;
	hmc_float im;
};
struct Matrixsu2 {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e10;
	hmc_complex e11;
} ;
struct Matrixsu2_pauli {
	hmc_float e00;
	hmc_float e01;
	hmc_float e10;
	hmc_float e11;
} ;
#endif

#ifdef _INKERNEL_
__constant hmc_complex hmc_complex_one = {1., 0.};
__constant hmc_complex hmc_complex_zero = {0., 0.};
__constant hmc_complex hmc_complex_minusone = { -1., 0.};
__constant hmc_complex hmc_complex_i = {0., 1.};
#else
/** A complex 1 */
hmc_complex const hmc_complex_one = {1., 0.};
/** A complex -1 */
hmc_complex const hmc_complex_minusone = { -1., 0.};
/** A complex 0 */
hmc_complex const hmc_complex_zero = {0., 0.};
/** A complex i */
hmc_complex const hmc_complex_i = {0., 1.};
#endif

//matrix definitions
#ifdef _INKERNEL_
//a generic 3x3 matrix
typedef struct {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e02;
	hmc_complex e10;
	hmc_complex e11;
	hmc_complex e12;
	hmc_complex e20;
	hmc_complex e21;
	hmc_complex e22;
} Matrix3x3;
//an su3 matrix
typedef struct {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e02;
	hmc_complex e10;
	hmc_complex e11;
	hmc_complex e12;
	hmc_complex e20;
	hmc_complex e21;
	hmc_complex e22;
} Matrixsu3;
#else
//an su3 matrix_
struct Matrixsu3 {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e02;
	hmc_complex e10;
	hmc_complex e11;
	hmc_complex e12;
	hmc_complex e20;
	hmc_complex e21;
	hmc_complex e22;
};
#endif // ifdef _INKERNEL_

#ifdef _USE_SOA_
typedef hmc_complex Matrixsu3StorageType;
#else
typedef Matrixsu3 Matrixsu3StorageType;
#endif

#ifndef _INKERNEL_
/** Storage type for state of the random number generator */
typedef cl_ulong4 hmc_ocl_ran;
#endif /* _INKERNEL_ */

//CP: this an algebraelement
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

#ifdef _USE_SOA_
typedef hmc_float aeStorageType;
#else
typedef ae aeStorageType;
#endif

#endif /* _TYPESH_ */
