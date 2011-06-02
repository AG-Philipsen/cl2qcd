/** @file
 * Common types used by HMC, both on host and OpenCL device.
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
typedef double hmc_float;
#else
typedef float hmc_float;
#endif

/** An OpenCL-compatible constant 1.
 * @todo this is rediculous, 1.0(f) is way more readable, let the
 *       compiler do the optimizations.
 */
#ifdef _INKERNEL_
__constant hmc_float hmc_one_f = 1.0f;
#else
hmc_float const hmc_one_f = static_cast<hmc_float>(1);
#endif

/** Complex number type, precision is the same as for hmc_float */
#ifdef _INKERNEL_
typedef struct {
  hmc_float re;
  hmc_float im;
} hmc_complex;
#else
struct hmc_complex {
  hmc_float re;
  hmc_float im;
};
#endif

#ifdef _INKERNEL_
__constant hmc_complex hmc_complex_one={1., 0.};
__constant hmc_complex hmc_complex_zero = {0., 0.};
__constant hmc_complex hmc_complex_minusone = {-1., 0.};
__constant hmc_complex hmc_complex_i = {0., 1.};
#else
/** A complex 1 */
hmc_complex const hmc_complex_one = {1., 0.};
/** A complex -1 */
hmc_complex const hmc_complex_minusone = {-1., 0.};
/** A complex 0 */
hmc_complex const hmc_complex_zero = {0., 0.};
/** A complex i */
hmc_complex const hmc_complex_i = {0., 1.};
#endif

#ifndef _INKERNEL_
//define a spinor field:  spinor[spin-color][coord3d][coord_time]
typedef hmc_complex hmc_full_spinor [NSPIN*NC];
typedef hmc_complex hmc_full_spinor_field [NSPIN*NC][VOLSPACE][NTIME];

//define a gauge field: gauge[su3][mu][coord3d][coord_time]
#ifdef _RECONSTRUCT_TWELVE_
typedef hmc_complex hmc_su3matrix [NC*(NC-1)];
typedef hmc_complex hmc_staplematrix [NC*NC];
typedef hmc_complex hmc_3x3matrix[3][3];
typedef hmc_complex hmc_gaugefield [NC*(NC-1)][NDIM][VOLSPACE][NTIME];
typedef hmc_float ildg_gaugefield[2*NC*(NC-1)*NDIM*VOLSPACE*NTIME];
typedef hmc_float hmc_gauge_momentum;
#else
/** A generic SU3 matrix */
typedef hmc_complex hmc_su3matrix [NC][NC];
/** A matrix representing a staple */
typedef hmc_su3matrix hmc_staplematrix;
/** A generic 3x3 matrix */
typedef hmc_complex hmc_3x3matrix[3][3];
/** The full gaugefield */
typedef hmc_complex hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME];
typedef hmc_float ildg_gaugefield[2*NC*NC*NDIM*VOLSPACE*NTIME];

#endif

typedef hmc_float hmc_gauge_momentum;
typedef hmc_float hmc_algebraelement [NC*NC-1];
typedef struct {
  hmc_float e0;
 	hmc_float e1;
  hmc_float e2;
  hmc_float e3;
  hmc_float e4;
  hmc_float e5;
  hmc_float e6;
  hmc_float e7;
} hmc_algebraelement2;

typedef hmc_complex hmc_su3vector[NC];

#endif // ifndef _INKERNEL_

typedef hmc_float hmc_ocl_spinor;
typedef hmc_complex hmc_ocl_su3matrix;
typedef hmc_complex hmc_ocl_staplematrix;
typedef hmc_float hmc_ocl_gaugefield;

//define a spinor field:  spinor_field[spin-color*coord3d*coord_time]
typedef hmc_complex hmc_color_vector;
typedef hmc_complex hmc_spinor;
typedef hmc_complex hmc_spinor_field;
typedef hmc_complex hmc_eoprec_spinor_field;

#ifdef _USEDOUBLEPREC_
hmc_float CONST projectioneps = 10.e-12;
int CONST iter_refresh = 50;
hmc_float CONST epssquare=1e-14;
int CONST use_eo = 0;
#else
hmc_float CONST projectioneps = 10.e-6;
int CONST iter_refresh = 50;
hmc_float CONST epssquare=1e-12;
int CONST use_eo = 0;
#endif

#ifndef _INKERNEL_
/**
 * Work-group size for OpenCL kernels.
 * @bug The proper work-group size is kernel dependent and cannot be
        specified globally.
 */
const size_t local_work_size  = NUMTHREADS;
/**
 * Global number of threads for OpenCL kernels.
 * @bug The proper number of threads size is kernel dependent and cannot be
        specified globally.
 */
const size_t global_work_size = NUMTHREADS;
#endif

#ifndef _INKERNEL_
/** Storage type for state of the random number generator */
typedef cl_ulong4 hmc_ocl_ran;
/** The array of random number generator states for usage by an OpenCL device
 * @todo dynamically size according to requirements by kernels / devices
 * @warning some kernel use NUMTHREADS threads, make sure this is always bigger!
 */
#ifdef _USE_GPU_
#define NUMRNDSTATES 5120
#else
#define NUMRNDSTATES 64
#endif
typedef hmc_ocl_ran hmc_rndarray[NUMRNDSTATES];
#endif /* _INKERNEL_ */

#endif /* _TYPESH_ */

