#ifndef _TYPESH_
#define _TYPESH_

#include "globaldefs.h"

//typedef double hmc_float;
typedef float hmc_float;

#ifdef _INKERNEL_
__constant hmc_float hmc_one_f = 1.0f;
#else
hmc_float const hmc_one_f = static_cast<hmc_float>(1);
#endif

struct hmc_complex {
  hmc_float re;
  hmc_float im;
};

#ifdef _INKERNEL_
__constant struct hmc_complex hmc_complex_one={1., 0.};
__constant struct hmc_complex hmc_complex_zero = {0., 0.};
__constant struct hmc_complex hmc_complex_i = {0., 1.};
#else
hmc_complex const hmc_complex_one = {1., 0.};
hmc_complex const hmc_complex_zero = {0., 0.};
hmc_complex const hmc_complex_i = {0., 1.};
#endif

//define a spinor field:  spinor[spin-color][coord3d][coord_time]
#ifdef _INKERNEL_
typedef struct hmc_complex hmc_full_spinor [NSPIN*NC];
typedef struct hmc_complex hmc_full_spinor_field [NSPIN*NC][VOLSPACE][NTIME];
#else
typedef hmc_complex hmc_full_spinor [NSPIN*NC];
typedef hmc_complex hmc_full_spinor_field [NSPIN*NC][VOLSPACE][NTIME];
#endif

#ifndef _INKERNEL_

//define a gauge field: gauge[su3][mu][coord3d][coord_time]
#ifdef _RECONSTRUCT_TWELVE_
typedef hmc_complex hmc_su3matrix [NC*(NC-1)];
typedef hmc_complex hmc_staplematrix [NC*NC];
typedef hmc_complex hmc_gaugefield [NC*(NC-1)][NDIM][VOLSPACE][NTIME];
#else
typedef hmc_complex hmc_su3matrix [NC][NC];
typedef hmc_su3matrix hmc_staplematrix;
typedef hmc_complex hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME];
#endif

#endif

typedef hmc_float hmc_ocl_su3matrix;
typedef hmc_float hmc_ocl_gaugefield;


#endif
