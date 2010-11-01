#ifndef _TYPESH_
#define _TYPESH_

#include "globals.h"

typedef double hmc_float;

hmc_float const hmc_one_f = static_cast<hmc_float>(1);

struct hmc_complex {
  hmc_float re;
  hmc_float im;
};

//define a spinor field:  spinor[spin-color][coord3d][coord_time]
typedef hmc_complex hmc_full_spinor [NSPIN*NC];
typedef hmc_complex hmc_full_spinor_field [NSPIN*NC][VOLSPACE][NTIME];

//define a gauge field: gauge[su3][mu][coord3d][coord_time]
typedef hmc_complex hmc_su3matrix [NC][NC];
typedef hmc_complex hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME];


#endif
