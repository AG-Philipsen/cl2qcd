#ifndef _GAUGEOBSERVABLESH_
#define _GAUGEOBSERVABLESH_
#include "globals.h"
#include "types.h"
#include "hmcerrs.h"
#include "operations.h"
#include "geometry.h"

hmc_float plaquette(hmc_gaugefield * field, hmc_float* tplaq, hmc_float* splaq);
hmc_float plaquette(hmc_gaugefield * field);
hmc_complex polyakov(hmc_gaugefield * field);
// the spatial polyakov loop functions have not been checked
// spatial_polyakov(u,1) gives different result than polyakov_x(u)
hmc_complex spatial_polyakov(hmc_gaugefield * field, int dir);
hmc_complex polyakov_x(hmc_gaugefield * field);
hmc_complex polyakov_y(hmc_gaugefield * field);
hmc_complex polyakov_z(hmc_gaugefield * field);

#endif
