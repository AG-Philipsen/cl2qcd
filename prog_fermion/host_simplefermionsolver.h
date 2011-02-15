#ifndef _SIMPLEFERMIONSOLVERH_
#define _SIMPLEFERMIONSOLVERH_

#include <cstdlib>
#include <cstring>
#include <iostream>

#include "host_operations.h"
#include "host_geometry.h"
#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"

extern int const iter_refresh;
extern double const epssquare;

hmc_error simple_correlator(hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, int cgmax);

hmc_error simple_solver(hmc_spinor_field* in, hmc_spinor_field* out, hmc_spinor_field* b, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, int cgmax);

hmc_error simple_solve_bicgstab(hmc_spinor_field* inout,hmc_spinor_field* source,hmc_gaugefield* gaugefield,hmc_float kappa,hmc_float mu,int cgmax);
hmc_error create_point_source(hmc_spinor_field* b,int i,int spacepos,int timepos,hmc_float kappa, hmc_float mu, hmc_gaugefield* gaugefield);
hmc_error M(hmc_spinor_field* in, hmc_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu);
hmc_error dslash(hmc_spinor_field* in, hmc_spinor_field* out, hmc_gaugefield* gaugefield);


#endif
