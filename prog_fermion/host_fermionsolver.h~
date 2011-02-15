#ifndef _FERMIONSOLVERH_
#define _FERMIONSOLVERH_

#include <cstdlib>
#include <cstring>
#include <iostream>

#include "operations.h"
#include "geometry.h"
#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"

extern double const epssquare;
extern int const iter_refresh;

hmc_error solver(hmc_spinor_field* in, hmc_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, int cgmax);

hmc_error solve_bicgstab(hmc_eoprec_spinor_field* out,hmc_eoprec_spinor_field* in,hmc_eoprec_spinor_field* source,hmc_gaugefield* gaugefield,hmc_float kappa,hmc_float mu,int cgmax);
hmc_error create_point_source(hmc_eoprec_spinor_field* be,hmc_eoprec_spinor_field* bo,int i,int spacepos,int timepos,hmc_float kappa, hmc_float mu, hmc_gaugefield* gaugefield);
hmc_error Aee(hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu);
hmc_error dslash(hmc_eoprec_spinor_field* out, hmc_eoprec_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float kappa, int evenodd);
hmc_error sitediagonal(hmc_eoprec_spinor_field* out, hmc_eoprec_spinor_field* in, hmc_float kappa, hmc_float mu);
hmc_error inverse_sitediagonal(hmc_eoprec_spinor_field* out, hmc_eoprec_spinor_field* in, hmc_float kappa, hmc_float mu);

hmc_error print_spinor(hmc_spinor* in);

#endif
