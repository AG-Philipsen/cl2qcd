#ifndef _OPERATIONS_GAUGEFIELDH_
#define _OPERATIONS_GAUGEFIELDH_
#include <iostream>
#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"
#include "host_geometry.h"
#include "host_operations_complex.h"
#include <cmath>

//gaugefield and su3 operations, global and local
hmc_error adjoin_su3matrix(hmc_su3matrix * mat);
hmc_complex trace_su3matrix(hmc_su3matrix * mat);
hmc_complex det_su3matrix(hmc_su3matrix * U);
hmc_error copy_su3matrix(hmc_su3matrix *out, hmc_su3matrix *in); 
hmc_error copy_staplematrix(hmc_staplematrix *out, hmc_staplematrix *in);
hmc_error unit_su3matrix(hmc_su3matrix * u); 
hmc_error random_su3matrix(hmc_su3matrix * u);
hmc_error zero_su3matrix(hmc_su3matrix * u); 
hmc_error zero_staplematrix(hmc_staplematrix * u);
hmc_error accumulate_su3matrix_prod(hmc_su3matrix *acc, hmc_su3matrix *multiplicator);
hmc_error multiply_su3matrices(hmc_su3matrix *out, hmc_su3matrix *p, hmc_su3matrix *q);
hmc_error multiply_staplematrix(hmc_staplematrix *out, hmc_su3matrix *p, hmc_staplematrix *q);
#ifdef _RECONSTRUCT_TWELVE_
hmc_complex reconstruct_su3(hmc_su3matrix *in, int ncomp);
#endif

hmc_error set_gaugefield_cold(hmc_gaugefield * field);
hmc_error set_gaugefield_hot(hmc_gaugefield * field);

hmc_error copy_gaugefield_from_ildg_format(hmc_gaugefield * gaugefield, hmc_float * gaugefield_tmp, int check);
hmc_error copy_gaugefield_to_ildg_format(ildg_gaugefield * dest, hmc_gaugefield * source);
hmc_error copy_to_ocl_format(hmc_ocl_gaugefield* host_gaugefield,hmc_gaugefield* gaugefield);
hmc_error copy_from_ocl_format(hmc_gaugefield* gaugefield,hmc_ocl_gaugefield* host_gaugefield);

hmc_error get_su3matrix(hmc_su3matrix* out, hmc_gaugefield * in, int spacepos, int timepos, int mu); //cl
hmc_error put_su3matrix(hmc_gaugefield * field, hmc_su3matrix * in, int spacepos, int timepos, int mu); //cl

hmc_error adjoin_su3(hmc_gaugefield * in, hmc_gaugefield * out);
hmc_complex global_trace_su3(hmc_gaugefield * field, int mu);
hmc_error accumulate_su3matrices_add(hmc_staplematrix *p, hmc_su3matrix *q);
void reduction (hmc_complex dest[su2_entries], hmc_staplematrix src, const int rand);
void extend (hmc_su3matrix * dest, const int random, hmc_complex src[su2_entries]);
hmc_error project_su3(hmc_su3matrix *U);

void gaugefield_apply_bc(hmc_su3matrix * in, hmc_float theta);
void gaugefield_apply_chem_pot(hmc_su3matrix * u, hmc_su3matrix * udagger, hmc_float chem_pot_re, hmc_float chem_pot_im);

void local_polyakov(hmc_gaugefield * field, hmc_su3matrix * prod, int n);
void local_plaquette(hmc_gaugefield * field, hmc_su3matrix * prod, int n, int t, int mu, int nu );

// copy-functions within cpu memory, gaugefield-related layers
hmc_error copy_gaugefield(hmc_gaugefield * source, hmc_gaugefield * dest);
hmc_error copy_gaugemomenta(hmc_gauge_momentum * source, hmc_gauge_momentum * dest);

#endif
