#ifndef _OPERATIONSH_
#define _OPERATIONSH_
#include "globals.h"
#include "types.h"
#include "hmcerrs.h"

//local operations
hmc_complex complexconj(hmc_complex *in);
hmc_complex complexmult(hmc_complex *a, hmc_complex *b);
hmc_complex complexadd(hmc_complex *a, hmc_complex *b);
hmc_complex complexsubtract(hmc_complex *a, hmc_complex *b);
hmc_error complexaccumulate(hmc_complex *inout, hmc_complex *incr);

hmc_error adjoin_su3matrix(hmc_su3matrix * mat);
hmc_complex trace_su3matrix(hmc_su3matrix * mat);
hmc_error copy_su3matrix(hmc_su3matrix *out, hmc_su3matrix *in);
hmc_error unit_su3matrix(hmc_su3matrix * u);
hmc_error accumulate_su3matrix_prod(hmc_su3matrix *acc, hmc_su3matrix *multiplicator);
hmc_error multiply_su3matrices(hmc_su3matrix *out, hmc_su3matrix *p, hmc_su3matrix *q);
#ifdef _RECONSTRUCT_TWELVE_
hmc_complex reconstruct_su3(hmc_su3matrix *in, int ncomp);
#endif

//spinor operations, global and local
hmc_error set_zero_spinor(hmc_full_spinor_field * field);
hmc_float local_squarenorm(hmc_full_spinor_field *field, int spacepos, int timepos);
hmc_float global_squarenorm(hmc_full_spinor_field *field);
hmc_error fill_with_one(hmc_full_spinor_field *field, int spacepos, int timepos, int j);


//gaugefield and su3 operations, global and local
hmc_error set_gaugefield_cold(hmc_gaugefield * field);

hmc_error get_su3matrix(hmc_su3matrix* out, hmc_gaugefield * in, int spacepos, int timepos, int mu);
hmc_error put_su3matrix(hmc_gaugefield * field, hmc_su3matrix * in, int spacepos, int timepos, int mu);

hmc_error adjoin_su3(hmc_gaugefield * in, hmc_gaugefield * out);
hmc_complex global_trace_su3(hmc_gaugefield * field, int mu);


#endif
