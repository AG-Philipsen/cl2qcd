#ifndef _KAPPAH_
#define _KAPPAH_

#include <cstdlib>
#include <cstdio>

#include "hmcerrs.h"
#include "globaldefs.h"
#include "types.h"
#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_operations_matrix.h"
#include "host_operations_gaugefield.h"

void kappa_karsch (hmc_gaugefield* field, hmc_float & kappa, const hmc_float beta);

void local_Q_plaquette(hmc_gaugefield * field, hmc_staplematrix * prod, int n, int t, int mu, int nu );

void testing_Qplak (hmc_gaugefield * field, hmc_float* plaq, hmc_float* tplaq, hmc_float* splaq);

void kappa_clover (hmc_gaugefield* field, hmc_float & kappa, const hmc_float beta);



#endif
