#ifndef _TESTINGH_
#define _TESTINGH_

#include <cstdlib>
#include <cstdio>

#include "hmcerrs.h"
#include "globaldefs.h"
#include "types.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_geometry.h"
#include "host_update_heatbath.h"
#include "host_fermionobservables.h"

void testing_correlator(hmc_gaugefield* gf, inputparameters* parameters);
void testing_fermionmatrix();
void testing_eoprec_spinor();
void print_fullspinorfield(hmc_spinor* in);

void testing_spinor();
void print_su3mat(hmc_su3matrix* A);
void print_staplemat(hmc_staplematrix* A);
void print_linkplusstaplematrix(hmc_gaugefield * in, int pos, int t, int mu);
void testing_su3mat();
void testing_gaugefield();
void testing_geometry();
void testing_su3matrix(hmc_gaugefield * in, int spacepos, int timepos);
void testing_adjoin(hmc_gaugefield * in, int spacepos, int timepos);
void testing_det_global(hmc_gaugefield * in);
void testing_matrix_ops(hmc_gaugefield * in);

void testing_heatbath_norandommat_no123(hmc_su3matrix * su3_in, hmc_staplematrix * staple_in, hmc_su3matrix * out, hmc_float beta);

void testing_heatbath_no123(hmc_su3matrix * su3_in, hmc_staplematrix * staple_in, hmc_su3matrix * out, hmc_float beta, hmc_float * rnd_array, int * cter_out);

void testing_heatbath(hmc_su3matrix * su3_in, hmc_staplematrix * staple_in, hmc_su3matrix * out, hmc_float beta, hmc_float * rnd_array, int * cter_out);

#endif
