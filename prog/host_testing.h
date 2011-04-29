/** @file
 * A bunch of functions useful for testing and debugging.
 */
#ifndef _TESTINGH_
#define _TESTINGH_

#include <cstdlib>
#include <cstdio>

#include "hmcerrs.h"
#include "globaldefs.h"
#include "types.h"
#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_operations_spinorfield.h"
#include "host_operations_fermionmatrix.h"
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

void testing_heatbath_norandommat_no123(hmc_su3matrix * su3_in, hmc_staplematrix * staple_in, hmc_su3matrix * out);

void testing_heatbath_no123(hmc_su3matrix * su3_in, hmc_staplematrix * staple_in, hmc_su3matrix * out, hmc_float beta, hmc_float * rnd_array, int * cter_out);

void testing_heatbath(hmc_su3matrix * su3_in, hmc_staplematrix * staple_in, hmc_su3matrix * out, hmc_float beta, hmc_float * rnd_array, int * cter_out);

void testing_colorvector_ops();
void testing_matrix_spinor_ops();
void testing_matrix_spinor_functions();
void testing_fermionmatrix_functions();

void unit_spinor(hmc_spinor * in);
void i_spinor(hmc_spinor * in);
void print_spinor(hmc_spinor * in);
void set_comp_to_one_spinor(hmc_spinor * in, int comp);
void set_comp_to_i_spinor(hmc_spinor * in, int comp);
void print_colorvector(hmc_color_vector * in);
void unit_colorvector(hmc_color_vector * in);
void zero_colorvector(hmc_color_vector * in);
void  i_colorvector(hmc_color_vector * in);
void acc_colorvector(hmc_color_vector * inout, hmc_color_vector * incr);
void mult_colorvector(hmc_color_vector * inout, hmc_complex * factor);
void fill_su3matrix_one(hmc_su3matrix * u);
void fill_su3matrix_i(hmc_su3matrix * u);

#endif
