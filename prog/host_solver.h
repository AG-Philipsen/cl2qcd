//!!CP: LZ should update this...

#ifndef _SOLVERH_
#define _SOLVERH_

#include <cstdlib>
#include <cstring>
#include <iostream>

#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"
#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_operations_spinorfield.h"
#include "host_operations_fermionmatrix.h"

//extern int const iter_refresh;
extern int const iter_refresh;
extern hmc_float const epssquare;

extern int const use_eo;

hmc_error solver(hmc_spinor_field* in, hmc_spinor_field* out, hmc_spinor_field* b, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax);

hmc_error solver(hmc_spinor_field* in, hmc_spinor_field* out, hmc_eoprec_spinor_field* be, hmc_eoprec_spinor_field* bo, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax);

hmc_error bicgstab(hmc_spinor_field* inout,hmc_spinor_field* source,hmc_gaugefield* gaugefield,hmc_float kappa,hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax);
hmc_error bicgstab_eoprec(hmc_eoprec_spinor_field* inout,hmc_eoprec_spinor_field* source,hmc_gaugefield* gaugefield,hmc_float kappa,hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im,int cgmax);

hmc_error cg(hmc_spinor_field* inout, hmc_eoprec_spinor_field* source, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax);
hmc_error cg_eoprec(hmc_eoprec_spinor_field* inout, hmc_eoprec_spinor_field* source, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax);

#endif
