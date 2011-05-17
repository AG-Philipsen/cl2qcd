/** @file
 * CB and BiCG based solver
 *
 * @todo CP: LZ should update this...
 */

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

/** globally avaiable switch to choose between normal and eoprec
 * todo CP: is this needed? There is a flag for that..
 */
extern int const use_eo;

/**
 * solves the sparse matrix system $f/ M \psi = b /f$ for given matrix M, initial field $f/ \psi /f$, source b and returns a solution $f/ \psi_{out} /f$
 * @param[in] parameters provides parameters needed
 * @param[in] in initial spinorfield
 * @param[in] b source field
 * @param[in] gaugefield input gaugefield
 * @param[in] use_cg switch to choose between CG and BiCGStab solver
 * @param[out] out output spinorfield
 * @remark tested by CP
 */
hmc_error solver(hmc_spinor_field* in, hmc_spinor_field* out, hmc_spinor_field* b, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax, int use_cg);

/**
 * solves the sparse matrix system $f/ M \psi = b /f$ for given matrix M, initial field $f/ \psi /f$, source b and returns a solution $f/ \psi_{out} /f$. For speedup, the EVEN-ODD preconditiong is used
 * @param[in] parameters provides parameters needed
 * @param[in] in initial spinorfield
 * @param[in] be even source field
 * @param[in] bo odd source field
 * @param[in] gaugefield input gaugefield
 * @param[in] use_cg switch to choose between CG and BiCGStab solver
 * @param[out] out output spinorfield
 * @remark tested by CP
 */
hmc_error solver_eoprec(hmc_spinor_field* in, hmc_spinor_field* out, hmc_eoprec_spinor_field* be, hmc_eoprec_spinor_field* bo, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax);

/**
 * Applies the BiCGStab-algorithm
 * @param[in] parameters provides parameters needed
 * @param[in,out] in initial/output spinorfield
 * @param[in] be even source field
 * @param[in] bo odd source field
 * @param[in] gaugefield input gaugefield
 * @remark tested by CP
 */
hmc_error bicgstab(hmc_spinor_field* inout,hmc_spinor_field* source,hmc_gaugefield* gaugefield,hmc_float kappa,hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax);

/**
 * Applies the eoprec-BiCGStab-algorithm
 * @param[in] parameters provides parameters needed
 * @param[in,out] in initial/output spinorfield
 * @param[in] be even source field
 * @param[in] bo odd source field
 * @param[in] gaugefield input gaugefield
 * @remark tested by CP
 */
hmc_error bicgstab_eoprec(hmc_eoprec_spinor_field* inout,hmc_eoprec_spinor_field* source,hmc_gaugefield* gaugefield,hmc_float kappa,hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im,int cgmax);

/**
 * Applies the CG-algorithm
 * @param[in] parameters provides parameters needed
 * @param[in,out] in initial/output spinorfield
 * @param[in] be even source field
 * @param[in] bo odd source field
 * @param[in] gaugefield input gaugefield
 * @remark tested by CP
 */
hmc_error cg(hmc_spinor_field* inout, hmc_eoprec_spinor_field* source, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax);

/**
 * Applies the eoprec CG-algorithm
 * @param[in] parameters provides parameters needed
 * @param[in,out] in initial/output spinorfield
 * @param[in] be even source field
 * @param[in] bo odd source field
 * @param[in] gaugefield input gaugefield
 * @remark tested by CP
 */
hmc_error cg_eoprec(hmc_eoprec_spinor_field* inout, hmc_eoprec_spinor_field* source, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax);

#endif
