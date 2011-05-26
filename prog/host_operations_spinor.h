/** @file
 * Mathematical operations on spinors.
 * A spinor is a NDIM*NC-component complex vector, living in Color- and Diracspace.
 * /f[
 * \psi = \ldots
 * /f]
 * Each of the NDIM NC-component complex 'subvector' is called a Colorvector, living in SU(NC)-Colorspace.
 * (Note: It is assumed throughout the program that NC = 3 and NDIM = 4)
 */

#ifndef _OPERATIONS_SPINORH_
#define _OPERATIONS_SPINORH_

#include <iostream>
#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"
#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include <cmath>

/**
 * Mutliply a SU(3)-Colorvector with a SU(3)-Matrix.
 *
 * @param[in] u The SU(3)-matrix
 * @param[in] in The Colorvector to be multiplied
 * @param[in] out The Colorvector the result is written to
 * @return Error code as defined in hmcerrs.h
 * @remark tested by CP
 */
hmc_error su3matrix_times_colorvector(hmc_su3matrix* u, hmc_color_vector* in, hmc_color_vector* out);//CP: checked

/**
 * Set a specific spinor to zero.
 *
 * @param[in,out] inout The spinor to be zeroed
 * @return Error code as defined in hmcerrs.h
 * @remark tested by CP
 */
hmc_error set_local_zero_spinor(hmc_spinor* inout);//CP: checked

/**
 * Calculate the squarenorm of a spinor.
 *
 * @param[in] in A spinor
 * @return The squarenorm (hmc_float)
 * @remark tested by CP
 */
hmc_float spinor_squarenorm(hmc_spinor* in);//CP: checked

/**
 * Multiply a spinor with a hmc_float.
 *
 * @param[in,out] inout Spinor to be multiplied
 * @param[in] factor A hmc_float
 * @return Error code as defined in hmcerrs.h
 * @remark tested by CP
 */
hmc_error real_multiply_spinor(hmc_spinor* inout, hmc_float factor);//CP: checked, corrected

/**
 * Accumulate one spinor to another (componentwise).
 *
 * @param[in,out] inout Spinor to be incremented
 * @param[in] incr Spinor to incremented
 * @return Error code as defined in hmcerrs.h
 * @remark tested by CP
 */
hmc_error spinors_accumulate(hmc_spinor* inout, hmc_spinor* incr);//CP: checked

/**
 * Multiplies each Colorvector of a Spinor with a SU(3)-Matrix
 *
 * @param[in] u A SU(3)-Matrix
 * @param[in] in A Spinor
 * @param[out] out The spinor to write the result to
 * @return Error code as defined in hmcerrs.h
 * @remark tested by CP
 */
hmc_error su3matrix_times_spinor(hmc_su3matrix* u, hmc_spinor* in, hmc_spinor* out);//CP: checked, corrected

/**
 * Multiplies each of the NDIM-vectors of a spinor (of which there are NC) by the Dirac-Matrix gamma_0.
 *
 * @param[in] in Input spinor
 * @param[out] out Output spinor
 * @return Error code as defined in hmcerrs.h
 * @remark tested by CP
 */
hmc_error multiply_spinor_gamma0(hmc_spinor* in,hmc_spinor* out);//CP: checked explicitly

/**
 * Multiplies each of the NDIM-vectors of a spinor (of which there are NC) by the Dirac-Matrix gamma_1.
 *
 * @param[in] in Input spinor
 * @param[out] out Output spinor
 * @return Error code as defined in hmcerrs.h
 * @remark tested by CP
 */
hmc_error multiply_spinor_gamma1(hmc_spinor* in,hmc_spinor* out);//CP: checked explicitly

/**
 * Multiplies each of the NDIM-vectors of a spinor (of which there are NC) by the Dirac-Matrix gamma_2.
 *
 * @param[in] in Input spinor
 * @param[out] out Output spinor
 * @return Error code as defined in hmcerrs.h
 * @remark tested by CP
 */
hmc_error multiply_spinor_gamma2(hmc_spinor* in,hmc_spinor* out);//CP: checked explicitly

/**
 * Multiplies each of the NDIM-vectors of a spinor (of which there are NC) by the Dirac-Matrix gamma_3.
 *
 * @param[in] in Input spinor
 * @param[out] out Output spinor
 * @return Error code as defined in hmcerrs.h
 * @remark tested by CP
 */
hmc_error multiply_spinor_gamma3(hmc_spinor* in,hmc_spinor* out);//CP: checked explicitly

/**
 * Mutliplies a spinor with i*factor and a gamma_5 in Diracspace.
 *
 * @param[in] in Input Spinor
 * @param[out] out Output Spinor
 * @param[in] factor Multiplying factor
 * @return Error code as defined in hmcerrs.h
 * @remark tested by CP
 */
hmc_error multiply_spinor_i_factor_gamma5(hmc_spinor* in,hmc_spinor* out, hmc_float factor);//CP: checked explicitly

/**
 * Calculates the expression
 * /f[
 * U \left( 1 \pm \gamma_0 \right) \psi \;, 
 * /f]
 * where /f$\psi/f$ is a spinor and /f$U/f$ is a SU(3)-Matrix. This expression allows to save some calculations ('spinprojection').
 * @param[in] u SU(3)-Matrix
 * @param[in,out] spin Input/Output spinor
 * @param[in] sign /f$\pm/f$ one
 * @return Error code as defined in hmcerrs.h
 * @remark tested by CP
 */
hmc_error spinprojectproduct_gamma0(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign);//CP: checked explicitly

/**
 * Calculates the expression
 * /f[
 * U \left( 1 \pm \gamma_1 \right) \psi \;, 
 * /f]
 * where /f$\psi/f$ is a spinor and /f$U/f$ is a SU(3)-Matrix. This expression allows to save some calculations ('spinprojection').
 * @param[in] u SU(3)-Matrix
 * @param[in,out] spin Input/Output spinor
 * @param[in] sign /f$\pm/f$ one
 * @return Error code as defined in hmcerrs.h
 * @remark tested by CP
 */
hmc_error spinprojectproduct_gamma1(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign);//CP: checked explicitly

/**
 * Calculates the expression
 * /f[
 * U \left( 1 \pm \gamma_2 \right) \psi \;, 
 * /f]
 * where /f$\psi/f$ is a spinor and /f$U/f$ is a SU(3)-Matrix. This expression allows to save some calculations ('spinprojection').
 * @param[in] u SU(3)-Matrix
 * @param[in,out] spin Input/Output spinor
 * @param[in] sign /f$\pm/f$ one
 * @return Error code as defined in hmcerrs.h
 * @remark tested by CP
 */
hmc_error spinprojectproduct_gamma2(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign);//CP: checked explicitly

/**
 * Calculates the expression
 * /f[
 * U \left( 1 \pm \gamma_3 \right) \psi \;, 
 * /f]
 * where /f$\psi/f$ is a spinor and /f$U/f$ is a SU(3)-Matrix. This expression allows to save some calculations ('spinprojection').
 * @param[in] u SU(3)-Matrix
 * @param[in,out] spin Input/Output spinor
 * @param[in] sign /f$\pm/f$ one
 * @return Error code as defined in hmcerrs.h
 * @remark tested by CP
 */
hmc_error spinprojectproduct_gamma3(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign);//CP: checked explicitly

/**
 * Apply Boundary Condition to a spinor. 
 * This corresponds to multiplying the spinor by a (complex) factor of /f$\exp(i*\theta) /f$.
 * /f$\theta = 0/f$ are 'periodic' BC.
 * /f$\theta = \Pi/f$ are 'antiperiodic' BC.
 *
 * @param[in,out] in Spinor to be changed
 * @param[in] theta angle /f$ \theta /f$
 * @return void
 * @remark tested by CP
 * @todo the calculation involves sin- and cos-evaluations. Perhaps one should optimize this for the two special cases mentioned above.
 */
void spinor_apply_bc(hmc_spinor * in, hmc_float theta);//CP: checked explicitly

/**
 * Applies the diagonal part of /f$ M^{\dagger}M/f$ to a spinor, with /f$M/f$ being the fermionmatrix.
 *
 * @param[in,out] spininout Spinor to be changed
 * @param[in] mubar /f$2\kappa\mu/f$ of twisted-mass fermions.
 * @return void
 * @remark tested by CP
 * @todo at the moment, this is only twisted-mass!!
 */
void M_diag_local(hmc_spinor* spininout, hmc_float mubar);//CP: checked explicitly

/**
 * Calculates the /f$\mu=0/f$-component of /f$ \notD/f$ to a spinor.
 *
 * @param[in] spinnext Neighbouring spinor in positive direction
 * @param[in] spinprev Neighbouring spinor in negative direction
 * @param[in] spinout spinor to add contribution to
 * @param[in] u SU(3)-Matrix linking to spinnext
 * @param[in] udagger (adjoined) SU(3)-Matrix linking to spinprev
 * @return void
 * @remark tested by CP
 * @todo at the moment, this is only twisted-mass!!
 */
void dslash_0(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger);//CP: checked explicitly

/**
 * Calculates the /f$\mu=1/f$-component of /f$ \notD/f$ to a spinor.
 *
 * @param[in] spinnext Neighbouring spinor in positive direction
 * @param[in] spinprev Neighbouring spinor in negative direction
 * @param[in] spinout spinor to add contribution to
 * @param[in] u SU(3)-Matrix linking to spinnext
 * @param[in] udagger (adjoined) SU(3)-Matrix linking to spinprev
 * @return void
 * @remark tested by CP
 * @todo at the moment, this is only twisted-mass!!
 */
void dslash_1(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger);//CP: checked explicitly

/**
 * Calculates the /f$\mu=2/f$-component of /f$ \notD/f$ to a spinor.
 *
 * @param[in] spinnext Neighbouring spinor in positive direction
 * @param[in] spinprev Neighbouring spinor in negative direction
 * @param[in] spinout spinor to add contribution to
 * @param[in] u SU(3)-Matrix linking to spinnext
 * @param[in] udagger (adjoined) SU(3)-Matrix linking to spinprev
 * @return void
 * @remark tested by CP
 * @todo at the moment, this is only twisted-mass!!
 */
void dslash_2(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger);//CP: checked explicitly

/**
 * Calculates the /f$\mu=3/f$-component of /f$ \notD/f$ to a spinor.
 *
 * @param[in] spinnext Neighbouring spinor in positive direction
 * @param[in] spinprev Neighbouring spinor in negative direction
 * @param[in] spinout spinor to add contribution to
 * @param[in] u SU(3)-Matrix linking to spinnext
 * @param[in] udagger (adjoined) SU(3)-Matrix linking to spinprev
 * @return void
 * @remark tested by CP
 * @todo at the moment, this is only twisted-mass!!
 */
void dslash_3(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger);//CP: checked explicitly

/**
 * Applies /$f \gamma_5 f$/ to a spinor
 * 
 * @param[in] inout spinor to be multiplied by the Dirac-Matrix
 */
void gamma_5_spinor(hmc_full_spinor inout);

void su3_vector_times_minusone(hmc_su3vector inout);

void su3_vector_acc(hmc_su3vector in, hmc_su3vector out);

void su3_vector_multiple_add(hmc_su3vector in1, hmc_su3vector in2, hmc_su3vector out);

void spinproj_gamma0_a(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign);

void spinproj_gamma0_b(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign);

//calculates the trace of generator times 3x3-matrix and stores this in a su3-algebraelement
void tr_lambda_u(hmc_3x3matrix in, hmc_algebraelement out);

//calculates the Dirac-Trace of the matrix resulting from multiplying v*u^dagger + w*x^dagger, where u, v, w, x are SU(3)-vectors
//	the result is a 3x3-matrix
void tr_v_times_u_dagger(hmc_su3vector v, hmc_su3vector u, hmc_su3vector w, hmc_su3vector x, hmc_3x3matrix out);

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// deprecated functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//CP: this was usefull because M_daggerM can also be calculated by successivly using M and gamma_5

/**
 * Applies the diagonal part of /f$ M^{\dagger}M/f$ to a spinor, with /f$M/f$ being the fermionmatrix.
 *
 * @param[in,out] spininout Spinor to be changed
 * @param[in] mubar /f$2\kappa\mu/f$ of twisted-mass fermions.
 * @return void
 * @remark tested by CP
 * @todo at the moment, this is only twisted-mass!!
 * @remark deprecated by CP
 */
void MdaggerM_diag_local(hmc_spinor* spininout, hmc_float kappa, hmc_float mu);//CP: not checked

/**
 * Calculates the /f$\mu=0/f$-component of /f$ \notD^{\dagger}/f$ to a spinor.
 *
 * @param[in] spinnext Neighbouring spinor in positive direction
 * @param[in] spinprev Neighbouring spinor in negative direction
 * @param[in] spinout spinor to add contribution to
 * @param[in] u SU(3)-Matrix linking to spinnext
 * @param[in] udagger (adjoined) SU(3)-Matrix linking to spinprev
 * @return void
 * @remark tested by CP
 * @remark deprecated by CP
 * @todo at the moment, this is only twisted-mass!!
 */
void ddaggerslash_0(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger);//CP: not checked

/**
 * Calculates the /f$\mu=1/f$-component of /f$ \notD^{\dagger}/f$ to a spinor.
 *
 * @param[in] spinnext Neighbouring spinor in positive direction
 * @param[in] spinprev Neighbouring spinor in negative direction
 * @param[in] spinout spinor to add contribution to
 * @param[in] u SU(3)-Matrix linking to spinnext
 * @param[in] udagger (adjoined) SU(3)-Matrix linking to spinprev
 * @return void
 * @remark tested by CP
 * @remark deprecated by CP
 * @todo at the moment, this is only twisted-mass!!
 */
void ddaggerslash_1(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger);//CP: not checked

/**
 * Calculates the /f$\mu=2/f$-component of /f$ \notD^{\dagger}/f$ to a spinor.
 *
 * @param[in] spinnext Neighbouring spinor in positive direction
 * @param[in] spinprev Neighbouring spinor in negative direction
 * @param[in] spinout spinor to add contribution to
 * @param[in] u SU(3)-Matrix linking to spinnext
 * @param[in] udagger (adjoined) SU(3)-Matrix linking to spinprev
 * @return void
 * @remark tested by CP
 * @remark deprecated by CP
 * @todo at the moment, this is only twisted-mass!!
 */
void ddaggerslash_2(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger);//CP: notchecked

/**
 * Calculates the /f$\mu=3/f$-component of /f$ \notD^{\dagger}/f$ to a spinor.
 *
 * @param[in] spinnext Neighbouring spinor in positive direction
 * @param[in] spinprev Neighbouring spinor in negative direction
 * @param[in] spinout spinor to add contribution to
 * @param[in] u SU(3)-Matrix linking to spinnext
 * @param[in] udagger (adjoined) SU(3)-Matrix linking to spinprev
 * @return void
 * @remark tested by CP
 * @remark deprecated by CP
 * @todo at the moment, this is only twisted-mass!!
 */
void ddaggerslash_3(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger);//CP: not checked


#endif /* _OPERATIONS_SPINORH_ */


