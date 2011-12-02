/** @file
 * Operations on SU3 matrices.
 *
 * @todo A lot of these operations use pointers and in-place operation
 *       without need. These might prevent compiler optimizations.
 */
#ifndef _OPERATIONS_SU3MATRIXH_
#define _OPERATIONS_SU3MATRIXH_
#include <iostream>
#include "globaldefs.h"
#include "types.h"
#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_operations_matrix.h"
#include <cmath>

/**
 * Replace the given matrix by it's adjoint.
 *
 * @param[in,out] mat The matrix to replace by it's adjoint
 */
void adjoin_su3matrix(hmc_su3matrix * mat);
/**
 * Calculate the trace of the given SU3 matrix.
 *
 * @param mat An SU3 matrix
 * @return The trace
 */
hmc_complex trace_su3matrix(hmc_su3matrix * mat);
/**
 * Calculate the determinant of the given SU3 matrix.
 *
 * @param mat An SU3 matrix
 * @return The determinant
 */
hmc_complex det_su3matrix(hmc_su3matrix * U);
/**
 * Copy the contents of an SU3 matrix to another one.
 *
 * @param[out] The SU3 matrix to copy to
 * @param[in] The SU3 matrix to copy from
 */
void copy_su3matrix(hmc_su3matrix *out, hmc_su3matrix *in);
/**
 * Copy the contents of an staple matrix to another one.
 *
 * @param[out] The staple matrix to copy to
 * @param[in] The staple matrix to copy from
 */
void copy_staplematrix(hmc_staplematrix *out, hmc_staplematrix *in);
/**
 * Replace the given matrix by a unit matrix.
 *
 * @param[out] mat The matrix to replace by a unit one
 */
void unit_su3matrix(hmc_su3matrix * u);

Matrixsu3 unit_matrixsu3();

/**
 * Replace the given matrix by a random matrix.
 *
 * @param[out] mat The matrix to replace by a random one
 */
void random_su3matrix(hmc_su3matrix * u);
/**
 * Replace the given matrix by a zero matrix.
 *
 * @param[out] mat The matrix to replace by a zero one
 */
void zero_su3matrix(hmc_su3matrix * u);
/**
 * Replace the given matrix by a zero matrix.
 *
 * @param[out] mat The matrix to replace by a zero one
 */
void zero_staplematrix(hmc_staplematrix * u);
/**
 * Multiply two SU3 matrices in place
 *
 * @param[in,out] acc The matrix to be replaced by the multiplication result
 * @param[in] multiplicator The matrix to multiply witht the acc matrix
 */
void accumulate_su3matrix_prod(hmc_su3matrix *acc, hmc_su3matrix *multiplicator);
/**
 * Multiply two SU3 matrices
 *
 * @param[out] out The matrix into which to store the multiplication result
 * @param[in] p Left matrix for the multiplication
 * @param[in] q Right matrix for the multiplication
 */
void multiply_su3matrices(hmc_su3matrix *out, hmc_su3matrix *p, hmc_su3matrix *q);
/**
 * Multiply two staple matrices
 *
 * @param[out] out The matrix into which to store the multiplication result
 * @param[in] p Left matrix for the multiplication
 * @param[in] q Right matrix for the multiplication
 */
void multiply_staplematrix(hmc_staplematrix *out, hmc_su3matrix *p, hmc_staplematrix *q);

/**
 * Accumulate SU3 matrices into a staple-matrix.
 *
 * @param[in,out] p The staple to accumulate into. Needs to be null element on first invocation.
 * @param[in] q The SU3 matrix to add to the accumulation
 */
void accumulate_su3matrices_add(hmc_staplematrix *p, hmc_su3matrix *q);

#endif
