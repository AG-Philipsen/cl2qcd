/** @file
 * Operations on 3x3 matrices
 */
#ifndef _OPERATIONS_MATRIXH_
#define _OPERATIONS_MATRIXH_

#include "types.h"
#include "host_operations_complex.h"

/**
 * Multiplies two 3x3 matrices
 * @param[out] out The matrix into which to store the multiplication result
 * @param[in] p Left matrix for the multiplication
 * @param[in] q Right matrix for the multiplication
 */
void multiply_3x3matrix (hmc_3x3matrix *out, const hmc_3x3matrix *p, const hmc_3x3matrix *q);
/**
 * Adds two 3x3 matrices
 * @param[out] out The matrix into which to store the sum result
 * @param[in] p Left matrix for the summation
 * @param[in] q Right matrix for the summation
 */
void add_3x3matrix (hmc_3x3matrix *out, hmc_3x3matrix *p, hmc_3x3matrix *q);
/**
 * Subtracts two 3x3 matrices
 * @param[out] out The matrix into which to store the difference result
 * @param[in] p Left matrix for the subtraction
 * @param[in] q Right matrix for the subtraction
 */
void subtract_3x3matrix (hmc_3x3matrix *out, hmc_3x3matrix *p, hmc_3x3matrix *q);

/**
* Sets a 3x3 matrix to identity
* @param[out] mat The matrix into which to write the identity
*/
void set_to_3x3_identity(hmc_3x3matrix *mat);

/**
* Copy a 3x3 matrix into another (implements: *dest = *src, basically)
* @param[out] dest The matrix into which to copy
* @param[in] src The matrix to copy from
*/
void copy_3x3_matrix(hmc_3x3matrix *dest, hmc_3x3matrix *src);

/**
* Scale a 3x3 matrix by a real factor, in-place
* @param[out] mat The matrix to scale
* @param[in] factor A real scale factor
*/
void multiply_3x3matrix_by_real(hmc_3x3matrix *mat, hmc_float factor);

/**
* Copies a SU(3) matrix into a 3x3 matrix
* @param[out] out 3x3 matrix out
* @param[in]  in SU(3) matrix in
* 
*/
void su3matrix_to_3x3matrix (hmc_3x3matrix * out, hmc_su3matrix * in);

/**
 * Accumulates a SU3-Matrix into an 3x3-Matrix
 * @param[out] out The matrix into which to store the accumulation
 * @param[in]  q matrix for the accumulation
 */
void accumulate_su3matrix_3x3_add(hmc_3x3matrix *out, hmc_su3matrix *q);

/**
 * Computes the trace of an 3x3-Matrix
 * @param[out] out The result of tracing
 * @param[in]  q matrix for tracing
 */
void trace_3x3matrix (hmc_complex *out, hmc_3x3matrix *q);

/**
 * Computes the adjoint of an 3x3-Matrix
 * @param[out] out The adjoint matrix
 * @param[in]  q matrix to adjoin
 */
void adjoint_3x3matrix (hmc_3x3matrix * out, hmc_3x3matrix *q);



/**
* evaluates the sum of abs of difference element-by-element between two matrices
* test it
* @param[out] result the final ``absolute difference'' of the two matrices
* @param[in] mat2 matrix1 to compare
* @param[in] mat2 matrix2 to compare
*/
void absoluteDifference_3x3_matrix(hmc_float *result, hmc_3x3matrix *mat1, hmc_3x3matrix *mat2);

#endif

