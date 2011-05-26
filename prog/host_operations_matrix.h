/** @file
 * Operations on 3x3 matrices
 */
#ifndef _OPERATIONS_MATRIXH_
#define _OPERATIONS_MATRIXH_

#include "types.h"
#include "hmcerrs.h"
#include "host_operations_complex.h"

/**
 * Multiplies two 3x3 matrices
 * @param[out] out The matrix into which to store the multiplication result
 * @param[in] p Left matrix for the multiplication
 * @param[in] q Right matrix for the multiplication
 * @return Error code as defined in hmcerrs.h
 */
hmc_error multiply_3x3matrix (hmc_3x3matrix *out, const hmc_3x3matrix *p, const hmc_3x3matrix *q);
/**
 * Adds two 3x3 matrices
 * @param[out] out The matrix into which to store the sum result
 * @param[in] p Left matrix for the summation
 * @param[in] q Right matrix for the summation
 * @return Error code as defined in hmcerrs.h
 */
hmc_error add_3x3matrix (hmc_3x3matrix *out, hmc_3x3matrix *p, hmc_3x3matrix *q);
/**
 * Subtracts two 3x3 matrices
 * @param[out] out The matrix into which to store the difference result
 * @param[in] p Left matrix for the subtraction
 * @param[in] q Right matrix for the subtraction
 * @return Error code as defined in hmcerrs.h
 */
hmc_error subtract_3x3matrix (hmc_3x3matrix *out, hmc_3x3matrix *p, hmc_3x3matrix *q);

/**
* Sets a 3x3 matrix to identity
* @param[out] mat The matrix into which to write the identity
* @return Error code as defined in hmcerrs.h
*/
hmc_error set_to_3x3_identity(hmc_3x3matrix *mat);

/**
* Copy a 3x3 matrix into another (implements: *dest = *src, basically)
* @param[out] dest The matrix into which to copy
* @param[in] src The matrix to copy from
* @return Error code as defined in hmcerrs.h
*/
hmc_error copy_3x3_matrix(hmc_3x3matrix *dest, hmc_3x3matrix *src);

/**
* Scale a 3x3 matrix by a real factor, in-place
* @param[out] mat The matrix to scale
* @param[in] factor A real scale factor
* @return Error code as defined in hmcerrs.h
*/
hmc_error multiply_3x3matrix_by_real(hmc_3x3matrix *mat, hmc_float factor);

/**
* Copies a SU(3) matrix into a 3x3 matrix
* @param[out] out 3x3 matrix out
* @param[in]  in SU(3) matrix in
* 
* @return Error code as defined in hmcerrs.h
*/
hmc_error su3matrix_to_3x3matrix (hmc_3x3matrix * out, hmc_su3matrix * in);

/**
 * Accumulates a SU3-Matrix into an 3x3-Matrix
 * @param[out] out The matrix into which to store the accumulation
 * @param[in]  q matrix for the accumulation
 * @return Error code as defined in hmcerrs.h
 */
hmc_error accumulate_su3matrix_3x3_add(hmc_3x3matrix *out, hmc_su3matrix *q);

/**
 * Computes the trace of an 3x3-Matrix
 * @param[out] out The result of tracing
 * @param[in]  q matrix for tracing
 * @return Error code as defined in hmcerrs.h
 */
hmc_error trace_3x3matrix (hmc_complex *out, hmc_3x3matrix *q);

/**
 * Computes the adjoint of an 3x3-Matrix
 * @param[out] out The adjoint matrix
 * @param[in]  q matrix to adjoin
 * @return Error code as defined in hmcerrs.h
 */
hmc_error adjoint_3x3matrix (hmc_3x3matrix * out, hmc_3x3matrix *q);



/**
* evaluates the sum of abs of difference element-by-element between two matrices
* test it
* @param[out] result the final ``absolute difference'' of the two matrices
* @param[in] mat2 matrix1 to compare
* @param[in] mat2 matrix2 to compare
* @return Error code as defined in hmcerrs.h
*/
hmc_error absoluteDifference_3x3_matrix(hmc_float *result, hmc_3x3matrix *mat1, hmc_3x3matrix *mat2);
/**
* Multiplies a T_i SU(3) generator (i.e. 1/2 lambda_i) by a generic 3x3 matrix
* WARNING: the generator index here runs from ONE to EIGHT !
* Hard-coded operations on single components!
* test it
* @param[out] out the result, T_i * M
* @param[in]  gen_index index i (1--8) of the generator T_i
* @param[in]  in input matrix M
* @return Error code as defined in hmcerrs.h
*/
hmc_error multiply_generator_3x3matrix (hmc_3x3matrix * out, int gen_index, hmc_3x3matrix *in);

/**
* Multiplies a generic 3x3 matrix by a T_i SU(3) generator (i.e. 1/2 lambda_i)
* WARNING: the generator index here runs from ONE to EIGHT !
* Hard-coded operations on single components!
* test it
* @param[out] out the result, M * T_i
* @param[in]  in input matrix M
* @param[in]  gen_index index i (1--8) of the generator T_i
* @return Error code as defined in hmcerrs.h
*/
hmc_error multiply_3x3matrix_generator (hmc_3x3matrix * out, hmc_3x3matrix *in, int gen_index);

/**
* Multiplies a T_i SU(3) generator (i.e. 1/2 lambda_i) by a SU(3) matrix
* WARNING: the generator index here runs from ONE to EIGHT !
* Relies on the "general 3x3 matrix" corresponding function
* The result is in any case a GENERIC 3x3 matrix
* test it
* @param[out] out the result, T_i * M
* @param[in]  gen_index index i (1--8) of the generator T_i
* @param[in]  in input matrix M
* @return Error code as defined in hmcerrs.h
*/
hmc_error multiply_generator_su3matrix (hmc_3x3matrix * out, int gen_index, hmc_su3matrix *in);

/**
* Multiplies a SU(3) matrix by a T_i SU(3) generator (i.e. 1/2 lambda_i)
* WARNING: the generator index here runs from ONE to EIGHT !
* Relies on the "general 3x3 matrix" corresponding function
* The result is in any case a GENERIC 3x3 matrix
* test it
* @param[out] out the result, T_i * M
* @param[in]  in input matrix M
* @param[in]  gen_index index i (1--8) of the generator T_i
* @return Error code as defined in hmcerrs.h
*/
hmc_error multiply_su3matrix_generator (hmc_3x3matrix * out, hmc_su3matrix *in, int gen_index);

/**
 * Called by build_su3matrix_by_exponentiation in case of "smart" approach.
 * Takes the 2*(8+1) real parameters beta_0, gamma_0, beta[8], gamma[8] and compiles all components
 * of the generic 3x3 complex matrix that is the linear combination of identity+generators:
 * /f[
 *     \exp(i\epsilon Q) = (\beta_0+i\gamma_0)\cdot 1 + \sum_\ell T_\ell (\beta_\ell + i \gamma_\ell)
 * /f]
 * the coefficients gamm and beta being calculated by the caller (as in Steo's notes) by contractions with 
 * the SU(3) symmetric structure constants with the 8 real parameters defining the algebra element
 * test it
 *
 * @param[in] beta_0 ... 
 * @param[in] gamma_0 ...
 * @param[in] beta ...
 * @param[in] gamma ...
 * @param[out] out output 3x3-Matrix
 * @return Error code as defined in hmcerrs.h
 * @todo needs testing
 */
hmc_error construct_3x3_combination(hmc_float beta_0, hmc_float gamma_0, hmc_float beta[], hmc_float gamma[], hmc_3x3matrix out); 

#endif

