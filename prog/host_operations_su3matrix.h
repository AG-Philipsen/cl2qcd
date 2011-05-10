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
#include "hmcerrs.h"
#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_operations_matrix.h"
#include <cmath>

/**
 * Replace the given matrix by it's adjoint.
 *
 * @param[in,out] mat The matrix to replace by it's adjoint
 * @return Error code as defined in hmcerrs.h
 */
hmc_error adjoin_su3matrix(hmc_su3matrix * mat);
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
 * @return Error code as defined in hmcerrs.h
 */
hmc_error copy_su3matrix(hmc_su3matrix *out, hmc_su3matrix *in); 
/**
 * Copy the contents of an staple matrix to another one.
 *
 * @param[out] The staple matrix to copy to
 * @param[in] The staple matrix to copy from
 * @return Error code as defined in hmcerrs.h
 */
hmc_error copy_staplematrix(hmc_staplematrix *out, hmc_staplematrix *in);
/**
 * Replace the given matrix by a unit matrix.
 *
 * @param[out] mat The matrix to replace by a unit one
 * @return Error code as defined in hmcerrs.h
 */
hmc_error unit_su3matrix(hmc_su3matrix * u); 
/**
 * Replace the given matrix by a random matrix.
 *
 * @param[out] mat The matrix to replace by a random one
 * @return Error code as defined in hmcerrs.h
 */
hmc_error random_su3matrix(hmc_su3matrix * u);
/**
 * Replace the given matrix by a zero matrix.
 *
 * @param[out] mat The matrix to replace by a zero one
 * @return Error code as defined in hmcerrs.h
 */
hmc_error zero_su3matrix(hmc_su3matrix * u); 
/**
 * Replace the given matrix by a zero matrix.
 *
 * @param[out] mat The matrix to replace by a zero one
 * @return Error code as defined in hmcerrs.h
 */
hmc_error zero_staplematrix(hmc_staplematrix * u);
/**
 * Multiply two SU3 matrices in place
 *
 * @param[in,out] acc The matrix to be replaced by the multiplication result
 * @param[in] multiplicator The matrix to multiply witht the acc matrix
 * @return Error code as defined in hmcerrs.h
 */
hmc_error accumulate_su3matrix_prod(hmc_su3matrix *acc, hmc_su3matrix *multiplicator);
/**
 * Multiply two SU3 matrices
 *
 * @param[out] out The matrix into which to store the multiplication result
 * @param[in] p Left matrix for the multiplication
 * @param[in] q Right matrix for the multiplication
 * @return Error code as defined in hmcerrs.h
 */
hmc_error multiply_su3matrices(hmc_su3matrix *out, hmc_su3matrix *p, hmc_su3matrix *q);
/**
 * Multiply two staple matrices
 *
 * @param[out] out The matrix into which to store the multiplication result
 * @param[in] p Left matrix for the multiplication
 * @param[in] q Right matrix for the multiplication
 * @return Error code as defined in hmcerrs.h
 */
hmc_error multiply_staplematrix(hmc_staplematrix *out, hmc_su3matrix *p, hmc_staplematrix *q);
#ifdef _RECONSTRUCT_TWELVE_
/**
 * Reconstruct the third row of a compressed SU3 matrix
 *
 * @param in The compressed SU3 matrix
 * @param ncomp The component of the row to reconstruct
 * @return The value of the chosen component
 */
hmc_complex reconstruct_su3(hmc_su3matrix *in, int ncomp);
#endif

/**
 * Accumulate SU3 matrices into a staple-matrix.
 *
 * @param[in,out] p The staple to accumulate into. Needs to be null element on first invocation.
 * @param[in] q The SU3 matrix to add to the accumulation
 * @return Error code as defined in hmcerrs.h
 */
hmc_error accumulate_su3matrices_add(hmc_staplematrix *p, hmc_su3matrix *q);
/**
 * Reduce an SU3 matrix to an SU2 (as in Cabbibo-Marinari)
 *
 * @param[out] dest SU(2)-Matrix to store to
 * @param[in] src The SU3 matrix to be reduced
 * @param[in] rand Which part of the SU(3)-Matrix to use
 * (this can be 1,2 or 3)
 */
void reduction (hmc_complex dest[su2_entries], hmc_staplematrix src, const int rand);
/**
 * Expand an SU(2)-Matrix to SU(3)
 *
 * @param[out] dest SU(3)-Matrix to store to
 * @param[in] rand This can be 1,2 or 3 and essentially gives the position of the 1 in the SU(3)-Matrix
 * @param[in] src The SU(2)-Matrix to extend
 */
void extend (hmc_su3matrix * dest, const int random, hmc_complex src[su2_entries]);
/**
 * Projects a 3x3-Matrix back to SU(3) using the Gram-Schmidt-Procedure.
 *
 * @param[in,out] U The SU3 matrix to project.
 * @return Error code as defined in hmcerrs.h
 */
hmc_error project_su3(hmc_su3matrix *U);
hmc_error project_su3_old(hmc_su3matrix *U);

/**
 * Apply Boundary Conditions to a SU(3)-Matrix. 
 * This corresponds to multiplying each component of the matrix by a (complex) factor of /f$\exp(i*\theta) /f$.
 * /f$\theta = 0/f$ are 'periodic' BC.
 * /f$\theta = \Pi/f$ are 'antiperiodic' BC.
 *
 * @param[in,out] in SU(3)-Matrix to be changed
 * @param[in] theta angle /f$\theta/f$
 * @return void
 * @todo the calculation involves sin- and cos-evaluations. Perhaps one should optimize this for the two special cases mentioned above.
 * 
 */
void gaugefield_apply_bc(hmc_su3matrix * in, hmc_float theta);

/**
 * Apply (complex) chemical Potential /f$\mu/f$ simultaneously to two SU(3)-Matrices. 
 * This corresponds to multiplying each component of the matrix by a (complex) factor of /f$\exp(\mu) /f$.
 *
 * @param[in,out] u SU(3)-Matrix to be changed
 * @param[in,out] udagger SU(3)-Matrix to be changed
 * @param[in] chem_pot_re real part of /f$\mu/f$
 * @param[in] chem_pot_im imaginary part of /f$\mu/f$
 * @return void
 * @remark In the OpenCL-part this function is explicitly splitted into real and imaginary chemical potential.
 * @remark Two Matrices are updated simultaneously since /f$\mu/f$ is usually applied in the /f$\notD/f$-operation involving two links.
 */
void gaugefield_apply_chem_pot(hmc_su3matrix * u, hmc_su3matrix * udagger, hmc_float chem_pot_re, hmc_float chem_pot_im);

/**
 * Calculates the SU(3)-Matrix /f$\exp(i\epsilon Q)/f$, where Q is a su(3)-algebra element (practically 8 real numbers). 
 * This can be done either by calculating the exponential series to order 2,3 or all orders or by applying the algorithm
 * provided by Morningstar-Peardon.
 * @param[in] in input su(3)-algebra element (8 real numbers)
 * @param[out] out output SU(3)-Matrix
 * @param[in] epsilon input parameter
 * @return Error code as defined in hmcerrs.h
 * @todo needs testing
 * @todo implement Morningstar-Peardon
 * @todo in the end this should be moved elsewhere since it is not specific to he hmc-algorithm
 */
hmc_error build_su3matrix_by_exponentiation(hmc_algebraelement in, hmc_su3matrix *out, hmc_float epsilon); 

#endif
