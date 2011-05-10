/** @file
 * Operations on SU3 matrices and the gauge field.
 *
 * @todo A lot of these operations use pointers and in-place operation
 *       without need. These might prevent compiler optimizations.
 */
#ifndef _OPERATIONS_GAUGEFIELDH_
#define _OPERATIONS_GAUGEFIELDH_
#include <iostream>
#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"
#include "host_geometry.h"
#include "host_operations_complex.h"
#include <cmath>

//gaugefield and su3 operations, global and local

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
 * Set all matrices in the gaugefield to unit matrices, representing a cold state.
 *
 * @param[out] field Pointer to a field of SU3 matrices
 * @return Error code as defined in hmcerrs.h
 */
hmc_error set_gaugefield_cold(hmc_gaugefield * field);
/**
 * Set all matrices in the gaugefield gausian distributed random matrices, representing a hot state.
 *
 * @param[out] field Pointer to a field of SU3 matrices
 * @return Error code as defined in hmcerrs.h
 */
hmc_error set_gaugefield_hot(hmc_gaugefield * field);

/**
 * Create the gaugefield from an array of floats as used by by ILDG.
 *
 * @todo Why doesn't this use ildg_gaugefield?
 *
 * @param[out] gaugefield Pointer to the new storage location.
 * @param[in] gaugefield_tmp Field in IDLG format
 * @param[in] check Size of the ILDG field.
 * @return Error code as defined in hmcerrs.h
 */
hmc_error copy_gaugefield_from_ildg_format(hmc_gaugefield * gaugefield, hmc_float * gaugefield_tmp, int check);
/**
 * Create the IDLG representation of the given gaugefield.
 *
 * @param[out] dest The location to store the ILDG representation to
 * @param[in] source The gaugefield in the internal representation
 * @return Error code as defined in hmcerrs.h
 */
hmc_error copy_gaugefield_to_ildg_format(ildg_gaugefield * dest, hmc_gaugefield * source);
/**
 * Create a representation of the gaugefield usable by the OpenCL kernels.
 *
 * @param[in] host_gaugefield The gaugefield in the internal representation
 * @param[out] gaugefield The location to store the OpenCL kernel compatible representation to
 * @return Error code as defined in hmcerrs.h
 */
hmc_error copy_to_ocl_format(hmc_ocl_gaugefield* host_gaugefield,hmc_gaugefield* gaugefield);
/**
 * Transform the gaugefield representation used by the OpenCL kernels into the normal one.
 *
 * @param[in] gaugefield The gaugefield in the representation used by the OpenCL kernels
 * @param[out] host_gaugefield The location to store the gaugefield to
 * @return Error code as defined in hmcerrs.h
 */
hmc_error copy_from_ocl_format(hmc_gaugefield* gaugefield,hmc_ocl_gaugefield* host_gaugefield);

/**
 * Retrieve an SU3 matrix form the gaugefield
 *
 * @param[out] out SU3 matrix to write the result to
 * @param[in] gaugefield from which to retrieve the SU3 matrix
 * @param[in] spacepos Spatial index of the matrix to retrieve
 * @param[in] timepos Temporal index of the matrix to retrieve
 * @param[in] mu Direction of the matrix to retrieve
 * @return Error code as defined in hmcerrs.h
 */
hmc_error get_su3matrix(hmc_su3matrix* out, hmc_gaugefield * in, int spacepos, int timepos, int mu); //cl
/**
 * Stores an SU3 matrix to the gaugefield
 *
 * @param[in,out] field from which to retrieve the SU3 matrix
 * @param[in] int SU3 matrix to store
 * @param[in] spacepos Spatial index of the matrix to store
 * @param[in] timepos Temporal index of the matrix to store
 * @param[in] mu Direction of the matrix to store
 * @return Error code as defined in hmcerrs.h
 */
hmc_error put_su3matrix(hmc_gaugefield * field, hmc_su3matrix * in, int spacepos, int timepos, int mu); //cl

/**
 * Adjoin all SU3 matrices in a gaugefield.
 *
 * @param[in] in Gaugefield to read the SU3 matrices from
 * @param[out] out Gaugefield to store the matrices to
 * @return Error code as defined in hmcerrs.h
 */
hmc_error adjoin_su3(hmc_gaugefield * in, hmc_gaugefield * out);
/**
 * Sum up the traces of all SU3 matrices in on direction of the gaugefield.
 *
 * @param field The gaugefield to sum over
 * @param mu The direction to sum over
 * @return Error code as defined in hmcerrs.h
 */
hmc_complex global_trace_su3(hmc_gaugefield * field, int mu);
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
 * Calculate the part of the (temporal) polyakov loop local to the given spatial index.
 *
 * @param[in] field The gaugefield to use
 * @param[out] The local part of the polyakov
 * @param[in] n The spatial index to use
 */
void local_polyakov(hmc_gaugefield * field, hmc_su3matrix * prod, int n);
/**
 * Calculate the part of the plaquette local to the given coordinates.
 *
 * @param[in] field The gaugefield to use
 * @param[out] The local part of the plaquette
 * @param[in] n The spatial index to use
 * @param[in] n The temporal index to use
 * @param[in] mu The first direction to use
 * @param[in] nu The second direction to use
 */
void local_plaquette(hmc_gaugefield * field, hmc_su3matrix * prod, int n, int t, int mu, int nu );

// copy-functions within cpu memory, gaugefield-related layers
/**
 * Create a copy of the gaugefield.
 *
 * @param[in] source The gaugefield to copy
 * @param[ou] dest The storage location for the copy
 * @return Error code as defined in hmcerrs.h
 */
hmc_error copy_gaugefield(hmc_gaugefield * source, hmc_gaugefield * dest);

//gauge-momenta operations
//TODO CP: these should go into a seperate file like host_operations_gaugemomenta.cpp

/**
 * Create a copy of the gauge momenta.
 *
 * @param[in] source The gauge momenta to copy
 * @param[out] dest The storage location for the copy
 * @return Error code as defined in hmcerrs.h
 */
hmc_error copy_gaugemomenta(hmc_gauge_momentum * source, hmc_gauge_momentum * dest);
/**
 * Calculate the squarenorm of the gauge momenta.
 *
 * @param[in] in The gauge momenta to use.
 * @param[out] result The square norm
 * @return Error code as defined in hmcerrs.h
 */
hmc_error gaugemomenta_squarenorm(hmc_gauge_momentum * in, hmc_float * result);
/**
 * Set gaugemomenta to zero.
 *
 * @param[out] in The gauge momenta to set to zero.
 * @return Error code as defined in hmcerrs.h
 */
hmc_error set_zero_gaugemomenta(hmc_gauge_momentum * in);

#endif
