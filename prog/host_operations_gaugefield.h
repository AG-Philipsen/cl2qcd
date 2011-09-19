/** @file
 * Operations on the gauge field.
 *
 * @todo A lot of these operations use pointers and in-place operation
 *       without need. These might prevent compiler optimizations.
 */
#ifndef _OPERATIONS_GAUGEFIELDH_
#define _OPERATIONS_GAUGEFIELDH_
#include <iostream>
#include "globaldefs.h"
#include "types.h"
#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_operations_su3matrix.h"
// #include "host_operations_gaugemomentum.h"
#include "inputparameters.h"
#include <cmath>

/**
 * Create the gaugefield from an array of floats as used by by ILDG.
 *
 * @param[out] gaugefield Pointer to the new storage location.
 * @param[in] gaugefield_tmp Field in IDLG format
 * @param[in] check Size of the ILDG field.
 * @todo Replace hmc_gaugefield type by s_gaugefield type (LZ)
 */
void copy_gaugefield_from_ildg_format(hmc_complex * gaugefield, hmc_float * gaugefield_tmp, int check, const inputparameters * const params);
/**
 * Create the IDLG representation of the given gaugefield.
 *
 * @param[out] dest The location to store the ILDG representation to
 * @param[in] source The gaugefield in the internal representation
 * @todo Replace hmc_gaugefield type by s_gaugefield type (LZ)
 */
void copy_gaugefield_to_ildg_format(hmc_float * dest, hmc_complex * source, const inputparameters * const params);
/**
 * Create a representation of the gaugefield usable by the OpenCL kernels.
 *
 * @param[in] host_gaugefield The gaugefield in the internal representation
 * @param[out] gaugefield The location to store the OpenCL kernel compatible representation to
 */
void copy_to_ocl_format(ocl_s_gaugefield* host_gaugefield, Matrixsu3* gaugefield, const inputparameters * const params);

/**
 * Transform the gaugefield representation used by the OpenCL kernels into the normal one.
 *
 * @param[in] gaugefield The gaugefield in the representation used by the OpenCL kernels
 * @param[out] host_gaugefield The location to store the gaugefield to
 */
void copy_from_ocl_format(Matrixsu3* gaugefield, ocl_s_gaugefield* host_gaugefield, const inputparameters * const params);


/* *****************************************************************************************************
LZ: Note that the following section provides functions that work on the old hmc_gaugefield format.
    For the new format (s_gaugefield) they can still be used in combination with the corresponding copy to/from
    functions that exist in the Gaugefield class.
****************************************************************************************************/


/**
 * Retrieve an SU3 matrix form the gaugefield
 *
 * @param[in] spacepos Spatial index of the matrix to retrieve
 * @param[in] timepos Temporal index of the matrix to retrieve
 * @param[in] mu Direction of the matrix to retrieve
 * @return The index to be applied on the gaugefield in [NC][NC][NDIM][VOLSPACE][NTIME]
 *         or [NC*(NC-1)][NDIM][VOLSPACE][NTIME] format, depending on whether REC12 is enabled.
 */
#ifdef _RECONSTRUCT_TWELVE_
size_t get_hmc_gaugefield_index(size_t ncolor, size_t spacepos, size_t timepos, size_t mu, const inputparameters * const parameters);
#else
size_t get_hmc_gaugefield_index(size_t m, size_t n, size_t spacepos, size_t timepos, size_t mu, const inputparameters * const parameters);
#endif

/**
 * Retrieve an SU3 matrix form the gaugefield
 *
 * @param[out] out SU3 matrix to write the result to
 * @param[in] gaugefield from which to retrieve the SU3 matrix
 * @param[in] spacepos Spatial index of the matrix to retrieve
 * @param[in] timepos Temporal index of the matrix to retrieve
 * @param[in] mu Direction of the matrix to retrieve
 * @TODO CP: perhaps a similar function would be more efficient in same cases that just returns a pointer to the specific matrix in the "big" gaugefield-array
 */
void get_su3matrix(hmc_su3matrix* out, hmc_complex * in, int spacepos, int timepos, int mu, const inputparameters * const parameters); //cl
/**
 * Stores an SU3 matrix to the gaugefield
 *
 * @param[in,out] field from which to retrieve the SU3 matrix
 * @param[in] int SU3 matrix to store
 * @param[in] spacepos Spatial index of the matrix to store
 * @param[in] timepos Temporal index of the matrix to store
 * @param[in] mu Direction of the matrix to store
 */
void put_su3matrix(hmc_complex * field, hmc_su3matrix * in, int spacepos, int timepos, int mu, const inputparameters * const parameters); //cl

/**
 * Adjoin all SU3 matrices in a gaugefield.
 *
 * @param[in] in Gaugefield to read the SU3 matrices from
 * @param[out] out Gaugefield to store the matrices to
 */
void adjoin_su3(hmc_complex * in, hmc_complex * out, const inputparameters * const parameters);
/**
 * Sum up the traces of all SU3 matrices in on direction of the gaugefield.
 *
 * @param field The gaugefield to sum over
 * @param mu The direction to sum over
 */
hmc_complex global_trace_su3(hmc_complex * field, int mu, const inputparameters * const parameters);



/**
 * Calculate the part of the (temporal) polyakov loop local to the given spatial index.
 *
 * @param[in] field The gaugefield to use
 * @param[out] The local part of the polyakov
 * @param[in] n The spatial index to use
 */
void local_polyakov(hmc_complex * field, hmc_su3matrix * prod, int n, const inputparameters * const parameters);
/**
 * Calculate the part of the plaquette local to the given coordinates.
 *
 * @param[in] field The gaugefield to use
 * @param[out] The local part of the plaquette
 * @param[in] n The spatial index to use
 * @param[in] t The temporal index to use
 * @param[in] mu The first direction to use
 * @param[in] nu The second direction to use
 */
void local_plaquette(hmc_complex * field, hmc_su3matrix * prod, int n, int t, int mu, int nu, const inputparameters * const params);

// copy-functions within cpu memory, gaugefield-related layers
/**
 * Create a copy of the gaugefield.
 *
 * @param[in] source The gaugefield to copy
 * @param[ou] dest The storage location for the copy
 */
void copy_gaugefield(hmc_complex * source, hmc_complex * dest, const inputparameters * const params);

/**
 * Calculates the Q-plaquette (Clover discretization) at the given spacetime-point and directions
 *  --   --
 * |  | |  |
 *  --   --
 *  --   --
 * |  | |  |
 *  --   --
 * @param[in] field The gaugefield to use
 * @param[out] out The local part of the Qplaquette
 * @param[in] n The spatial index to use
 * @param[in] t The temporal index to use
 * @param[in] mu The first direction to use
 * @param[in] nu The second direction to use
 */
void local_Q_plaquette(hmc_3x3matrix * out, hmc_complex * field, int n, int t, int mu, int nu, const inputparameters * const params);

#endif
