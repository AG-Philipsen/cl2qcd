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
#include "hmcerrs.h"
#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_operations_su3matrix.h"
// #include "host_operations_gaugemomentum.h"
#include <cmath>

/**
 * Create the gaugefield from an array of floats as used by by ILDG.
 *
 * @todo Why doesn't this use ildg_gaugefield?
 *
 * @param[out] gaugefield Pointer to the new storage location.
 * @param[in] gaugefield_tmp Field in IDLG format
 * @param[in] check Size of the ILDG field.
 * @return Error code as defined in hmcerrs.h
 * @todo Replace hmc_gaugefield type by s_gaugefield type (LZ)
 */
hmc_error copy_gaugefield_from_ildg_format(hmc_gaugefield * gaugefield, hmc_float * gaugefield_tmp, int check);
/**
 * Create the IDLG representation of the given gaugefield.
 *
 * @param[out] dest The location to store the ILDG representation to
 * @param[in] source The gaugefield in the internal representation
 * @return Error code as defined in hmcerrs.h
 * @todo Replace hmc_gaugefield type by s_gaugefield type (LZ) 
 */
hmc_error copy_gaugefield_to_ildg_format(ildg_gaugefield * dest, hmc_gaugefield * source);
/**
 * Create a representation of the gaugefield usable by the OpenCL kernels.
 *
 * @param[in] host_gaugefield The gaugefield in the internal representation
 * @param[out] gaugefield The location to store the OpenCL kernel compatible representation to
 * @return Error code as defined in hmcerrs.h
 */
hmc_error copy_to_ocl_format(ocl_s_gaugefield* host_gaugefield, s_gaugefield* gaugefield);

/**
 * Transform the gaugefield representation used by the OpenCL kernels into the normal one.
 *
 * @param[in] gaugefield The gaugefield in the representation used by the OpenCL kernels
 * @param[out] host_gaugefield The location to store the gaugefield to
 * @return Error code as defined in hmcerrs.h
 */
hmc_error copy_from_ocl_format(s_gaugefield* gaugefield, ocl_s_gaugefield* host_gaugefield);


/* *****************************************************************************************************
LZ: Note that the following section provides functions that work on the old hmc_gaugefield format.
    For the new format (s_gaugefield) they can still be used in combination with the corresponding copy to/from
    functions that exist in the Gaugefield class.
****************************************************************************************************/


/**
 * Retrieve an SU3 matrix form the gaugefield
 *
 * @param[out] out SU3 matrix to write the result to
 * @param[in] gaugefield from which to retrieve the SU3 matrix
 * @param[in] spacepos Spatial index of the matrix to retrieve
 * @param[in] timepos Temporal index of the matrix to retrieve
 * @param[in] mu Direction of the matrix to retrieve
 * @return Error code as defined in hmcerrs.h
 * @TODO CP: perhaps a similar function would be more efficient in same cases that just returns a pointer to the specific matrix in the "big" gaugefield-array
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
 * @param[in] t The temporal index to use
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
void local_Q_plaquette(hmc_3x3matrix * out, hmc_gaugefield * field, int n, int t, int mu, int nu );

#endif
