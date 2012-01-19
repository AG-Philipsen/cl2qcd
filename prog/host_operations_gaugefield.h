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
#include "inputparameters.h"
#include <cmath>

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
 * Calculate the part of the (temporal) polyakov loop local to the given spatial index.
 *
 * @param[in] field The gaugefield to use
 * @return The local part of the polyakov
 * @param[in] n The spatial index to use
 */
Matrixsu3 local_polyakov(Matrixsu3 * field, int n, const inputparameters * const parameters);
/**
 * Calculate the part of the plaquette local to the given coordinates.
 *
 * @param[in] field The gaugefield to use
 * @return The local part of the plaquette
 * @param[in] n The spatial index to use
 * @param[in] t The temporal index to use
 * @param[in] mu The first direction to use
 * @param[in] nu The second direction to use
 */
Matrixsu3 local_plaquette(Matrixsu3 * field, int coord_in[NDIM], int mu, int nu, const inputparameters * const parameters);

void put_matrixsu3(Matrixsu3 * field, Matrixsu3 in, int spacepos, int timepos, int mu, const inputparameters * const parameters);

Matrixsu3 convert_hmc_matrixsu3_to_Matrixsu3(hmc_su3matrix in);

Matrixsu3 unit_matrixsu3();

#endif
