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
#include "host_random.h"
#include "meta/inputparameters.hpp"
#include <cmath>

/**
 * Calculate the part of the (temporal) polyakov loop local to the given spatial index.
 *
 * @param[in] field The gaugefield to use
 * @return The local part of the polyakov
 * @param[in] n The spatial index to use
 */
Matrixsu3 local_polyakov(Matrixsu3 * field, int n, const meta::Inputparameters parameters);

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
Matrixsu3 local_plaquette(Matrixsu3 * field, int coord_in[NDIM], int mu, int nu, const meta::Inputparameters parameters);

/**
 * Stores an SU3 matrix to the gaugefield
 *
 * @param[in] field from which to retrieve the SU3 matrix
 * @param[in] SU3 matrix to store
 * @param[in] spacepos Spatial index of the matrix to store
 * @param[in] timepos Temporal index of the matrix to store
 * @param[in] mu Direction of the matrix to store
 */
void put_matrixsu3(Matrixsu3 * field, Matrixsu3 in, int spacepos, int timepos, int mu, const meta::Inputparameters parameters);

/**
 * Returns an SU3 matrix form the gaugefield
 *
 * @param[in] gaugefield from which to retrieve the SU3 matrix
 * @param[in] spacepos Spatial index of the matrix to retrieve
 * @param[in] timepos Temporal index of the matrix to retrieve
 * @param[in] mu Direction of the matrix to retrieve
 */
Matrixsu3 get_matrixsu3(Matrixsu3 * in, int spacepos, int timepos, int mu, const meta::Inputparameters parameters);

Matrixsu3 unit_matrixsu3();

Matrixsu3 random_matrixsu3();

#endif
