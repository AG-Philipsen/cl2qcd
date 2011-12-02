/** @file
 * Handling of lattice geometry.
 *
 * @todo The conventions used here must be the same as in the opencl code. How can this be done automatically?
 */
#ifndef _GEOMETRYH_
#define _GEOMETRYH_

#include "globaldefs.h"
#include <cmath>
#include "inputparameters.h"

//coord[0] = t
//coord[1] = x
//coord[2] = y
//coord[3] = z

/**
 * Calculate the spatial index of the given cartesian coordinates.
 *
 * @param coord Pointer to NDIM integers representing cartesian coordinates.
 *              time is expected in index 0 of this array.
 * @param Spatial index
 */
int get_nspace(int* coord, const inputparameters * const params);

/**
 * Get the non-even-odd-preconditioned index of a site based on the spatial and temporal
 * index.
 *
 * @param spacepos Spatial index
 * @param t Temporal index
 * @return Global index
 */
int get_global_pos(int spacepos, int t, const inputparameters * const params);

/**
 * Get the non-even-odd-preconditioned index link based on the spatial, temporal
 * and Dirac index.
 *
 * @param spacepos Spatial index
 * @param t Temporal index
 * @return Global index
 */
int get_global_link_pos(int mu, int spacepos, int t, const inputparameters * const params);

/**
 * Retrieve an SU3 matrix form the gaugefield
 *
 * @param[in] spacepos Spatial index of the matrix to retrieve
 * @param[in] timepos Temporal index of the matrix to retrieve
 * @param[in] mu Direction of the matrix to retrieve
 * @return The index to be applied on the gaugefield in [NC][NC][NDIM][VOLSPACE][NTIME]
 *         or [NC*(NC-1)][NDIM][VOLSPACE][NTIME] format, depending on whether REC12 is enabled.
 */
size_t get_hmc_gaugefield_index(size_t m, size_t n, size_t spacepos, size_t timepos, size_t mu, const inputparameters * const parameters);

#endif /* _GEOMETRYH_ */
