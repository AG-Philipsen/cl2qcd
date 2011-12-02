/** @file
 * Handling of lattice geometry.
 *
 * @todo The conventions used here must be the same as in the opencl code.
 * 	Therefore, include the same file here!!!
 */
#ifndef _GEOMETRYH_
#define _GEOMETRYH_

#include "globaldefs.h"
#include <cstdlib>
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
 * Calculate the cartesian coordinates of the given spatial index.
 *
 * @param nspace Spatial index
 * @param dir The dimension for which to retrieve the cartisian coordinate
 * @param Cartesian coordinate of nspace in dimension dir.
 */
int get_spacecoord(int nspace, int dir, const inputparameters * const params);

/**
 * Calculate the next spatial index in a given direction.
 *
 * @param nspace The spatial index to start with
 * @param dir The dimension in which to move
 * @return Index of the next site in the given direction
 */
int get_neighbor(int nspace, int dir, const inputparameters * const params);

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

#endif /* _GEOMETRYH_ */
