/** @file
 * Handling of lattice geometry.
 *
 * @todo The conventions used here must be the same as in the opencl code. How can this be done automatically?
 */
#ifndef _GEOMETRYH_
#define _GEOMETRYH_

#include "globaldefs.h"
#include <cmath>
#include "meta/inputparameters.hpp"

/** Identify each spacetime direction */
#define TDIR 0
#define XDIR 1
#define YDIR 2
#define ZDIR 3

/**
 * Calculate the spatial index of the given cartesian coordinates.
 *
 * @param coord Pointer to NDIM integers representing cartesian coordinates.
 *              time is expected in index 0 of this array.
 * @param Spatial index
 */
int get_nspace(int* coord, const meta::Inputparameters& params);

/**
 * Get the non-even-odd-preconditioned index of a site based on the spatial and temporal
 * index.
 *
 * @param spacepos Spatial index
 * @param t Temporal index
 * @return Global index
 */
int get_global_pos(int spacepos, int t, const meta::Inputparameters& params);

/**
 * Get the non-even-odd-preconditioned index link based on the spatial, temporal
 * and Dirac index.
 *
 * @param spacepos Spatial index
 * @param t Temporal index
 * @return Global index
 */
int get_global_link_pos(int mu, int spacepos, int t, const meta::Inputparameters& params);

/**
 * This returns the index of a single su3 matrix entry in the ildg format
 * which is [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
 */
size_t get_su3_idx_ildg_format(size_t n, size_t m, size_t x, size_t y, size_t z, size_t t, size_t mu, const meta::Inputparameters& parameters);

/**
 * This returns the link index in the ildg format
 * which is [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
 */
size_t get_link_idx_ildg_format(size_t x, size_t y, size_t z, size_t t, size_t mu, const meta::Inputparameters& parameters);

int get_source_pos_spatial(const meta::Inputparameters& params);
#endif /* _GEOMETRYH_ */
