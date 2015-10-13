/** @file
 * Handling of lattice geometry.
 *
 * @todo The conventions used here must be the same as in the opencl code. How can this be done automatically?
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _GEOMETRYH_
#define _GEOMETRYH_

#include "../common_header_files/globaldefs.h"
#include <cmath>
#include "../meta/inputparameters.hpp"
#include "../hardware/size_4.hpp"

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
int get_nspace(int* coord, const int nt, const int ns);

/**
 * Get the non-even-odd-preconditioned index of a site based on the spatial and temporal
 * index.
 *
 * @param spacepos Spatial index
 * @param t Temporal index
 * @return Global index
 */
int get_global_pos(int spacepos, int t, const meta::Inputparameters& params);
int get_global_pos(int spacepos, int t, const int nt, const int ns);

/**
 * Get the non-even-odd-preconditioned index based on cartesian coordinates.
 *
 * @param cart Cartisian coordinates
 * @return Global index
 */
int get_global_pos(size_4 cart, const meta::Inputparameters& params);
int get_global_pos(size_4 cart, const int nt, const int n);

/**
 * Get the non-even-odd-preconditioned index link based on the spatial, temporal
 * and Dirac index.
 *
 * @param spacepos Spatial index
 * @param t Temporal index
 * @return Global index
 */
int get_global_link_pos(int mu, int spacepos, int t, const meta::Inputparameters& params);
int get_global_link_pos(int mu, int spacepos, int t, const int nt, const int ns);

/**
 * Get the non-even-odd-preconditioned link index cartesian coordinates
 * and Dirac index.
 *
 * @param mu Dirac index
 * @param cart Cartisian coordinates
 * @return Global index
 */
int get_global_link_pos(int mu, size_4 cart, const meta::Inputparameters& params);
int get_global_link_pos(int mu, size_4 cart, const int nt, const int ns);

/**
 * This returns the index of a single su3 matrix entry in the ildg format
 * which is [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
 */
size_t get_su3_idx_ildg_format(size_t n, size_t m, size_t x, size_t y, size_t z, size_t t, size_t mu, const int nt, const int ns);


/**
 * This returns the link index in the ildg format
 * which is [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
 */
size_t get_link_idx_ildg_format(size_t x, size_t y, size_t z, size_t t, size_t mu, const int nt, const int ns);

int get_source_pos_spatial(const meta::Inputparameters& params);
#endif /* _GEOMETRYH_ */
