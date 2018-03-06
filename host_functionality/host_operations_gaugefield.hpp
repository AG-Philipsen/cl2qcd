/** @file
 * Operations on the gauge field.
 *
 * @todo A lot of these operations use pointers and in-place operation
 *       without need. These might prevent compiler optimizations.
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
#ifndef _OPERATIONS_GAUGEFIELDH_
#define _OPERATIONS_GAUGEFIELDH_
#include <iostream>
#include "../common_header_files/globaldefs.h"
#include "../common_header_files/types.h"
#include "../geometry/index.hpp"
#include "../meta/inputparameters.hpp"
#include <cmath>

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
Matrixsu3 local_plaquette(Matrixsu3 * field, int coord_in[NDIM], int mu, int nu, const meta::Inputparameters& parameters);

/**
 * Returns an SU3 matrix form the gaugefield
 *
 * @param[in] gaugefield from which to retrieve the SU3 matrix
 * @param[in] spacepos Spatial index of the matrix to retrieve
 * @param[in] timepos Temporal index of the matrix to retrieve
 * @param[in] mu Direction of the matrix to retrieve
 */
Matrixsu3 get_matrixsu3(Matrixsu3 * in, int spacepos, int timepos, int mu, const meta::Inputparameters& parameters);

Matrixsu3 unit_matrixsu3();
Matrixsu3 nonTrivialSu3Matrix();

#endif
