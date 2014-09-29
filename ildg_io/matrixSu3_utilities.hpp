/** @file
 * Utilities for su3 matrices.
 *
 * Copyright 2014, Christopher Pinke
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

#ifndef _MATRIXSU3UTILITIES_HPP_
#define _MATRIXSU3UTILITIES_HPP_

#include "../common_header_files/types.h"
#include "../common_header_files/operations_complex.h"

#include <vector>

namespace Matrixsu3_utilities
{
	enum FillType {ZERO, ONE, DIAGONAL, RANDOM, FILLED};

	void fillMatrixSu3Array_constantMatrix(std::vector<Matrixsu3> & in, FillType fillType);
	void fillMatrixSu3Array_randomMatrix(std::vector<Matrixsu3> & in);
	hmc_complex sumUpAllMatrixElements(const std::vector<Matrixsu3> & in);
	hmc_complex sumUpDiagonalMatrixElements(const std::vector<Matrixsu3> & in);
	hmc_complex sumUpOffDiagonalMatrixElements(const std::vector<Matrixsu3> & in);
	Matrixsu3 getUnitMatrix();
	Matrixsu3 getZeroMatrix();
	Matrixsu3 getFilledMatrix();
}

#endif