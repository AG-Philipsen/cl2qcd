/** @file
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

#include "matrixSu3_utilities.hpp"
#include <stdexcept>
#include <string>

using namespace Matrixsu3_utilities;

Matrixsu3 createMatrixSu3BasedOnFillType(FillType fillTypeIn)
{
	if (fillTypeIn == ONE)
	{
		return { {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0} };
	}
	else if (fillTypeIn == ZERO)
	{
		return { {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0} };
	}
	else
	{
		throw std::logic_error("Do not know fill type for matrixSu3!");
	}
}

void Matrixsu3_utilities::fillMatrixSu3Array_constantMatrix(Matrixsu3 * in, size_t numberOfElements, FillType fillType)
{
	for(int iteration = 0; iteration < (int) numberOfElements; iteration ++)
	{
		in[iteration] = createMatrixSu3BasedOnFillType(fillType);
	}
}

hmc_complex Matrixsu3_utilities::sumUpAllMatrixElements(Matrixsu3 * in, size_t numberOfElements)
{
	hmc_complex sum = {0., 0.};
	for(int iteration = 0; iteration < (int) numberOfElements; iteration ++)
	{
		hmc_complex local_sum = {0.,0.};
		local_sum = complexadd(local_sum, in[iteration].e00);
		local_sum = complexadd(local_sum, in[iteration].e01);
		local_sum = complexadd(local_sum, in[iteration].e02);
		local_sum = complexadd(local_sum, in[iteration].e10);
		local_sum = complexadd(local_sum, in[iteration].e11);
		local_sum = complexadd(local_sum, in[iteration].e12);
		local_sum = complexadd(local_sum, in[iteration].e20);
		local_sum = complexadd(local_sum, in[iteration].e21);
		local_sum = complexadd(local_sum, in[iteration].e22);
		
		sum = complexadd(sum, local_sum);
	}
	
	return sum;
}
