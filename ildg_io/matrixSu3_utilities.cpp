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

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>

using namespace Matrixsu3_utilities;

Matrixsu3 Matrixsu3_utilities::getUnitMatrix()
{
	return { {1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0} };
}

Matrixsu3 Matrixsu3_utilities::getZeroMatrix()
{
	return { {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0} };
}

Matrixsu3 getRandomMatrix(boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution < > > generateRandomNumbers)
{
	Matrixsu3 tmp = getZeroMatrix();
	
	tmp.e00.re = generateRandomNumbers();
	tmp.e01.re = generateRandomNumbers();
	tmp.e02.re = generateRandomNumbers();
	tmp.e10.re = generateRandomNumbers();
	tmp.e11.re = generateRandomNumbers();
	tmp.e12.re = generateRandomNumbers();
	tmp.e20.re = generateRandomNumbers();
	tmp.e21.re = generateRandomNumbers();
	tmp.e22.re = generateRandomNumbers();
	
	tmp.e00.im = generateRandomNumbers();
	tmp.e01.im = generateRandomNumbers();
	tmp.e02.im = generateRandomNumbers();
	tmp.e10.im = generateRandomNumbers();
	tmp.e11.im = generateRandomNumbers();
	tmp.e12.im = generateRandomNumbers();
	tmp.e20.im = generateRandomNumbers();
	tmp.e21.im = generateRandomNumbers();
	tmp.e22.im = generateRandomNumbers();
	
	return tmp;
}

Matrixsu3 createMatrixSu3BasedOnFillType(FillType fillTypeIn)
{
	if (fillTypeIn == ONE)
	{
		return { {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0} };
	}
	else if (fillTypeIn == ZERO)
	{
		return getZeroMatrix();
	}
	else if (fillTypeIn == DIAGONAL)
	{
		return getUnitMatrix();
	}
	else
	{
		throw std::logic_error("Do not know fill type for matrixSu3!");
	}
}

void Matrixsu3_utilities::fillMatrixSu3Array_constantMatrix(std::vector<Matrixsu3> & in, FillType fillType)
{
	for(int iteration = 0; iteration < (int) in.size(); iteration ++)
	{
		in[iteration] = createMatrixSu3BasedOnFillType(fillType);
	}
}

void Matrixsu3_utilities::fillMatrixSu3Array_randomMatrix(std::vector<Matrixsu3> & in)
{
	boost::mt19937 randomNumbergenerator( time( 0 ) );
	boost::random::uniform_real_distribution< double > uniformDistribution( 0.0, 1.0 );
	boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution < > > 
	generateRandomNumbers( randomNumbergenerator, uniformDistribution );
	
	for(int iteration = 0; iteration < (int) in.size(); iteration ++)
	{
		in[iteration] = getRandomMatrix(generateRandomNumbers);
	}
}

hmc_complex Matrixsu3_utilities::sumUpAllMatrixElements(const std::vector<Matrixsu3> & in)
{
	hmc_complex sum = {0., 0.};
	for(int iteration = 0; iteration < (int) in.size(); iteration ++)
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

hmc_complex Matrixsu3_utilities::sumUpDiagonalMatrixElements(const std::vector<Matrixsu3> & in)
{
	hmc_complex sum = {0., 0.};
	for(int iteration = 0; iteration < (int) in.size(); iteration ++)
	{
		hmc_complex local_sum = {0.,0.};
		local_sum = complexadd(local_sum, in[iteration].e00);
		local_sum = complexadd(local_sum, in[iteration].e11);
		local_sum = complexadd(local_sum, in[iteration].e22);
		
		sum = complexadd(sum, local_sum);
	}
	
	return sum;
}

hmc_complex Matrixsu3_utilities::sumUpOffDiagonalMatrixElements(const std::vector<Matrixsu3> & in)
{
	hmc_complex sum = {0., 0.};
	for(int iteration = 0; iteration < (int) in.size(); iteration ++)
	{
		hmc_complex local_sum = {0.,0.};
		local_sum = complexadd(local_sum, in[iteration].e01);
		local_sum = complexadd(local_sum, in[iteration].e02);
		local_sum = complexadd(local_sum, in[iteration].e10);
		local_sum = complexadd(local_sum, in[iteration].e12);
		local_sum = complexadd(local_sum, in[iteration].e20);
		local_sum = complexadd(local_sum, in[iteration].e21);
		
		sum = complexadd(sum, local_sum);
	}
	
	return sum;
}

