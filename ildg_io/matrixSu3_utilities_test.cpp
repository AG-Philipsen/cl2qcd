/** @file
 * Testcases for the matrixsu3 utilities
 *
 * Copyright (c) 2014 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE matrixsu3_utilities
#include <boost/test/unit_test.hpp>

#include "matrixSu3_utilities.hpp"

using namespace Matrixsu3_utilities;

const int numMatrixEntries = 9;

class MatrixSu3Tester
{
public:
	MatrixSu3Tester(int ntimeIn, int nspaceIn) : 
		ntime(ntimeIn), nspace(nspaceIn)
	{
		vol4d = ntime * nspace * nspace * nspace;
		numberOfElements = vol4d * 4;
		gaugefield = std::vector<Matrixsu3>(numberOfElements);
	}
	
	int getVol4d() {return vol4d;}
	int getNumberOfElements() {return numberOfElements;}
private:
	int ntime;
	int nspace;
	int vol4d;
	int numberOfElements;
	//todo: make this a smart pointer?
	std::vector<Matrixsu3> gaugefield;
};

BOOST_AUTO_TEST_CASE(setZero)
{
	MatrixSu3Tester tester(4,4);
	
	const size_t numberOfElements = tester.getNumberOfElements();
	Matrixsu3 * gaugefield = new Matrixsu3[ numberOfElements ];

	Matrixsu3_utilities::fillMatrixSu3Array_constantMatrix(gaugefield, numberOfElements, ZERO);
	
	hmc_complex sum = Matrixsu3_utilities::sumUpAllMatrixElements(gaugefield, numberOfElements);
	
	BOOST_REQUIRE_EQUAL(sum.re, 0.);
	BOOST_REQUIRE_EQUAL(sum.im, 0.);
}

BOOST_AUTO_TEST_CASE(setOne)
{
	int ntime = 4;
	int nspace = 4;
	int vol4d = ntime * nspace * nspace * nspace;
	
	const size_t numberOfElements = vol4d * 4;
	Matrixsu3 * gaugefield = new Matrixsu3[ numberOfElements ];

	Matrixsu3_utilities::fillMatrixSu3Array_constantMatrix(gaugefield, numberOfElements, ONE);
	
	double expectedResult = numberOfElements * numMatrixEntries;
	
	hmc_complex sum = Matrixsu3_utilities::sumUpAllMatrixElements(gaugefield, numberOfElements);
	
	BOOST_REQUIRE_EQUAL(expectedResult, sum.re);
	BOOST_REQUIRE_EQUAL(expectedResult, sum.im);
}

