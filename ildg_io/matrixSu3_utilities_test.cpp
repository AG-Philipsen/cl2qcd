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
	MatrixSu3Tester(int ntimeIn, int nspaceIn, Matrixsu3_utilities::FillType filltype = ZERO) : 
		ntime(ntimeIn), nspace(nspaceIn)
	{
		vol4d = ntime * nspace * nspace * nspace;
		numberOfElements = vol4d * 4;
		gaugefield = std::vector<Matrixsu3>(numberOfElements);
		
		if (filltype == Matrixsu3_utilities::RANDOM)
		{
			Matrixsu3_utilities::fillMatrixSu3Array_randomMatrix(gaugefield);
		}
		else
		{
			Matrixsu3_utilities::fillMatrixSu3Array_constantMatrix(gaugefield, filltype);
		}
	}
	
	int getVol4d() {return vol4d;}
	int getNumberOfElements() {return numberOfElements;}
	
	const std::vector<Matrixsu3> &getGaugefield() const {return gaugefield;}
	
	Matrixsu3 getEntry(int position){ return gaugefield[position]; }
protected:
	int ntime;
	int nspace;
	int vol4d;
	int numberOfElements;
	//todo: make this a smart pointer?
	std::vector<Matrixsu3> gaugefield;
};

#include "../executables/exceptions.h"

#include <iostream>

class MatrixSu3Specific : public MatrixSu3Tester
{
public:
	MatrixSu3Specific(int ntime, int nspace, int position, FillType filltype = DIAGONAL) : MatrixSu3Tester(ntime, nspace, ZERO)
	{
		if ( position > numberOfElements)
		{
			throw std::logic_error("Got invalid argument to set specific MatrixSu3 entry. Aborting...");
		}
		if (filltype == FILLED)
		{
			gaugefield[position] = getFilledMatrix();
		}
		else
		{
			gaugefield[position] = getUnitMatrix();
		}
	}
};

BOOST_AUTO_TEST_CASE(setZero)
{
	MatrixSu3Tester tester(7,23);
	
	hmc_complex sum = Matrixsu3_utilities::sumUpAllMatrixElements(tester.getGaugefield() );
	
	BOOST_REQUIRE_EQUAL(sum.re, 0.);
	BOOST_REQUIRE_EQUAL(sum.im, 0.);
}

BOOST_AUTO_TEST_CASE(setOne)
{
	MatrixSu3Tester tester(15,31, ONE);

	hmc_complex sum = Matrixsu3_utilities::sumUpAllMatrixElements(tester.getGaugefield() );
	double expectedResult = tester.getNumberOfElements() * numMatrixEntries;
	
	BOOST_REQUIRE_EQUAL(expectedResult, sum.re);
	BOOST_REQUIRE_EQUAL(expectedResult, sum.im);
}

BOOST_AUTO_TEST_CASE(setFilled)
{
	MatrixSu3Tester tester(15,31, FILLED);

	hmc_complex sum = Matrixsu3_utilities::sumUpAllMatrixElements(tester.getGaugefield() );
	double expectedResult = tester.getNumberOfElements() * ( 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9);
	
	BOOST_REQUIRE_EQUAL(expectedResult, sum.re);
	BOOST_REQUIRE_EQUAL(expectedResult, sum.im);
}

BOOST_AUTO_TEST_CASE(setDiagonal)
{
	MatrixSu3Tester tester(3,8, DIAGONAL);

	hmc_complex sum = Matrixsu3_utilities::sumUpDiagonalMatrixElements(tester.getGaugefield() );
	double expectedResultForRealPart = tester.getNumberOfElements() * 3;
	
	BOOST_REQUIRE_EQUAL(expectedResultForRealPart, sum.re);
	BOOST_REQUIRE_EQUAL(0., sum.im);
	
	sum = Matrixsu3_utilities::sumUpOffDiagonalMatrixElements(tester.getGaugefield() );	
	
	BOOST_REQUIRE_EQUAL(0., sum.re);
	BOOST_REQUIRE_EQUAL(0., sum.im);	
}

BOOST_AUTO_TEST_CASE(setSpecific_throw)
{
	int positionToSet = 1000;
	BOOST_REQUIRE_THROW(MatrixSu3Specific tester(3,3, positionToSet), std::logic_error );
}

void checkMatrixSu3ForDiagonalType(Matrixsu3 in)
{
	BOOST_REQUIRE_EQUAL(1., in.e00.re);
	BOOST_REQUIRE_EQUAL(0., in.e01.re);
	BOOST_REQUIRE_EQUAL(0., in.e02.re);
	BOOST_REQUIRE_EQUAL(0., in.e10.re);
	BOOST_REQUIRE_EQUAL(1., in.e11.re);
	BOOST_REQUIRE_EQUAL(0., in.e12.re);
	BOOST_REQUIRE_EQUAL(0., in.e20.re);
	BOOST_REQUIRE_EQUAL(0., in.e21.re);
	BOOST_REQUIRE_EQUAL(1., in.e22.re);
	
	BOOST_REQUIRE_EQUAL(0., in.e00.im);
	BOOST_REQUIRE_EQUAL(0., in.e01.im);
	BOOST_REQUIRE_EQUAL(0., in.e02.im);
	BOOST_REQUIRE_EQUAL(0., in.e10.im);
	BOOST_REQUIRE_EQUAL(0., in.e11.im);
	BOOST_REQUIRE_EQUAL(0., in.e12.im);
	BOOST_REQUIRE_EQUAL(0., in.e20.im);
	BOOST_REQUIRE_EQUAL(0., in.e21.im);
	BOOST_REQUIRE_EQUAL(0., in.e22.im);	
}

void checkMatrixSu3ForFilledType(Matrixsu3 in)
{
	BOOST_REQUIRE_EQUAL(1., in.e00.re);
	BOOST_REQUIRE_EQUAL(2., in.e01.re);
	BOOST_REQUIRE_EQUAL(3., in.e02.re);
	BOOST_REQUIRE_EQUAL(4., in.e10.re);
	BOOST_REQUIRE_EQUAL(5., in.e11.re);
	BOOST_REQUIRE_EQUAL(6., in.e12.re);
	BOOST_REQUIRE_EQUAL(7., in.e20.re);
	BOOST_REQUIRE_EQUAL(8., in.e21.re);
	BOOST_REQUIRE_EQUAL(9., in.e22.re);
	
	BOOST_REQUIRE_EQUAL(1., in.e00.im);
	BOOST_REQUIRE_EQUAL(2., in.e01.im);
	BOOST_REQUIRE_EQUAL(3., in.e02.im);
	BOOST_REQUIRE_EQUAL(4., in.e10.im);
	BOOST_REQUIRE_EQUAL(5., in.e11.im);
	BOOST_REQUIRE_EQUAL(6., in.e12.im);
	BOOST_REQUIRE_EQUAL(7., in.e20.im);
	BOOST_REQUIRE_EQUAL(8., in.e21.im);
	BOOST_REQUIRE_EQUAL(9., in.e22.im);	
}

BOOST_AUTO_TEST_CASE(setSpecific_diagonal)
{
	int positionToSet = 24;
	
	MatrixSu3Specific tester(21,9, positionToSet);

	hmc_complex sum = Matrixsu3_utilities::sumUpAllMatrixElements(tester.getGaugefield() );
	double expectedResultForRealPart = 3;
	
	BOOST_REQUIRE_EQUAL(expectedResultForRealPart, sum.re);
	BOOST_REQUIRE_EQUAL(0., sum.im);
	
	Matrixsu3 set = tester.getEntry(positionToSet);
	checkMatrixSu3ForDiagonalType(set);	
}

BOOST_AUTO_TEST_CASE(setSpecific_filled)
{
	int positionToSet = 13;
	
	MatrixSu3Specific tester(6,5, positionToSet, FILLED);

	hmc_complex sum = Matrixsu3_utilities::sumUpAllMatrixElements(tester.getGaugefield() );
	double expectedResult = ( 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9);
	
	BOOST_REQUIRE_EQUAL(expectedResult, sum.re);
	BOOST_REQUIRE_EQUAL(expectedResult, sum.im);
	
	Matrixsu3 set = tester.getEntry(positionToSet);
	checkMatrixSu3ForFilledType(set);	
}

