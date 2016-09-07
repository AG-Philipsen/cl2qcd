/** @file
 * Testcases for the ildg I/O class
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
#define BOOST_TEST_MODULE ildg_read_gaugefield
#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>

#include "ildgIo_gaugefield.hpp"
#include "../meta/util.hpp"

using namespace ildgIo;

// TODO: Remove all occurences of meta::Inputparameters and Inputparameters. The tests should only use IlgdIoParameters_gaugefield!!

size_t getPrecisionOfDoubleInBits()
{
	return 8 * sizeof(double);
}

size_t getElementsOfGaugefield(int lx, int ly, int lz, int lt)
{
	const size_t spatialVolume = lx * ly * lz;
	return 2 * NC * NC * NDIM * spatialVolume * lt;
}

//todo: put this into a common file
std::string getFieldType_gaugefield()
{
	return "su3gauge";
}

Sourcefileparameters setSourceFileParametersToSpecificValuesForGaugefield()
{
	Sourcefileparameters srcFileParams;
	srcFileParams.lx = 3;
	srcFileParams.ly = 3;
	srcFileParams.lz = 3;
	srcFileParams.lt = 5;
	
	srcFileParams.prec = getPrecisionOfDoubleInBits();
	srcFileParams.num_entries = getElementsOfGaugefield(
		srcFileParams.lx, srcFileParams.ly, srcFileParams.lz, srcFileParams.lt);
	srcFileParams.trajectorynr = 1234567890;
	srcFileParams.time = -12345;
	srcFileParams.time_solver = -12345;
	srcFileParams.plaquettevalue = -12.34567833;
	srcFileParams.beta = -12.345678;
	srcFileParams.kappa = -56.7890;
	srcFileParams.mu = -56.7890;
	
	srcFileParams.field = getFieldType_gaugefield();
	srcFileParams.date = "someDate";
	srcFileParams.hmcversion = "0.0";
	srcFileParams.solvertype = "someSolver";
	srcFileParams.hmcversion_solver = "notImplemented";
	srcFileParams.date_solver = "someDate";
	
	return srcFileParams;
}

// one cannot expect that hmcVersion, date, time and time_solver will match..
// not implemented are fermion parameters: solvertype_source, hmcversion_solver_source, flavours_source, noiter_source, kappa_solver_source, mu_solver_source, epssq_source
BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( writeGaugefield_metaData, 12 )

void compareTwoSourcefileParameters(Sourcefileparameters toCheck1, Sourcefileparameters toCheck2)
{
  BOOST_REQUIRE_EQUAL(toCheck1.lx, toCheck2.lx);
  BOOST_REQUIRE_EQUAL(toCheck1.ly, toCheck2.ly);
  BOOST_REQUIRE_EQUAL(toCheck1.lz, toCheck2.lz);
  BOOST_REQUIRE_EQUAL(toCheck1.lt, toCheck2.lt);
  BOOST_REQUIRE_EQUAL(toCheck1.prec, toCheck2.prec);
  BOOST_REQUIRE_EQUAL(toCheck1.num_entries, toCheck2.num_entries);
  BOOST_CHECK_EQUAL(toCheck1.flavours, toCheck2.flavours);
  BOOST_REQUIRE_EQUAL(toCheck1.trajectorynr, toCheck2.trajectorynr);
  BOOST_CHECK_EQUAL(toCheck1.time, toCheck2.time);
  BOOST_CHECK_EQUAL(toCheck1.time_solver, toCheck2.time_solver);
  BOOST_CHECK_EQUAL(toCheck1.kappa_solver, toCheck2.kappa_solver);
  BOOST_CHECK_EQUAL(toCheck1.mu_solver, toCheck2.mu_solver);
  BOOST_CHECK_EQUAL(toCheck1.noiter, toCheck2.noiter);
  BOOST_CHECK_EQUAL(toCheck1.epssq, toCheck2.epssq);
  BOOST_REQUIRE_EQUAL(toCheck1.plaquettevalue, toCheck2.plaquettevalue);
  BOOST_REQUIRE_EQUAL(toCheck1.beta, toCheck2.beta);
  BOOST_REQUIRE_EQUAL(toCheck1.kappa, toCheck2.kappa);
  BOOST_REQUIRE_EQUAL(toCheck1.mu, toCheck2.mu);
  BOOST_REQUIRE_EQUAL(toCheck1.c2_rec, toCheck2.c2_rec);
  BOOST_REQUIRE_EQUAL(toCheck1.mubar, toCheck2.mubar);
  BOOST_REQUIRE_EQUAL(toCheck1.epsilonbar, toCheck2.epsilonbar);
	
	BOOST_REQUIRE_EQUAL(toCheck1.field, toCheck2.field);
	BOOST_CHECK_EQUAL(toCheck1.date, toCheck2.date);
	BOOST_CHECK_EQUAL(toCheck1.hmcversion, toCheck2.hmcversion);
	BOOST_CHECK_EQUAL(toCheck1.hmcversion_solver, toCheck2.hmcversion_solver);
	BOOST_CHECK_EQUAL(toCheck1.date_solver, toCheck2.date_solver);
	BOOST_CHECK_EQUAL(toCheck1.solvertype, toCheck2.solvertype);
}

void writeEmptyGaugefieldFromSourcefileParameters(const meta::Inputparameters * parameters, std::string configurationName)
{
	const n_uint64_t numberElements = getNumberOfElements_gaugefield(parameters);
	std::vector<Matrixsu3> gaugefield(numberElements);
	
	const physics::lattices::GaugefieldParametersImplementation tmp (parameters);
	Inputparameters test2( &tmp );
	IldgIoParameters_gaugefield test(&test2);
	IldgIoWriter_gaugefield writer( gaugefield, &test ,configurationName, 1234567890, -12.34567833);
}

BOOST_AUTO_TEST_CASE(writeGaugefield_metaData)
{
	Sourcefileparameters srcFileParams_1 = setSourceFileParametersToSpecificValuesForGaugefield();
	//do not consider the checksum in this test
	const char * tmp [] = {"foo", "--nspace=3", "--ntime=5", "--ignore_checksum_errors=true", "--beta=-12.345678", "--kappa=-56.7890", "--mu=-56.7890"};
	const meta::Inputparameters parameters(7, tmp);
	//TODO: test with single?
	std::string configurationName = "conf.test";
	
	writeEmptyGaugefieldFromSourcefileParameters(&parameters, configurationName);
	
	Matrixsu3 * readBinaryData = nullptr;
	const physics::lattices::GaugefieldParametersImplementation tmp2 (&parameters);
	Inputparameters test2( &tmp2 );
	IldgIoParameters_gaugefield test(&test2);
	IldgIoReader_gaugefield readGaugefield(configurationName, &test, &readBinaryData);
	delete readBinaryData;
	
	compareTwoSourcefileParameters(srcFileParams_1, readGaugefield.parameters);
}

#include "matrixSu3_utilities.hpp"

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

void checkSumForOneDiagonalMatrix(hmc_complex sum)
{
	BOOST_REQUIRE_EQUAL(sum.re, 3.);
	BOOST_REQUIRE_EQUAL(sum.im, 0.);
}

void checkSumForFilledMatrix(hmc_complex sum)
{
	BOOST_REQUIRE_EQUAL(sum.re, 45.);
	BOOST_REQUIRE_EQUAL(sum.im, 45.);
}

class MatrixSu3Field
{
public:
	MatrixSu3Field(uint ntime, uint nspace)
	{
		std::string nspaceString = "--nspace=" + boost::lexical_cast<std::string>(nspace);
		std::string ntimeString = "--ntime="  + boost::lexical_cast<std::string>(ntime);
		const char * tmp [] = {"foo", nspaceString.c_str(), ntimeString.c_str()};
		parameters = new meta::Inputparameters(3, tmp);
		numberOfElements = getNumberOfElements_gaugefield(parameters);
		gaugefield = std::vector<Matrixsu3>(numberOfElements);
		Matrixsu3_utilities::fillMatrixSu3Array_constantMatrix(gaugefield, Matrixsu3_utilities::FillType::ZERO);
	}
	
	n_uint64_t getNumberOfElements() { return numberOfElements; }
	void setSpecificEntry(Matrixsu3 in, int positionToSet)
	{
		if ( ((int) numberOfElements - positionToSet) < 0)
			{
				logger.fatal() << "Got bad coordinates to set matrix. Aborting...";
				BOOST_REQUIRE_EQUAL(0,1);
			}
			gaugefield[positionToSet] = in;
	}
	
	const std::vector<Matrixsu3> &getField() {return gaugefield; }
	const meta::Inputparameters * getParameters() {return parameters;}
	Matrixsu3 getEntry(int position) { return gaugefield[position];}
	Matrixsu3 * getPointerToField() { return &gaugefield[0]; }
	void setField(Matrixsu3 * in) { gaugefield.assign(in, in + getNumberOfElements()); }
	void fillField(Matrixsu3_utilities::FillType fillType) 
	{ 
		if (fillType == Matrixsu3_utilities::RANDOM)
		{
			Matrixsu3_utilities::fillMatrixSu3Array_randomMatrix(gaugefield);
		}
		else
		{
			Matrixsu3_utilities::fillMatrixSu3Array_constantMatrix(gaugefield, fillType);
		}
	}
private:
	const meta::Inputparameters * parameters;
	n_uint64_t numberOfElements;
	std::vector<Matrixsu3> gaugefield;
};

BOOST_AUTO_TEST_SUITE(conversionToAndFromIldgFormat)

	void convertGaugefieldToAndFromIldg(MatrixSu3Field & in)
	{
		n_uint64_t num_bytes = getSizeInBytes_gaugefield(in.getNumberOfElements());
		std::vector<char> binary_data(num_bytes);
		char * binary_data_ptr = &binary_data[0];

		Matrixsu3 * gaugefieldTmp = in.getPointerToField();

		const physics::lattices::GaugefieldParametersImplementation tmp (in.getParameters());
		Inputparameters test2( &tmp );
		IldgIoParameters_gaugefield test(&test2);
		copy_gaugefield_to_ildg_format(binary_data, in.getField(), test  );
		copy_gaugefield_from_ildg_format(gaugefieldTmp, binary_data_ptr, in.getNumberOfElements() * 9 * 2, test );

		in.setField(gaugefieldTmp);
	}

	void fillAndConvert(Matrixsu3_utilities::FillType fillType)
	{
		MatrixSu3Field gaugefield(4,4);

		gaugefield.fillField( fillType );

		hmc_complex sumBeforeConversion = Matrixsu3_utilities::sumUpAllMatrixElements(gaugefield.getField());

		convertGaugefieldToAndFromIldg(gaugefield);

		hmc_complex sumAfterConversion = Matrixsu3_utilities::sumUpAllMatrixElements(gaugefield.getField());

		BOOST_REQUIRE_EQUAL(sumBeforeConversion.re, sumAfterConversion.re);
		BOOST_REQUIRE_EQUAL(sumBeforeConversion.im, sumAfterConversion.im);
	}

	BOOST_AUTO_TEST_CASE(conversionToAndFromIldgFormat_global_constant)
	{
		fillAndConvert(Matrixsu3_utilities::ONE);
	}

	BOOST_AUTO_TEST_CASE(conversionToAndFromIldgFormat_global_random)
	{
		fillAndConvert(Matrixsu3_utilities::RANDOM);
	}

	void testConversion(uint nspace, uint ntime)
	{
		uint positionToSet = 0;
		MatrixSu3Field gaugefield(ntime, nspace);
		uint sitesVisited = 0;

		for(uint x = 0; x<nspace; x++)
		{
			for(uint y = 0; y<nspace; y++)
			{
				for(uint z = 0; z<nspace; z++)
				{
					for(uint t = 0; t<ntime; t++)
					{
						for(uint mu = 0; mu < NDIM; mu++)
						{
							positionToSet = uint(LinkIndex(Index(x,y,z,t,LatticeExtents(nspace,ntime)),static_cast<Direction>(mu)));
							gaugefield.setSpecificEntry(Matrixsu3_utilities::getUnitMatrix(), positionToSet);

							hmc_complex sumBeforeConversion = Matrixsu3_utilities::sumUpAllMatrixElements( gaugefield.getField() );
							checkSumForOneDiagonalMatrix(sumBeforeConversion);

							convertGaugefieldToAndFromIldg(gaugefield);

							hmc_complex sumAfterConversion = Matrixsu3_utilities::sumUpAllMatrixElements(gaugefield.getField());

							checkSumForOneDiagonalMatrix(sumAfterConversion);

							Matrixsu3 set = gaugefield.getEntry(positionToSet);
							checkMatrixSu3ForDiagonalType(set);

							gaugefield.setSpecificEntry(Matrixsu3_utilities::getZeroMatrix(), positionToSet);

							sitesVisited ++;
						}
					}
				}
			}
		}

		BOOST_REQUIRE_EQUAL(sitesVisited, gaugefield.getNumberOfElements() );
	}

	BOOST_AUTO_TEST_CASE(conversionToAndFromIldgFormat_specific_1)
	{
		testConversion(3,3);
	}

	BOOST_AUTO_TEST_CASE(conversionToAndFromIldgFormat_specific_2)
	{
		testConversion(3,5);
	}

	BOOST_AUTO_TEST_CASE(conversionToAndFromIldgFormat_specific_3)
	{
		testConversion(5,3);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(writeAndRead)

	void writeFieldToFile(MatrixSu3Field &in, std::string filename)
	{
		const physics::lattices::GaugefieldParametersImplementation tmp (in.getParameters());
		Inputparameters test2( &tmp );
		IldgIoParameters_gaugefield test(&test2);
		IldgIoWriter_gaugefield writer( in.getField(), &test , filename, 0, 0.);
	}

	void readFieldFromFile(MatrixSu3Field &out, std::string filename)
	{
		Matrixsu3 * gaugefieldTmp = NULL;
		const physics::lattices::GaugefieldParametersImplementation tmp ( out.getParameters());
		Inputparameters test2( &tmp );
		IldgIoParameters_gaugefield test(&test2);
		IldgIoReader_gaugefield reader(filename, &test, &gaugefieldTmp);
		out.setField(gaugefieldTmp);
	}

	void testWriteRead(uint nspace, uint ntime)
	{
		uint positionToSet = 0;
		MatrixSu3Field gaugefield(ntime, nspace);
		std::string filename = "conf.test";
		uint sitesVisited = 0;

		for(uint x = 0; x<nspace; x++)
		{
			for(uint y = 0; y<nspace; y++)
			{
				for(uint z = 0; z<nspace; z++)
				{
					for(uint t = 0; t<ntime; t++)
					{
						for(uint mu = 0; mu < NDIM; mu++)
						{
							positionToSet = uint(LinkIndex(Index(x,y,z,t,LatticeExtents(nspace,ntime)),static_cast<Direction>(mu)));
							gaugefield.setSpecificEntry(Matrixsu3_utilities::getFilledMatrix(), positionToSet);

							hmc_complex sumBeforeConversion = Matrixsu3_utilities::sumUpAllMatrixElements( gaugefield.getField() );
							checkSumForFilledMatrix(sumBeforeConversion);

							writeFieldToFile(gaugefield, filename);
							gaugefield.setSpecificEntry(Matrixsu3_utilities::getZeroMatrix(), positionToSet);

							readFieldFromFile(gaugefield, filename);

							hmc_complex sumAfterConversion = Matrixsu3_utilities::sumUpAllMatrixElements(gaugefield.getField());
							checkSumForFilledMatrix(sumAfterConversion);

							Matrixsu3 set = gaugefield.getEntry(positionToSet);
							checkMatrixSu3ForFilledType(set);

							gaugefield.setSpecificEntry(Matrixsu3_utilities::getZeroMatrix(), positionToSet);

							sitesVisited ++;
						}
					}
				}
			}
		}

		BOOST_REQUIRE_EQUAL(sitesVisited, gaugefield.getNumberOfElements() );
	}

	BOOST_AUTO_TEST_CASE(writeAndRead_1)
	{
		testWriteRead(2,3);
	}
	BOOST_AUTO_TEST_CASE(writeAndRead_2)
	{
		testWriteRead(3,2);
	}
	BOOST_AUTO_TEST_CASE(writeAndRead_3)
	{
		testWriteRead(2,1);
	}

BOOST_AUTO_TEST_SUITE_END()
