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

#include "ildgIo_gaugefield.hpp"
#include "../meta/util.hpp"

int expectedPrecision = 64;
std::string nameOfExistingGaugefieldFile = std::string(SOURCEDIR) + "/ildg_io/conf.00200";

using namespace ildgIo;

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
// not implemented or fermion parameters: solvertype_source, hmcversion_solver_source, flavours_source, noiter_source, kappa_solver_source, mu_solver_source, epssq_source
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
	Matrixsu3 * gaugefield = new Matrixsu3[ numberElements ];
	
	IldgIoWriter_gaugefield writer( gaugefield, parameters ,configurationName, 1234567890, -12.34567833);
	
	delete gaugefield;
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
	IldgIoReader_gaugefield readGaugefield(configurationName, &parameters, &readBinaryData);
	delete readBinaryData;
	
	compareTwoSourcefileParameters(srcFileParams_1, readGaugefield.parameters);
}

#include "matrixSu3_utilities.hpp"

BOOST_AUTO_TEST_CASE(conversionToAndFromIldgFormat)
{
	const char * tmp [] = {"foo", "--nspace=4", "--ntime=4"};
	const meta::Inputparameters parameters(3, tmp);
	
	const n_uint64_t numberOfElements = getNumberOfElements_gaugefield(&parameters);
	Matrixsu3 * gaugefield = new Matrixsu3[ numberOfElements ];
	Matrixsu3_utilities::fillMatrixSu3Array_constantMatrix(gaugefield, numberOfElements, Matrixsu3_utilities::FillType::ONE);
	
	hmc_complex sumBeforeConversion = Matrixsu3_utilities::sumUpAllMatrixElements(gaugefield, numberOfElements);
	
	n_uint64_t num_bytes = getSizeInBytes_gaugefield(numberOfElements);
	char * binary_data = new char[num_bytes];
	
	copy_gaugefield_to_ildg_format(binary_data, gaugefield, parameters);
	copy_gaugefield_from_ildg_format(gaugefield, binary_data, numberOfElements * 9 * 2, parameters);
	
	hmc_complex sumAfterConversion = Matrixsu3_utilities::sumUpAllMatrixElements(gaugefield, numberOfElements);
	
	BOOST_REQUIRE_EQUAL(sumBeforeConversion.re, sumAfterConversion.re);
	BOOST_REQUIRE_EQUAL(sumBeforeConversion.im, sumAfterConversion.im);
}
