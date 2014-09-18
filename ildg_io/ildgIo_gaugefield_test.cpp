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

int expectedPrecision = 64;
std::string nameOfExistingGaugefieldFile = std::string(SOURCEDIR) + "/ildg_io/conf.00200";

void checkMetadataOfSpecificGaugefieldFile(sourcefileparameters_values toCheck)
{
  BOOST_REQUIRE_EQUAL(toCheck.lx, 4);
  BOOST_REQUIRE_EQUAL(toCheck.ly, 4);
  BOOST_REQUIRE_EQUAL(toCheck.lz, 4);
  BOOST_REQUIRE_EQUAL(toCheck.lt, 4);
  BOOST_REQUIRE_EQUAL(toCheck.prec, 64);
  BOOST_REQUIRE_EQUAL(toCheck.num_entries, 18432);
  BOOST_REQUIRE_EQUAL(toCheck.flavours, 0);
  BOOST_REQUIRE_EQUAL(toCheck.trajectorynr, 200);
  /**
   * Some values (time!) are not properly written to file, or their values are not correctly gathered.
   */
  BOOST_REQUIRE_EQUAL(toCheck.time, -619635472);
  BOOST_REQUIRE_EQUAL(toCheck.time_solver, 0.);
  BOOST_REQUIRE_EQUAL(toCheck.kappa_solver, 0.);
  BOOST_REQUIRE_EQUAL(toCheck.mu_solver, 0);
  BOOST_REQUIRE_EQUAL(toCheck.noiter, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epssq, 0);
  BOOST_REQUIRE_EQUAL(toCheck.plaquettevalue, 0.571077);
  BOOST_REQUIRE_EQUAL(toCheck.beta, 5.69);
  BOOST_REQUIRE_EQUAL(toCheck.kappa, 0.125);
  BOOST_REQUIRE_EQUAL(toCheck.mu, 0.006);
  BOOST_REQUIRE_EQUAL(toCheck.c2_rec, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mubar, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epsilonbar, 0);
}

BOOST_AUTO_TEST_CASE(readInGaugefieldFailureWithFileException)
{
  std::string nameOfNonexistingGaugefieldFile = "thisfileshouldnotbethere";
  char * bufferToStoreGaugefield;
  BOOST_CHECK_THROW(IldgIoReader_gaugefield srcFileParams(nameOfNonexistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield), File_Exception);
}

BOOST_AUTO_TEST_CASE(readInGaugefieldFailureWithWrongPrecision)
{
  char * bufferToStoreGaugefield;
  int wrongPrecision = 27;
  BOOST_CHECK_THROW(IldgIoReader_gaugefield srcFileParams(nameOfExistingGaugefieldFile.c_str(), wrongPrecision, &bufferToStoreGaugefield), std::exception);
}

BOOST_AUTO_TEST_CASE(readInGaugefieldSuccess)
{
  char * bufferToStoreGaugefield;
  BOOST_REQUIRE_NO_THROW(IldgIoReader_gaugefield srcFileParams(nameOfExistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield));
}

BOOST_AUTO_TEST_CASE(readInGaugefieldCheckMetadata)
{
  char * bufferToStoreGaugefield;
  IldgIoReader_gaugefield srcFileParams(nameOfExistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield);
  checkMetadataOfSpecificGaugefieldFile(srcFileParams.parameters);
}

//todo: this test probably will never work!
BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( readInGaugefieldCheckBufferSize, 1 )
BOOST_AUTO_TEST_CASE(readInGaugefieldCheckBufferSize)
{
  char * bufferToStoreGaugefield;
	IldgIoReader_gaugefield srcFileParams(nameOfExistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield);
  size_t expectedSizeOfBuffer = srcFileParams.parameters.num_entries * sizeof(hmc_float);
  size_t actualSizeOfBuffer = sizeof(bufferToStoreGaugefield);
	BOOST_CHECK_EQUAL(expectedSizeOfBuffer, actualSizeOfBuffer);
}

BOOST_AUTO_TEST_CASE(readInGaugefieldCheckChecksum)
{
  char * bufferToStoreGaugefield;
  uint32_t referenceChecksumA = 171641288;
  uint32_t referenceChecksumB = 3618036129;
  IldgIoReader_gaugefield srcFileParams(nameOfExistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield);
  Checksum referenceChecksum(referenceChecksumA, referenceChecksumB);
  BOOST_REQUIRE_EQUAL(referenceChecksum == srcFileParams.parameters.checksum, true);
}

std::string nameOfExistingGaugefieldFileFromTmlqcd = std::string(SOURCEDIR) + "/ildg_io/conf.tmlqcd";

BOOST_AUTO_TEST_CASE(readInGaugefieldFromTmlqcd_CheckChecksum)
{
  char * bufferToStoreGaugefield;
  uint32_t referenceChecksumA = 398012545;
  uint32_t referenceChecksumB = 1610757546;
  IldgIoReader_gaugefield srcFileParams(nameOfExistingGaugefieldFileFromTmlqcd.c_str(), expectedPrecision, &bufferToStoreGaugefield);
  Checksum referenceChecksum(referenceChecksumA, referenceChecksumB);
  BOOST_REQUIRE_EQUAL(referenceChecksum == srcFileParams.parameters.checksum, true);
}

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

sourcefileparameters_values setSourceFileParametersToSpecificValuesForGaugefield()
{
	sourcefileparameters_values srcFileParams;
	srcFileParams.lx = 3;
	srcFileParams.ly = 3;
	srcFileParams.lz = 3;
	srcFileParams.lt = 5;
	
	srcFileParams.prec = getPrecisionOfDoubleInBits();
	srcFileParams.num_entries = getElementsOfGaugefield(
		srcFileParams.lx, srcFileParams.ly, srcFileParams.lz, srcFileParams.lt);
	srcFileParams.flavours = 33;
	srcFileParams.trajectorynr = 1234567890;
	srcFileParams.time = -12345;
	srcFileParams.time_solver = -12345;
	srcFileParams.noiter = -6789;
	srcFileParams.plaquettevalue = -12.34567833;
	srcFileParams.beta = -12.345678;
	srcFileParams.kappa = -56.7890;
	srcFileParams.mu = -56.7890;
	srcFileParams.c2_rec = -13.579;
	srcFileParams.mubar = -2.4680;
	srcFileParams.epsilonbar = -19.8273;
	srcFileParams.epssq = -19.87654;
	srcFileParams.kappa_solver = -45.6784;
	srcFileParams.mu_solver = -10.9283;
	
	srcFileParams.field = getFieldType_gaugefield();
	srcFileParams.date = "someDate";
	srcFileParams.hmcversion = "0.0";
	srcFileParams.solvertype = "someSolver";
	srcFileParams.hmcversion_solver = "notImplemented";
	srcFileParams.date_solver = "someDate";
	
	return srcFileParams;
}

// one cannot expect that date, time and time_solver will match..
// not implemented or fermion parameters: solvertype_source, hmcversion_solver_source, flavours_source, noiter_source, kappa_solver_source, mu_solver_source, epssq_source
BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( writeGaugefield_metaData, 11 )

void compareTwoSourcefileParameters(sourcefileparameters_values toCheck1, sourcefileparameters_values toCheck2)
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
	BOOST_REQUIRE_EQUAL(toCheck1.hmcversion, toCheck2.hmcversion);
	BOOST_CHECK_EQUAL(toCheck1.hmcversion_solver, toCheck2.hmcversion_solver);
	BOOST_CHECK_EQUAL(toCheck1.date_solver, toCheck2.date_solver);
	BOOST_CHECK_EQUAL(toCheck1.solvertype, toCheck2.solvertype);
}

void writeEmptyGaugefieldFromSourcefileParameters(sourcefileparameters_values srcFileParams, std::string configurationName)
{
	const n_uint64_t bufferSize_gaugefield = getElementsOfGaugefield(srcFileParams.lx, srcFileParams.ly, srcFileParams.lz, srcFileParams.lt) * sizeof(double);
	
	char * binaryData = new char[ bufferSize_gaugefield ];
	
	//TODO: hmc version currently can not be anything else than #.# !!
	IldgIoWriter_gaugefield writer( binaryData, bufferSize_gaugefield, srcFileParams ,configurationName);
	
	delete binaryData;
}

BOOST_AUTO_TEST_CASE(writeGaugefield_metaData)
{
	sourcefileparameters_values srcFileParams_1 = setSourceFileParametersToSpecificValuesForGaugefield();
	
	//TODO: test with single?
	std::string configurationName = "conf.test";
	
	writeEmptyGaugefieldFromSourcefileParameters(srcFileParams_1, configurationName);
	
	char * readBinaryData;
	IldgIoReader_gaugefield srcFileParams_2(configurationName.c_str(), srcFileParams_1.prec, &readBinaryData);
	delete readBinaryData;
	
	compareTwoSourcefileParameters(srcFileParams_1, srcFileParams_2.parameters);
}
