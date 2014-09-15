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

void checkDefaults(sourcefileparameters toCheck)
{
  BOOST_REQUIRE_EQUAL(toCheck.lx_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.ly_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.lz_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.lt_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.prec_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.num_entries_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.flavours_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.trajectorynr_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.time_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.time_solver_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.noiter_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.plaquettevalue_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.beta_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.kappa_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mu_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.c2_rec_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mubar_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epsilonbar_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epssq_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.kappa_solver_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mu_solver_source, 0);
}

void checkMetadataOfSpecificGaugefieldFile(sourcefileparameters toCheck)
{
  BOOST_REQUIRE_EQUAL(toCheck.lx_source, 4);
  BOOST_REQUIRE_EQUAL(toCheck.ly_source, 4);
  BOOST_REQUIRE_EQUAL(toCheck.lz_source, 4);
  BOOST_REQUIRE_EQUAL(toCheck.lt_source, 4);
  BOOST_REQUIRE_EQUAL(toCheck.prec_source, 64);
  BOOST_REQUIRE_EQUAL(toCheck.num_entries_source, 18432);
  BOOST_REQUIRE_EQUAL(toCheck.flavours_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.trajectorynr_source, 200);
  /**
   * Some values (time!) are not properly written to file, or their values are not correctly gathered.
   * TODO: repair!
   */
  BOOST_REQUIRE_EQUAL(toCheck.time_source, -619635472);
  BOOST_REQUIRE_EQUAL(toCheck.time_solver_source, 0.);
  BOOST_REQUIRE_EQUAL(toCheck.kappa_solver_source, 0.);
  BOOST_REQUIRE_EQUAL(toCheck.mu_solver_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.noiter_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epssq_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.plaquettevalue_source, 0.571077);
  BOOST_REQUIRE_EQUAL(toCheck.beta_source, 5.69);
  BOOST_REQUIRE_EQUAL(toCheck.kappa_source, 0.125);
  BOOST_REQUIRE_EQUAL(toCheck.mu_source, 0.006);
  BOOST_REQUIRE_EQUAL(toCheck.c2_rec_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mubar_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epsilonbar_source, 0);
}

BOOST_AUTO_TEST_CASE(defaults)
{
  sourcefileparameters srcFileParams;
  checkDefaults(srcFileParams);
}

BOOST_AUTO_TEST_CASE(readInGaugefieldFailureWithFileException)
{
  sourcefileparameters srcFileParams;
  std::string nameOfNonexistingGaugefieldFile = "thisfileshouldnotbethere";
  char * bufferToStoreGaugefield;
  BOOST_CHECK_THROW(srcFileParams.readsourcefile(nameOfNonexistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield), File_Exception);
}

BOOST_AUTO_TEST_CASE(readInGaugefieldFailureWithWrongPrecision)
{
  sourcefileparameters srcFileParams;
  char * bufferToStoreGaugefield;
  int wrongPrecision = 27;
  BOOST_CHECK_THROW(srcFileParams.readsourcefile(nameOfExistingGaugefieldFile.c_str(), wrongPrecision, &bufferToStoreGaugefield), std::exception);
}

BOOST_AUTO_TEST_CASE(readInGaugefieldSuccess)
{
  sourcefileparameters srcFileParams;
  char * bufferToStoreGaugefield;
  BOOST_REQUIRE_NO_THROW(srcFileParams.readsourcefile(nameOfExistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield));
}

BOOST_AUTO_TEST_CASE(readInGaugefieldCheckMetadata)
{
  sourcefileparameters srcFileParams;
  char * bufferToStoreGaugefield;
  srcFileParams.readsourcefile(nameOfExistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield);
  checkMetadataOfSpecificGaugefieldFile(srcFileParams);
}

//todo: this test probably will never work!
BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( readInGaugefieldCheckBufferSize, 1 )
BOOST_AUTO_TEST_CASE(readInGaugefieldCheckBufferSize)
{
  sourcefileparameters srcFileParams;
  char * bufferToStoreGaugefield;
	srcFileParams.readsourcefile(nameOfExistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield);
  size_t expectedSizeOfBuffer = srcFileParams.num_entries_source * sizeof(hmc_float);
  size_t actualSizeOfBuffer = sizeof(bufferToStoreGaugefield);
	BOOST_CHECK_EQUAL(expectedSizeOfBuffer, actualSizeOfBuffer);
}

BOOST_AUTO_TEST_CASE(readInGaugefieldCheckChecksum)
{
  sourcefileparameters srcFileParams;
  char * bufferToStoreGaugefield;
  uint32_t referenceChecksumA = 171641288;
  uint32_t referenceChecksumB = 3618036129;
  srcFileParams.readsourcefile(nameOfExistingGaugefieldFile.c_str(), expectedPrecision, &bufferToStoreGaugefield);
  Checksum referenceChecksum(referenceChecksumA, referenceChecksumB);
  BOOST_REQUIRE_EQUAL(referenceChecksum == srcFileParams.checksum, true);
}

std::string nameOfExistingGaugefieldFileFromTmlqcd = std::string(SOURCEDIR) + "/ildg_io/conf.tmlqcd";

BOOST_AUTO_TEST_CASE(readInGaugefieldFromTmlqcd_CheckChecksum)
{
  sourcefileparameters srcFileParams;
  char * bufferToStoreGaugefield;
  uint32_t referenceChecksumA = 398012545;
  uint32_t referenceChecksumB = 1610757546;
  srcFileParams.readsourcefile(nameOfExistingGaugefieldFileFromTmlqcd.c_str(), expectedPrecision, &bufferToStoreGaugefield);
  Checksum referenceChecksum(referenceChecksumA, referenceChecksumB);
  BOOST_REQUIRE_EQUAL(referenceChecksum == srcFileParams.checksum, true);
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

sourcefileparameters setSourceFileParametersToSpecificValuesForGaugefield()
{
	sourcefileparameters srcFileParams;
	srcFileParams.lx_source = 2;
	srcFileParams.ly_source = 2;
	srcFileParams.lz_source = 2;
	srcFileParams.lt_source = 2;
	
	srcFileParams.prec_source = getPrecisionOfDoubleInBits();
	srcFileParams.num_entries_source = getElementsOfGaugefield(
		srcFileParams.lx_source, srcFileParams.ly_source, srcFileParams.lz_source, srcFileParams.lt_source);
	srcFileParams.flavours_source = 0;
	srcFileParams.trajectorynr_source = 0;
	srcFileParams.time_source = 0;
	srcFileParams.time_solver_source = 0;
	srcFileParams.noiter_source = 0;
	srcFileParams.plaquettevalue_source = 0;
	srcFileParams.beta_source = 0;
	srcFileParams.kappa_source = 0;
	srcFileParams.mu_source = 0;
	srcFileParams.c2_rec_source = 0;
	srcFileParams.mubar_source = 0;
	srcFileParams.epsilonbar_source = 0;
	srcFileParams.epssq_source = 0;
	srcFileParams.kappa_solver_source = 0;
	srcFileParams.mu_solver_source = 0;
	
	srcFileParams.field_source = getFieldType_gaugefield();
	srcFileParams.date_source = "someDate";
	srcFileParams.hmcversion_source = "0.0";
	srcFileParams.solvertype_source = "someSolver";
	srcFileParams.hmcversion_solver_source = "notImplemented";
	srcFileParams.date_solver_source = "someDate";
	
	return srcFileParams;
}

// one cannot expect that date, time and time_solver will match..
// not implemented anyway: solvertype_source, hmcversion_solver_source
BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( writeGaugefield, 5 )

void compareTwoSourcefileParameters(sourcefileparameters toCheck1, sourcefileparameters toCheck2)
{
  BOOST_REQUIRE_EQUAL(toCheck1.lx_source, toCheck2.lx_source);
  BOOST_REQUIRE_EQUAL(toCheck1.ly_source, toCheck2.ly_source);
  BOOST_REQUIRE_EQUAL(toCheck1.lz_source, toCheck2.lz_source);
  BOOST_REQUIRE_EQUAL(toCheck1.lt_source, toCheck2.lt_source);
  BOOST_REQUIRE_EQUAL(toCheck1.prec_source, toCheck2.prec_source);
  BOOST_REQUIRE_EQUAL(toCheck1.num_entries_source, toCheck2.num_entries_source);
  BOOST_REQUIRE_EQUAL(toCheck1.flavours_source, toCheck2.flavours_source);
  BOOST_REQUIRE_EQUAL(toCheck1.trajectorynr_source, toCheck2.trajectorynr_source);
  BOOST_CHECK_EQUAL(toCheck1.time_source, toCheck2.time_source);
  BOOST_CHECK_EQUAL(toCheck1.time_solver_source, toCheck2.time_solver_source);
  BOOST_REQUIRE_EQUAL(toCheck1.kappa_solver_source, toCheck2.kappa_solver_source);
  BOOST_REQUIRE_EQUAL(toCheck1.mu_solver_source, toCheck2.mu_solver_source);
  BOOST_REQUIRE_EQUAL(toCheck1.noiter_source, toCheck2.noiter_source);
  BOOST_REQUIRE_EQUAL(toCheck1.epssq_source, toCheck2.epssq_source);
  BOOST_REQUIRE_EQUAL(toCheck1.plaquettevalue_source, toCheck2.plaquettevalue_source);
  BOOST_REQUIRE_EQUAL(toCheck1.beta_source, toCheck2.beta_source);
  BOOST_REQUIRE_EQUAL(toCheck1.kappa_source, toCheck2.kappa_source);
  BOOST_REQUIRE_EQUAL(toCheck1.mu_source, toCheck2.mu_source);
  BOOST_REQUIRE_EQUAL(toCheck1.c2_rec_source, toCheck2.c2_rec_source);
  BOOST_REQUIRE_EQUAL(toCheck1.mubar_source, toCheck2.mubar_source);
  BOOST_REQUIRE_EQUAL(toCheck1.epsilonbar_source, toCheck2.epsilonbar_source);
	
	BOOST_REQUIRE_EQUAL(toCheck1.field_source, toCheck2.field_source);
	BOOST_CHECK_EQUAL(toCheck1.date_source, toCheck2.date_source);
	BOOST_REQUIRE_EQUAL(toCheck1.hmcversion_source, toCheck2.hmcversion_source);
	BOOST_CHECK_EQUAL(toCheck1.hmcversion_solver_source, toCheck2.hmcversion_solver_source);
	BOOST_CHECK_EQUAL(toCheck1.date_solver_source, toCheck2.date_solver_source);
	BOOST_CHECK_EQUAL(toCheck1.solvertype_source, toCheck2.solvertype_source);
}

void writeEmptyGaugefieldFromSourcefileParameters(sourcefileparameters srcFileParams, std::string configurationName)
{
	Checksum checksum;
	
	int ns = srcFileParams.lx_source;
	int nt = srcFileParams.lt_source;
	size_t precision = srcFileParams.prec_source;
	
	const n_uint64_t bufferSize_gaugefield = getElementsOfGaugefield(srcFileParams.lx_source, srcFileParams.ly_source, srcFileParams.lz_source, srcFileParams.lt_source) * sizeof(double);
	
	char * binaryData = new char[ bufferSize_gaugefield ];
	
	//TODO: hmc version currently can not be anything else than #.# !!
	write_gaugefield (
		binaryData, bufferSize_gaugefield, checksum,
		ns, ns, ns, nt, precision,
		srcFileParams.trajectorynr_source, srcFileParams.plaquettevalue_source, srcFileParams.beta_source, srcFileParams.kappa_solver_source, srcFileParams.mu_source, srcFileParams.c2_rec_source, srcFileParams.epsilonbar_source, srcFileParams.mubar_source, srcFileParams.hmcversion_source.c_str() ,configurationName.c_str());
	
	delete binaryData;
}

BOOST_AUTO_TEST_CASE(writeGaugefield)
{
	sourcefileparameters srcFileParams_1 = setSourceFileParametersToSpecificValuesForGaugefield();
	
	//TODO: test with single?
	std::string configurationName = "conf.test";
	
	writeEmptyGaugefieldFromSourcefileParameters(srcFileParams_1, configurationName);
	
	sourcefileparameters srcFileParams_2;
	char * readBinaryData;
	srcFileParams_2.readsourcefile(configurationName.c_str(), srcFileParams_1.prec_source, &readBinaryData);
	delete readBinaryData;
	
	compareTwoSourcefileParameters(srcFileParams_1, srcFileParams_2);
}

