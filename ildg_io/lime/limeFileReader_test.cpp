/** @file
 * Testcases for lime utilities
 *
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE limeFileReader
#include <boost/test/unit_test.hpp>

#include "limeFileReader.hpp"

#include "../executables/exceptions.hpp"

int expectedPrecision = 64;
std::string nameOfExistingGaugefieldFile = std::string(SOURCEDIR) + "/ildg_io/conf.00200";

BOOST_AUTO_TEST_CASE(readInLimeFile_failureWithFileException)
{
  std::string nameOfNonexistingLimeFile = "thisfileshouldnotbethere";

	BOOST_CHECK_THROW(LimeFileReader srcFileParams(nameOfNonexistingLimeFile, expectedPrecision), File_Exception);
}

BOOST_AUTO_TEST_CASE(readInLimeFile_failureWithWrongPrecision)
{
  int wrongPrecision = 27;
  BOOST_CHECK_THROW(LimeFileReader srcFileParams(nameOfExistingGaugefieldFile, wrongPrecision), std::exception);
}

BOOST_AUTO_TEST_CASE(readInLimeFile_Success)
{
  BOOST_REQUIRE_NO_THROW(LimeFileReader srcFileParams(nameOfExistingGaugefieldFile, expectedPrecision));
}

void checkMetadataOfSpecificGaugefieldFile(Sourcefileparameters toCheck)
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

BOOST_AUTO_TEST_CASE(readInLimeFile_CheckMetadata)
{
  LimeFileReader srcFileParams(nameOfExistingGaugefieldFile, expectedPrecision);
  checkMetadataOfSpecificGaugefieldFile(srcFileParams.parameters);
}

BOOST_AUTO_TEST_CASE(readInLimeFile_CheckChecksum)
{
  uint32_t referenceChecksumA = 171641288;
  uint32_t referenceChecksumB = 3618036129;
  LimeFileReader srcFileParams(nameOfExistingGaugefieldFile, expectedPrecision);
  Checksum referenceChecksum(referenceChecksumA, referenceChecksumB);
  BOOST_REQUIRE_EQUAL(referenceChecksum == srcFileParams.parameters.checksum, true);
}

std::string nameOfExistingGaugefieldFileFromTmlqcd = std::string(SOURCEDIR) + "/ildg_io/conf.tmlqcd";

BOOST_AUTO_TEST_CASE(readInLimeFile_FromTmlqcd_CheckChecksum)
{
  uint32_t referenceChecksumA = 398012545;
  uint32_t referenceChecksumB = 1610757546;
  LimeFileReader srcFileParams(nameOfExistingGaugefieldFileFromTmlqcd, expectedPrecision);
  Checksum referenceChecksum(referenceChecksumA, referenceChecksumB);
  BOOST_REQUIRE_EQUAL(referenceChecksum == srcFileParams.parameters.checksum, true);
}
