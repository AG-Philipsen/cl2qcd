/** @file
 * Testcases for lime utilities
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
#define BOOST_TEST_MODULE limeUtilities
#include <boost/test/unit_test.hpp>

#include "limeUtilities.hpp"

#include "../executables/exceptions.h"

int expectedPrecision = 64;
std::string nameOfExistingGaugefieldFile = std::string(SOURCEDIR) + "/ildg_io/conf.00200";

BOOST_AUTO_TEST_CASE(readInLimeFile_failureWithFileException)
{
  std::string nameOfNonexistingLimeFile = "thisfileshouldnotbethere";
  char * bufferToStore;
	
	BOOST_CHECK_THROW(LimeFileReader srcFileParams(nameOfNonexistingLimeFile, expectedPrecision, &bufferToStore), File_Exception);
}

