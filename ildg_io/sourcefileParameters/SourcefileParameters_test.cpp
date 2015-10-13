/** @file
 * Testcases for the sourcefileparamters_values class
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
#define BOOST_TEST_MODULE Sourcefileparameters
#include <boost/test/unit_test.hpp>

#include "SourcefileParameters.hpp"

#include "../executables/exceptions.h"

// TODO: Remove all occurences of meta::Inputparameters and Inputparameters. The tests should only use IlgdIoParameters_gaugefield!!

void checkDefaults(Sourcefileparameters toCheck)
{
  BOOST_REQUIRE_EQUAL(toCheck.lx, 0);
  BOOST_REQUIRE_EQUAL(toCheck.ly, 0);
  BOOST_REQUIRE_EQUAL(toCheck.lz, 0);
  BOOST_REQUIRE_EQUAL(toCheck.lt, 0);
  BOOST_REQUIRE_EQUAL(toCheck.prec, 0);
  BOOST_REQUIRE_EQUAL(toCheck.num_entries, 0);
  BOOST_REQUIRE_EQUAL(toCheck.flavours, 0);
  BOOST_REQUIRE_EQUAL(toCheck.trajectorynr, 0);
  BOOST_REQUIRE_EQUAL(toCheck.time, 0);
  BOOST_REQUIRE_EQUAL(toCheck.time_solver, 0);
  BOOST_REQUIRE_EQUAL(toCheck.noiter, 0);
  BOOST_REQUIRE_EQUAL(toCheck.plaquettevalue, 0);
  BOOST_REQUIRE_EQUAL(toCheck.beta, 0);
  BOOST_REQUIRE_EQUAL(toCheck.kappa, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mu, 0);
  BOOST_REQUIRE_EQUAL(toCheck.c2_rec, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mubar, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epsilonbar, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epssq, 0);
  BOOST_REQUIRE_EQUAL(toCheck.kappa_solver, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mu_solver, 0);
	
  BOOST_REQUIRE_EQUAL(toCheck.numberOfFermionFieldsRead, 0);
	
	BOOST_REQUIRE_EQUAL(toCheck.field, "");
	BOOST_REQUIRE_EQUAL(toCheck.date, "");
	BOOST_REQUIRE_EQUAL(toCheck.hmcversion, "");
	BOOST_REQUIRE_EQUAL(toCheck.solvertype, "");
	BOOST_REQUIRE_EQUAL(toCheck.hmcversion_solver, "");
	BOOST_REQUIRE_EQUAL(toCheck.date_solver, "");
	
	Checksum checksum;
	BOOST_REQUIRE(toCheck.checksum == checksum);
}

BOOST_AUTO_TEST_CASE(defaults)
{
  Sourcefileparameters srcFileParams;
  checkDefaults(srcFileParams);
}

void checkSpecificParameters(Sourcefileparameters toCheck)
{
  BOOST_REQUIRE_EQUAL(toCheck.lx, 65);
  BOOST_REQUIRE_EQUAL(toCheck.ly, 65);
  BOOST_REQUIRE_EQUAL(toCheck.lz, 65);
  BOOST_REQUIRE_EQUAL(toCheck.lt, 41);
  BOOST_REQUIRE_EQUAL(toCheck.prec, 32);
  BOOST_REQUIRE_EQUAL(toCheck.num_entries, 0);
  BOOST_REQUIRE_EQUAL(toCheck.flavours, 0);
  BOOST_REQUIRE_EQUAL(toCheck.trajectorynr, -1982);
  BOOST_REQUIRE_EQUAL(toCheck.time, 0);
  BOOST_REQUIRE_EQUAL(toCheck.time_solver, 0);
  BOOST_REQUIRE_EQUAL(toCheck.noiter, 0);
  BOOST_REQUIRE_EQUAL(toCheck.plaquettevalue, -4.321);
  BOOST_REQUIRE_EQUAL(toCheck.beta, 4.5);
  BOOST_REQUIRE_EQUAL(toCheck.kappa, -12.345);
  BOOST_REQUIRE_EQUAL(toCheck.mu, 23.41);
  BOOST_REQUIRE_EQUAL(toCheck.c2_rec, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mubar, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epsilonbar, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epssq, 0);
  BOOST_REQUIRE_EQUAL(toCheck.kappa_solver, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mu_solver, 0);
	
  BOOST_REQUIRE_EQUAL(toCheck.numberOfFermionFieldsRead, 0);
	
	BOOST_REQUIRE_EQUAL(toCheck.field, "");
	BOOST_REQUIRE_EQUAL(toCheck.date, "");
	BOOST_REQUIRE_EQUAL(toCheck.hmcversion, "3.95");
	BOOST_REQUIRE_EQUAL(toCheck.solvertype, "");
	BOOST_REQUIRE_EQUAL(toCheck.hmcversion_solver, "");
	BOOST_REQUIRE_EQUAL(toCheck.date_solver, "");
	
	Checksum checksum(1,2);
	BOOST_REQUIRE(toCheck.checksum == checksum);	
}

BOOST_AUTO_TEST_CASE(initFromParameters)
{
	int trajectoryNumber = -1982;
	double plaquette = -4.321;
	std::string hmcVersion = "3.95";
	Checksum checksum(1,2);
	const char * _params[] = {"foo", 
		"--ntime=41", "--nspace=65", "--kappa=-12.345", "--prec=32", "--beta=4.5", "--mu=23.41"
	};
	meta::Inputparameters parameters(7, _params);
	const physics::lattices::GaugefieldParametersImplementation tmp (&parameters);
	Inputparameters test2( &tmp );
	IldgIoParameters_gaugefield test(&test2);
	
	Sourcefileparameters srcFileParams(&test, trajectoryNumber, plaquette, checksum, hmcVersion);
	checkSpecificParameters(srcFileParams);
}

BOOST_AUTO_TEST_CASE(checkAgainstParameters_exception)
{
	int trajectoryNumber = -1982;
	double plaquette = -4.321;
	std::string hmcVersion = "3.95";
	Checksum checksum(1,2);
	const char * _params[] = {"foo", 
		"--ntime=41", "--nspace=65", "--prec=32", "--beta=4.5",  "--kappa=-12.345", "--mu=23.41"
	};
	int numberOfParameters = 7;
	meta::Inputparameters parameters(numberOfParameters, _params);
	const physics::lattices::GaugefieldParametersImplementation tmp (&parameters);
	Inputparameters test2( &tmp );
	IldgIoParameters_gaugefield test(&test2);
	Sourcefileparameters srcFileParams(&test, trajectoryNumber, plaquette, checksum, hmcVersion);

	for(int iteration = 1; iteration < 3 + 1; iteration++)
	{
		const meta::Inputparameters standardParameters(iteration, _params);
		const physics::lattices::GaugefieldParametersImplementation tmp (&standardParameters);
		Inputparameters test2( &tmp );
		IldgIoParameters_gaugefield test(&test2);
		BOOST_REQUIRE_THROW(srcFileParams.checkAgainstInputparameters(&test), std::invalid_argument );
	}
}

BOOST_AUTO_TEST_CASE(checkAgainstChecksum)
{
	Checksum checksum;
	Sourcefileparameters values;
	
	BOOST_REQUIRE(checksum == values.checksum);
}

BOOST_AUTO_TEST_CASE(checkAgainstChecksum_exception)
{
	Checksum checksum(1,2);
	Sourcefileparameters values;
	
	BOOST_REQUIRE_THROW(values.checkAgainstChecksum(checksum, false, "nameOfFile"), File_Exception );
}

BOOST_AUTO_TEST_CASE(checkAgainstChecksum_noExceptionByParameters)
{
	Checksum checksum(1,2);
	Sourcefileparameters values;
	
	const char * _params[] = {"foo", "--ignore_checksum_errors=true" };

	meta::Inputparameters parameters(2, _params);
	const physics::lattices::GaugefieldParametersImplementation tmp (&parameters);
	Inputparameters test2( &tmp );
	IldgIoParameters_gaugefield test(&test2);
	Sourcefileparameters srcFileParams(&test, 123, 4.56, checksum, "8.9");
	
	BOOST_CHECK_NO_THROW(values.checkAgainstChecksum(checksum, parameters.get_ignore_checksum_errors(), "nameOfFile"));
}

BOOST_AUTO_TEST_CASE(checkAgainstChecksum_exceptionByParameters)
{
	Checksum checksum(1,2);
	Sourcefileparameters values;
	
	const char * _params[] = {"foo", "--ignore_checksum_errors=false" };
	meta::Inputparameters parameters(2, _params);
	const physics::lattices::GaugefieldParametersImplementation tmp (&parameters);
	Inputparameters test2( &tmp );
	IldgIoParameters_gaugefield test(&test2);
	Sourcefileparameters srcFileParams(&test, 123, 4.56, checksum, "8.9");
	
	BOOST_REQUIRE_THROW(values.checkAgainstChecksum(checksum, parameters.get_ignore_checksum_errors(), "nameOfFile"),  File_Exception);
}

BOOST_AUTO_TEST_CASE(checkAgainstChecksum_exceptionByParameters_defaultSetting)
{
	Checksum checksum(1,2);
	Sourcefileparameters values;
	
	const char * _params[] = {"foo"};
	meta::Inputparameters parameters(1, _params);
	const physics::lattices::GaugefieldParametersImplementation tmp (&parameters);
	Inputparameters test2( &tmp );
	IldgIoParameters_gaugefield test(&test2);
	Sourcefileparameters srcFileParams(&test, 123, 4.56, checksum, "8.9");
	
	BOOST_REQUIRE_THROW(values.checkAgainstChecksum(checksum, parameters.get_ignore_checksum_errors(), "nameOfFile"),  File_Exception);
}
