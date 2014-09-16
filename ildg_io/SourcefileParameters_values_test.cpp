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
#define BOOST_TEST_MODULE ildg_read_gaugefield
#include <boost/test/unit_test.hpp>

#include "SourcefileParameters_values.hpp"

void checkDefaults(sourcefileparameters_values toCheck)
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
	
  BOOST_REQUIRE_EQUAL(toCheck.numberOfFermionFieldsRead, 0);
	
	BOOST_REQUIRE_EQUAL(toCheck.field_source, "");
	BOOST_REQUIRE_EQUAL(toCheck.date_source, "");
	BOOST_REQUIRE_EQUAL(toCheck.hmcversion_source, "");
	BOOST_REQUIRE_EQUAL(toCheck.solvertype_source, "");
	BOOST_REQUIRE_EQUAL(toCheck.hmcversion_solver_source, "");
	BOOST_REQUIRE_EQUAL(toCheck.date_solver_source, "");
}

BOOST_AUTO_TEST_CASE(defaults)
{
  sourcefileparameters_values srcFileParams;
  checkDefaults(srcFileParams);
}

void checkSpecificParameters(sourcefileparameters_values toCheck)
{
  BOOST_REQUIRE_EQUAL(toCheck.lx_source, 65);
  BOOST_REQUIRE_EQUAL(toCheck.ly_source, 65);
  BOOST_REQUIRE_EQUAL(toCheck.lz_source, 65);
  BOOST_REQUIRE_EQUAL(toCheck.lt_source, 41);
  BOOST_REQUIRE_EQUAL(toCheck.prec_source, 32);
  BOOST_REQUIRE_EQUAL(toCheck.num_entries_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.flavours_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.trajectorynr_source, -1982);
  BOOST_REQUIRE_EQUAL(toCheck.time_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.time_solver_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.noiter_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.plaquettevalue_source, -4.321);
  BOOST_REQUIRE_EQUAL(toCheck.beta_source, 4.5);
  BOOST_REQUIRE_EQUAL(toCheck.kappa_source, -12.345);
  BOOST_REQUIRE_EQUAL(toCheck.mu_source, 23.41);
  BOOST_REQUIRE_EQUAL(toCheck.c2_rec_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mubar_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epsilonbar_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.epssq_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.kappa_solver_source, 0);
  BOOST_REQUIRE_EQUAL(toCheck.mu_solver_source, 0);
	
  BOOST_REQUIRE_EQUAL(toCheck.numberOfFermionFieldsRead, 0);
	
	BOOST_REQUIRE_EQUAL(toCheck.field_source, "");
	BOOST_REQUIRE_EQUAL(toCheck.date_source, "");
	BOOST_REQUIRE_EQUAL(toCheck.hmcversion_source, "3.95");
	BOOST_REQUIRE_EQUAL(toCheck.solvertype_source, "");
	BOOST_REQUIRE_EQUAL(toCheck.hmcversion_solver_source, "");
	BOOST_REQUIRE_EQUAL(toCheck.date_solver_source, "");
}

BOOST_AUTO_TEST_CASE(initFromParameters)
{
	int trajectoryNumber = -1982;
	double plaquette = -4.321;
	std::string hmcVersion = "3.95";
	const char * _params[] = {"foo", 
		"--ntime=41", "--nspace=65", "--kappa=-12.345", "--prec=32", "--beta=4.5", "--mu=23.41"
	};
	meta::Inputparameters parameters(7, _params);
	
  sourcefileparameters_values srcFileParams(&parameters, trajectoryNumber, plaquette, hmcVersion);
  checkSpecificParameters(srcFileParams);
}

