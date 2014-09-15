/** @file
 * class containing parameters of the sourcefile
 *
 * Copyright 2014, Christopher Pinke
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

#ifndef _SOURCEFILEPARAMETERS_HPP_
#define _SOURCEFILEPARAMETERS_HPP_

#include "../host_functionality/logger.hpp"
#include "checksum.h"
#include "../meta/inputparameters.hpp"

//TODO: separate between gaugefield and fermion field parameters
//TODO: remove "source" postfix
class sourcefileparameters_values {
public:
  sourcefileparameters_values();
	sourcefileparameters_values(const meta::Inputparameters * parameters, int trajectoryNumber, double plaquette);
	
	int lx_source, ly_source, lz_source, lt_source, prec_source, num_entries_source, flavours_source,
	    trajectorynr_source, time_source, time_solver_source, noiter_source;
	double plaquettevalue_source, beta_source, kappa_source, mu_source, c2_rec_source, mubar_source, epsilonbar_source, epssq_source, kappa_solver_source, mu_solver_source;
	Checksum checksum;
	
	std::string field_source;
	std::string date_source;
	std::string hmcversion_source;
	std::string solvertype_source;
	std::string hmcversion_solver_source;
	std::string date_solver_source;
	
	int numberOfFermionFieldsRead;
	
	void printMetaDataToScreen(std::string sourceFilename);
	
private:
	void set_defaults();
};
	
#endif