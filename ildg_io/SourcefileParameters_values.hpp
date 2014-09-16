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

#include "checksum.h"
#include "../meta/inputparameters.hpp"

//TODO: it may be advantageous to separate between gaugefield and fermion field parameters
class sourcefileparameters_values {
public:
  sourcefileparameters_values();
	sourcefileparameters_values(const meta::Inputparameters * parameters, int trajectoryNumber, double plaquette, Checksum checksumIn, std::string hmcVersion);
	
	std::string getInfo_ildgFormat_gaugefield();
	std::string getInfo_scidacChecksum();
	std::string getInfo_xlfInfo();
	
	int lx, ly, lz, lt, prec, num_entries, flavours, trajectorynr, time, time_solver, noiter;
	double plaquettevalue, beta, kappa, mu, c2_rec, mubar, epsilonbar, epssq, kappa_solver, mu_solver;
	Checksum checksum;
	std::string field, date, hmcversion, solvertype, hmcversion_solver, date_solver;
	
	int numberOfFermionFieldsRead;
	
	void printMetaDataToScreen(std::string sourceFilename);
	
	void checkAgainstInputparameters(const meta::Inputparameters * toCheck);
	
private:
	void set_defaults();
};
	
#endif
