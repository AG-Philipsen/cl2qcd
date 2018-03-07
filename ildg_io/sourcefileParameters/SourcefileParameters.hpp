/** @file
 * class containing parameters of the sourcefile
 *
 * Copyright (c) 2014,2015 Christopher Pinke
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

#ifndef _SOURCEFILEPARAMETERS_HPP_
#define _SOURCEFILEPARAMETERS_HPP_

#include "../checksum.hpp"
#include "../ildgIoParameters.hpp"

//TODO: it may be advantageous to separate between gaugefield and fermion field parameters
class Sourcefileparameters {
public:
  Sourcefileparameters();
	Sourcefileparameters(const IldgIoParameters * parameters, int trajectoryNumber, double plaquette, Checksum checksumIn, std::string hmcVersion);

	std::string getInfo_ildgFormat_gaugefield();
	std::string getInfo_scidacChecksum();
	std::string getInfo_xlfInfo();

	int lx, ly, lz, lt, prec, num_entries, flavours, trajectorynr, time, time_solver, noiter;
	double plaquettevalue, beta, kappa, mu, c2_rec, mubar, epsilonbar, epssq, kappa_solver, mu_solver;
	Checksum checksum;
	std::string field, date, hmcversion, solvertype, hmcversion_solver, date_solver;

	int numberOfFermionFieldsRead;

	void printMetaDataToScreen(std::string sourceFilename);

	void checkAgainstInputparameters(const IldgIoParameters * toCheck);
	void checkAgainstChecksum(Checksum checksum, bool ignoreChecksumErrors = false, std::string filename = "");

	size_t getSizeInBytes() noexcept;

private:
	void set_defaults();
};



#endif
