/** @file
 * ildg IO utilities
 *
 * Copyright 2014 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
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

#include "ildgIo.hpp"
#include "ildgIoParameters.hpp"
#include "ildgIo_gaugefield.hpp"

Matrixsu3 * ildgIo::readGaugefieldFromSourcefile(std::string ildgfile, const meta::Inputparameters * parameters, int & trajectoryNumberAtInit, double & plaq)
{
	Matrixsu3 * gf_host = nullptr; // will be allocated by the following line, init to 0 to ensure error in case we fuck up.

	InputparametersTmp test2( parameters );
	IldgIoParameters_gaugefield test(&test2);
	IldgIoReader_gaugefield reader(ildgfile, &test, &gf_host);

	trajectoryNumberAtInit = reader.getReadTrajectoryNumber();
	plaq = reader.getReadPlaquetteValue();

	return gf_host;
}

void ildgIo::writeGaugefieldToFile(std::string outputfile, std::vector<Matrixsu3> & host_buf, const meta::Inputparameters * parameters, int trajectoryNumber, double plaquetteValue)
{
	InputparametersTmp test2( parameters );
	IldgIoParameters_gaugefield test(&test2);
	IldgIoWriter_gaugefield writer(host_buf, &test, outputfile, trajectoryNumber, plaquetteValue);
}
