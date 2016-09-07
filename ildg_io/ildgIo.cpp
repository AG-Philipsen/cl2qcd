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

Matrixsu3 * ildgIo::readGaugefieldFromSourcefile(std::string ildgfile, const physics::lattices::GaugefieldParametersInterface * parameters, int & trajectoryNumberAtInit)
{
	Matrixsu3 * gf_host = nullptr; // will be allocated by the following line, init to 0 to ensure error in case we fuck up.

	//NOTE: this is a workaround, because this call:
	//IldgIoParameters_gaugefield ildgIoParameters = createIldgIoParameters( parameters );
	// does not work as the Inputparameters instance created in this fct. goes out of scope, causing a segfault.

	Inputparameters parameters2( parameters );
	IldgIoParameters_gaugefield ildgIoParameters(&parameters2);

	IldgIoReader_gaugefield reader(ildgfile, &ildgIoParameters, &gf_host);

	trajectoryNumberAtInit = reader.getReadTrajectoryNumber();


	//TODO: insert here the plaquette check, calculating it on the host using gf_host and comparing to plaq
	//hmc_float plaqSourcefile = reader.getReadPlaquetteValue();
	//hmc_float plaq = calculatePlaquette(gf_host); <--- this function does not exist yet!
	//check_plaq(plaq, plaqSourcefile);


	return gf_host;
}

void ildgIo::writeGaugefieldToFile(std::string outputfile, std::vector<Matrixsu3> & host_buf, const physics::lattices::GaugefieldParametersInterface * parameters, int trajectoryNumber)
{
	//NOTE: this is a workaround, because this call:
	//IldgIoParameters_gaugefield ildgIoParameters = createIldgIoParameters( parameters );
	// does not work as the Inputparameters instance created in this fct. goes out of scope, causing a segfault.

	Inputparameters parameters2( parameters );
	IldgIoParameters_gaugefield ildgIoParameters(&parameters2);


	//TODO: calculate here the plaquette on the host, and then give its value to the writer to put it as metainformation in the lime file
	//hmc_float plaquetteValue = calculatePlaquette(gf_host); <--- this function does not exist yet!
	//For the moment we write a fake number like 0, it is only temporary!
	hmc_float plaquetteValue = 0.0;

	IldgIoWriter_gaugefield writer(host_buf, &ildgIoParameters, outputfile, trajectoryNumber, plaquetteValue);
}


static void check_plaq(const hmc_float plaquette, double plaqSourcefile)
{
    logger.info() << "Checking plaquette against sourcefile value...";
    std::string msg = "Minor parameters do not match: ";
    hmc_float float1, float2;
    std::string testobj = msg + "plaquette";
    float1 = plaquette;
    float2 = plaqSourcefile;
    if(float1 != float2) {
        logger.warn() << testobj;
        logger.warn() << "\tExpected: " << float1 << "\tFound: " << float2;
    }

    logger.info() << "...done";
    return;
}
