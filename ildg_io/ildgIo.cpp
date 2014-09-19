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

#include "ildgIo_gaugefield.hpp"
#include "../meta/util.hpp"
#include "checksum.h"
#include <cassert>
#include "../executables/exceptions.h"

Matrixsu3 * ildgIo::readGaugefieldFromSourcefile(std::string ildgfile, const meta::Inputparameters * parameters, int & trajectoryNumberAtInit, double & plaq)
{
	Matrixsu3 * gf_host;
	//todo: this should not be that explicit here!	
	gf_host = new Matrixsu3[meta::get_vol4d(*parameters) * 4];

	IldgIoReader_gaugefield reader(ildgfile, parameters->get_precision(), parameters, gf_host);

	//todo: make access member fcts.
	trajectoryNumberAtInit = reader.parameters.trajectorynr;
	plaq = reader.parameters.plaquettevalue;

	return gf_host;
}

static size_t getBufferSize_gaugefield(const meta::Inputparameters * parameters) noexcept
{
	return 2 * NC * NC * NDIM * meta::get_volspace(*parameters) * parameters->get_ntime() * sizeof(hmc_float);
}

void ildgIo::writeGaugefieldToFile(std::string outputfile, Matrixsu3 * host_buf, const meta::Inputparameters * parameters, int trajectoryNumber, double plaquetteValue)
{
	//todo: move this into the writer class!
	const size_t gaugefield_buf_size = getBufferSize_gaugefield(parameters);
	IldgIoWriter_gaugefield writer(host_buf, gaugefield_buf_size, parameters, outputfile, trajectoryNumber, plaquetteValue);
}
