/** @file
 * Reading of gauge field from files.
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#ifndef _READGAUGEH_
#define _READGAUGEH_

#include "../common_header_files/types.h"
#include "../executables/exceptions.h"

//todo: remove this eventually
#include "limeFileReader.hpp"
#include "limeFileWriter.hpp"
#include "SourcefileParameters.hpp"

namespace ildgIo {

	/**
	* ILDG compatible reader class for gaugefield.
	*
	* Contains metadata of the parsed gaugefield as member.
	* TODO: change this.
	*
	* @param[in] file      The file to read the gauge configuration from
	* @param[in] precision The precision expected for the gaugefield.
	* @param[in] parameters The input parameters @Todo: remove this
	* @param[out] data    The loaded gaugefield
	*/
	class IldgIoReader_gaugefield : public LimeFileReader
	{
	public:
		IldgIoReader_gaugefield(std::string file, int precision, const meta::Inputparameters * parameters, Matrixsu3 * data);
	};

	/**
	* ILDG compatible writer class for gaugefield.
	*
	* \param data The gaugefield to write
	* \param num_bytes The number of bytes to be written.
	* \param srcFileParameters Collection of parameters associated with the gaugefield.
	*/
	class IldgIoWriter_gaugefield: public LimeFileWriter
	{
	public:
		IldgIoWriter_gaugefield(Matrixsu3 * data, n_uint64_t num_bytes, const meta::Inputparameters * parameters, std::string filenameIn, int trajectoryNumber, double plaquetteValue);
	};
	
	Checksum calculate_ildg_checksum(const char * buf, size_t nbytes, const meta::Inputparameters& inputparameters);
	void copy_gaugefield_from_ildg_format(Matrixsu3 * gaugefield, char * gaugefield_tmp, int check, const 	meta::Inputparameters& parameters);
	void copy_gaugefield_to_ildg_format(char * dest, Matrixsu3 * source_in, const meta::Inputparameters& parameters);
}

#endif /* _READGAUGEH_ */
