/*
 * Copyright 2012, 2013, 2014 Lars Zeidlewicz, Christopher Pinke,
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

#include "ildgIo_gaugefield.hpp"

#include "limeUtilities.hpp"
#include "../host_functionality/logger.hpp"

using namespace ildgIo;

void checkLimeFileForFieldType(std::string fieldTypeIn) throw(std::logic_error)
{
	if (fieldTypeIn != "su3gauge")
	{
		throw std::logic_error("LIME file does not seem to include gaugefield data. Aborting...");
	}
}

//todo: make char ** std::vector<char*>
IldgIoReader_gaugefield::IldgIoReader_gaugefield(std::string sourceFilenameIn, int desiredPrecisionIn, char ** destination) : LimeFileReader(sourceFilenameIn, desiredPrecisionIn)
{
	if ( limeFileProp.numberOfBinaryDataEntries >= 1 )
	{
		checkLimeFileForFieldType(parameters.field);
		extractDataFromLimeFile(destination);
	}
	else 
	{
		throw std::logic_error("LIME file does not seem to include binary data. Aborting...");
	}
}

void* createVoidPointerFromString(std::string stringIn) noexcept
{
	return (void*) stringIn.c_str();
}

IldgIoWriter_gaugefield::IldgIoWriter_gaugefield(char * binary_data, n_uint64_t num_bytes, Sourcefileparameters srcFileParameters, std::string filenameIn): LimeFileWriter(filenameIn)
{
	logger.info() << "writing gaugefield to lime-file \""  + filenameIn + "\"";
	
	std::string xlfInfo = srcFileParameters.getInfo_xlfInfo();
	std::string scidac_checksum = srcFileParameters.getInfo_scidacChecksum();
	std::string ildgFormat = srcFileParameters.getInfo_ildgFormat_gaugefield();
	
	writeMemoryToLimeFile( createVoidPointerFromString(xlfInfo), xlfInfo.size(), limeEntryTypes["xlf"]);
	writeMemoryToLimeFile( createVoidPointerFromString(ildgFormat), ildgFormat.size(), limeEntryTypes["ildg"]);
	writeMemoryToLimeFile( binary_data, num_bytes, limeEntryTypes["ildg binary data"]);
	writeMemoryToLimeFile( createVoidPointerFromString(scidac_checksum), scidac_checksum.size(), limeEntryTypes["scidac checksum"]);
}


