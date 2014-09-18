/** @file
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

#include "limeFileReader.hpp"

#include "../host_functionality/logger.hpp"
#include <boost/lexical_cast.hpp>
#include "../executables/exceptions.h"

void checkIfFileExists(std::string file) throw(File_Exception)
{
	FILE * checker;
	checker = fopen(file.c_str(), "r");
	if(checker == 0) {
		throw File_Exception(file);
	}
	fclose(checker);
	return;
}

LimeFileReader::LimeFileReader(std::string sourceFilenameIn, int precision, char ** data)
{
	sourceFilename = sourceFilenameIn;
	desiredPrecision = precision;
	
	checkIfFileExists(sourceFilename);
	
// 	extractMetadataFromLimeFile();
}

LimeFileReader::~LimeFileReader()
{
}

void LimeFileReader::openFile()
{
	limeFileOpenedForReading = fopen (sourceFilename.c_str(), "r");
	limeReader = limeCreateReader(limeFileOpenedForReading);
}

void LimeFileReader::closeFile()
{
	limeDestroyReader(limeReader);
	fclose(limeFileOpenedForReading);
}
