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

#include "limeFileWriter.hpp"

#include "../host_functionality/logger.hpp"
#include <boost/lexical_cast.hpp>
#include "../executables/exceptions.h"
//http://www.ridgesolutions.ie/index.php/2013/05/30/boost-link-error-undefined-reference-to-boostfilesystemdetailcopy_file/
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>

LimeFileWriter::LimeFileWriter(std::string filenameIn) : LimeFile_basic(filenameIn)
{
	if ( boost::filesystem::exists( filenameIn ) )
	{
		const std::string backupFilename = filenameIn + "_backup";
		logger.warn() << "Found existing file of name \"" << filenameIn << "\". Store this to \"" << backupFilename << "\"...";
		boost::filesystem::copy_file( filenameIn, backupFilename, boost::filesystem::copy_option::overwrite_if_exists);
	}

	MB_flag = 1;
	ME_flag = 1;
	writtenBytes = 0;
	writer = NULL;
	
	outputfile = fopen(filename.c_str(), "w");
	writer = limeCreateWriter(outputfile);
}

void LimeFileWriter::closeLimeFile()
{
	fclose(outputfile);
	limeDestroyWriter(writer);
	writer = NULL;
	logger.info() << "  " << (float) ( (float) (writtenBytes) / 1024 / 1024 ) << " MBytes were written to the lime file " << filename;
}

LimeFileWriter::~LimeFileWriter()
{
	if (writer != NULL)
	{
		closeLimeFile();
	}
}

void writeLimeHeaderToLimeFile(LimeRecordHeader * header, LimeWriter * writer)
{
	int returnCode = 0;
	
	returnCode = limeWriteRecordHeader(header, writer);
	if ( returnCode != LIME_SUCCESS )
	{
		throw Print_Error_Message( "Could not write header to LIME file. Return code: " + boost::lexical_cast<std::string>(returnCode), __FILE__, __LINE__);
	}
}

void LimeFileWriter::writeMemoryToLimeFile(void * memoryPointer, n_uint64_t bytes, std::string description)
{
	logger.debug() << "writing \"" + description + "\" to lime file...";

	n_uint64_t bytesToBeWritten = bytes;
	int returnCode = 0;
	
	LimeRecordHeader * header = limeCreateHeader(this->MB_flag, this->ME_flag, (char*) description.c_str(), bytesToBeWritten);
	this->ME_flag++;
	writeLimeHeaderToLimeFile(header, this->writer);
	limeDestroyHeader(header);
	
	returnCode = limeWriteRecordData( memoryPointer, &bytesToBeWritten, this->writer);
	if ( returnCode != LIME_SUCCESS )
	{
		throw Print_Error_Message( "Could not write to LIME file. Return code: " + boost::lexical_cast<std::string>(returnCode), __FILE__, __LINE__);
	}
	else if ( bytesToBeWritten != bytes )
	{
		throw Print_Error_Message( "There was an error writing to Lime file...", __FILE__, __LINE__);
	}

	this->writtenBytes += bytesToBeWritten;
}
