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
#include "limeUtilities.hpp"
#include "limeFileReader.hpp"
#include "limeFileWriter.hpp"
#include "SourcefileParameters.hpp"

#include "checksum.h"
extern "C" {
#include <lime.h>
}

//todo: add namespace ildg_io

/**
 * ILDG compatible reader class for gaugefield.
 *
 * Contains metadata of the parsed gaugefield as member.
 * TODO: change this.
 *
 * @param[in] file      The file to read the gauge configuration from
 * @param[in] precision The precision expected for the gaugefield.
 * @param[out] data    The loaded gaugefield
 */
class IldgIoReader_gaugefield : public LimeFileReader
{
public:
	IldgIoReader_gaugefield(std::string file, int precision, char ** data);
	//todo: make this private?
	Sourcefileparameters parameters;
private:
	void readMetaDataFromLimeFile();
	void readDataFromLimeFile(char ** destination);
	int calcNumberOfEntriesBasedOnFieldType(std::string fieldType);
	LimeFileProperties extractMetaDataFromLimeEntry(LimeReader * r, LimeHeaderData limeHeaderData);
	size_t	sizeOfGaugefieldBuffer();
	void extractBinaryDataFromLimeEntry(LimeReader * r, LimeHeaderData limeHeaderData, char ** destination);
	void readLimeFile(char ** destination);
	void extractMetadataFromLimeFile();
	void extractDataFromLimeFile(char ** destination);
	void extractInformationFromLimeEntry(LimeReader * r, char ** destination);
	void goThroughLimeRecords(LimeReader * r, char ** destination);

	int checkLimeEntryForFermionInformations(std::string lime_type);
	bool checkLimeEntryForBinaryData(std::string lime_type);

	LimeFilePropertiesCollector limeFileProp;
	LimeEntryTypes limeEntryTypes;
};

/**
 * ILDG compatible writer class for gaugefield.
 *
 * \param binary_data The gaugefield in binary format.
 * \param num_bytes The number of bytes to be written.
 * \param srcFileParameters Collection of parameters associated with the gaugefield.
 */
class IldgIoWriter_gaugefield: public LimeFileWriter
{
public:
	IldgIoWriter_gaugefield(char * binary_data, n_uint64_t num_bytes, Sourcefileparameters srcFileParameters, std::string filenameIn);
};


#endif /* _READGAUGEH_ */
