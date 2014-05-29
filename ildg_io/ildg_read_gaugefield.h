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
#include "SourcefileParameters_values.hpp"

#include "checksum.h"
extern "C" {
#include <lime.h>
}

//todo: add namespace ildg_io

/**
 * Parser class for a stored gaugefield.
 *
 * Contains metadata of the parsed gaugefield as members.
 */
class sourcefileparameters : public sourcefileparameters_values {
public:
	/**
	 * Read gauge configuration from the given file into the given array.
	 *
	 * @param[in] file      The file to read the gauge configuration from
	 * @param[in] precision The precision expected for the gaugefield.
	 * @param[out] array    The loaded gaugefield
	 */
  void readsourcefile(std::string file, int precision, char ** data);
	
 private:
	void readMetaDataFromLimeFile();
	void readDataFromLimeFile(char ** destination);
	int calcNumberOfEntriesBasedOnFieldType(std::string fieldType);
	int calcNumberOfEntriesForDiracFermionfield();
	int calcNumberOfEntriesForGaugefield();
	void checkLimeEntryForInverterInfos(std::string lime_type, LimeReader *r, size_t nbytes);
	void checkLimeEntryForXlfInfos(std::string lime_type, LimeReader *r, size_t nbytes);
	void checkLimeEntryForXlmInfos(std::string lime_type, LimeReader *r, size_t nbytes);
	void checkLimeEntryForScidacChecksum(std::string lime_type, LimeReader *r, size_t nbytes);
	LimeFileProperties extractMetaDataFromLimeEntry(LimeReader * r, LimeHeaderData limeHeaderData);
	size_t	sizeOfGaugefieldBuffer();
	void extractBinaryDataFromLimeEntry_NeedsDifferentName(LimeReader * r, LimeHeaderData limeHeaderData, char ** destination);
	void extractBinaryDataFromLimeEntry(LimeReader * r, char ** destination, LimeHeaderData limeHeaderData);
	void readLimeFile(char ** destination);
	void extractMetadataFromLimeFile();
	void extractDataFromLimeFile(char ** destination);
	void extractInformationFromLimeEntry(LimeReader * r, char ** destination);
	void goThroughLimeRecords(LimeReader * r, char ** destination);

	LimeFilePropertiesCollector limeFileProp;
	
	std::string sourceFilename;
	int desiredPrecision;
};

#endif /* _READGAUGEH_ */
