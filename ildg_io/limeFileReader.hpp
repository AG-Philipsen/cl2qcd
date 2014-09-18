/** @file
 * Reader for LIME files
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

#ifndef _LIMEFILEREADER_HPP_
#define _LIMEFILEREADER_HPP_

#include "limeUtilities.hpp"
#include "SourcefileParameters.hpp"

class LimeFileReader
{
public:
	//todo: remove precision?
	LimeFileReader(std::string sourceFilenameIn, int precision, char ** data);
	~LimeFileReader();
	//todo: make this private or remove/find better solution?
	Sourcefileparameters parameters;
protected:
	void openFile();
	void closeFile();
	void readMetaDataFromLimeFile();
	void extractMetadataFromLimeFile();
	void readLimeFile(char ** destination);
	void extractDataFromLimeFile(char ** destination);
	void readDataFromLimeFile(char ** destination);
	void goThroughLimeRecords(char ** destination);
	void extractInformationFromLimeEntry(char ** destination);
	LimeFileProperties extractMetaDataFromLimeEntry(LimeHeaderData limeHeaderData);
	void extractBinaryDataFromLimeEntry(LimeHeaderData limeHeaderData, char ** destination);
	int calcNumberOfEntriesBasedOnFieldType(std::string fieldType);
	
	void handleLimeEntry_xlf(Sourcefileparameters & parameters, char * buffer, std::string lime_type);
	void handleLimeEntry_ildg(Sourcefileparameters & parameters, char * buffer, std::string lime_type, size_t numberOfBytes);
	
	size_t	sizeOfGaugefieldBuffer();
	
	int checkLimeEntryForFermionInformations(std::string lime_type);
	bool checkLimeEntryForBinaryData(std::string lime_type);
	
	std::string sourceFilename;
	int desiredPrecision;
	LimeReader * limeReader;
	FILE *limeFileOpenedForReading;
	
	LimeFilePropertiesCollector limeFileProp;
	LimeEntryTypes limeEntryTypes;
};

#endif