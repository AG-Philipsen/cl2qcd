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
#include "../SourcefileParameters.hpp"

class LimeFileReader : public LimeFileReader_basic
{
public:
	//todo: remove precision?
	LimeFileReader(std::string sourceFilenameIn, int precision);
	//todo: make this private or remove/find better solution?
	Sourcefileparameters parameters;
	
	int getReadTrajectoryNumber() noexcept;
	double getReadPlaquetteValue() noexcept;
protected:
	void readMetaDataFromLimeFile();
	void extractMetadataFromLimeFile();
	void readLimeFile(char ** destination, size_t expectedNumberOfBytes = 0);
	void extractDataFromLimeFile(char ** destination, size_t expectedNumberOfBytes);
	void goThroughLimeRecords(char ** destination, size_t expectedNumberOfBytes);
	void extractInformationFromLimeEntry(char ** destination, size_t expectedNumberOfBytes);
	LimeFileProperties extractMetaDataFromLimeEntry(LimeHeaderData limeHeaderData);
	void extractBinaryDataFromLimeEntry(LimeHeaderData limeHeaderData, char ** destination, size_t expectedNumberOfBytes);
	
	void handleLimeEntry_xlf(Sourcefileparameters & parameters, char * buffer, std::string lime_type);
	void handleLimeEntry_ildg(Sourcefileparameters & parameters, char * buffer, std::string lime_type, size_t numberOfBytes);
	void handleLimeEntry_scidacChecksum(char * buffer, std::string lime_type, size_t numberOfBytes);
	void handleLimeEntry_inverter(std::string lime_type) throw(std::logic_error);
	void handleLimeEntry_etmcPropagator(std::string lime_type) throw(std::logic_error);
	
	int desiredPrecision;
};

#endif