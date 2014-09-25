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

#include "../../host_functionality/logger.hpp"
#include <boost/lexical_cast.hpp>
#include "../../executables/exceptions.h"
#include "../sourcefileParameters/SourcefileParameters_utilities.hpp"

int checkLimeEntryForFermionInformations(std::string lime_type, LimeEntryTypes limeEntryTypes)
{
	return ( 
		limeEntryTypes["propagator"] == lime_type || 
		limeEntryTypes["inverter"] == lime_type || 
		limeEntryTypes["etmc-propagator"] == lime_type
		) ? 1 : 0;
}

bool checkLimeEntryForBinaryData(std::string lime_type, LimeEntryTypes limeEntryTypes)
{
	return	( 
		limeEntryTypes["scidac binary data"] == lime_type || 
		limeEntryTypes["ildg binary data"] == lime_type
		) ? true : false;
}

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

LimeFileReader::LimeFileReader(std::string filenameIn, int precision) : 
	LimeFileReader_basic(filenameIn), desiredPrecision(precision)
{
	checkIfFileExists(filename);
	
	extractMetadataFromLimeFile();
}

void LimeFileReader::extractDataFromLimeFile(char ** destination, size_t expectedNumberOfBytes)
{
	logger.trace() << "Reading data from LIME file \"" << filename << "\"...";
	readLimeFile(destination, expectedNumberOfBytes);
	logger.trace() << "\tsuccesfully read data from LIME file " << filename;
}

char* createBuffer(size_t datasize)
{
	char * buffer = new char[datasize];
	return buffer;
}

static void checkBufferSize(size_t actualSize, size_t expectedSize)
{
	if(actualSize != expectedSize) {
		throw Invalid_Parameters("Binary data does not have expected size.", expectedSize, actualSize);
	}
}

void LimeFileReader::extractBinaryDataFromLimeEntry( LimeHeaderData limeHeaderData, char ** destination, size_t expectedNumberOfBytes)
{
	if( checkLimeEntryForBinaryData(limeHeaderData.limeEntryType, limeEntryTypes) )
	{
		if (limeFileProp.numberOfBinaryDataEntries == 1)
		{
			//todo: generalize for diff. field types..
			checkBufferSize(limeHeaderData.numberOfBytes, expectedNumberOfBytes);
			destination[limeFileProp.numberOfBinaryDataEntries-1] = createBuffer(limeHeaderData.numberOfBytes);
			limeReaderReadData(destination[limeFileProp.numberOfBinaryDataEntries-1], &limeHeaderData.numberOfBytes, limeReader);
		}
		else
		{
			logger.fatal() << "Reading of more than one binary data entry from LIME file not yet implemented...";
		}
	}
}

void checkPrecision(int desiredPrecision, int actualPrecision) throw(Print_Error_Message)
{
	if(desiredPrecision != actualPrecision) 
		throw Print_Error_Message("\nThe desired precision and the one from the sourcefile do not match. Aborting", __FILE__, __LINE__);
}

void LimeFileReader::extractMetadataFromLimeFile()
{
	readMetaDataFromLimeFile();
	
	this->parameters.printMetaDataToScreen(filename);
	
	//todo: this may be unified with a check against the inputparameters..
	checkPrecision(desiredPrecision, this->parameters.prec);
}

void LimeFileReader::readMetaDataFromLimeFile()
{
	logger.trace() << "Reading metadata from LIME file \"" << filename << "\"...";
	readLimeFile(NULL);
	limeFileProp.readMetaData = true;
	logger.trace() << "\tsuccesfully read metadata from LIME file " << filename;
}

void LimeFileReader::readLimeFile(char ** destination, size_t expectedNumberOfBytes)
{
	//TODO: this construction is not nice, but currently necessary as the file is potentially read multiple times
	openFile();
	
	goThroughLimeRecords(destination, expectedNumberOfBytes);
	
	closeFile();
}

void checkLimeRecordReadForFailure(int returnValueFromLimeRecordRead)
{
	if( returnValueFromLimeRecordRead != LIME_SUCCESS ) {
		std::ostringstream errorMessage;
		errorMessage << "\t\tlimeReaderNextRecord returned status = "  << returnValueFromLimeRecordRead;
		throw Print_Error_Message( errorMessage.str(), __FILE__, __LINE__);
	}
}

void LimeFileReader::goThroughLimeRecords(char ** destination, size_t expectedNumberOfBytes)
{
	int statusOfLimeReader = LIME_SUCCESS;
	while( (statusOfLimeReader = limeReaderNextRecord(limeReader)) != LIME_EOF ) 
	{
		checkLimeRecordReadForFailure(statusOfLimeReader);
		extractInformationFromLimeEntry(destination, expectedNumberOfBytes);
	}
}

void LimeFileReader::extractInformationFromLimeEntry(char ** destination, size_t expectedNumberOfBytes)
{
	LimeHeaderData limeHeaderData(limeReader);
	logger.trace() << "Found entry in LIME file of type \"" + limeHeaderData.limeEntryType + "\"";
	if (limeHeaderData.MB_flag != LIME_ERR_MBME) 
	{
		if( ! limeFileProp.readMetaData )
		{
			this->limeFileProp += extractMetaDataFromLimeEntry(limeHeaderData);
		}
		else 
		{
			extractBinaryDataFromLimeEntry(limeHeaderData, destination, expectedNumberOfBytes);
		}
	}
	else
	{
		logger.fatal() << "Error while reading LIME entry. MB flag is \"" + boost::lexical_cast<std::string>(limeHeaderData.MB_flag) + "\"";
	}
}

static void logger_readLimeEntry(std::string type)
{
	logger.trace() << "\tfound \"" << type << "\" entry...";
}

static void logger_readLimeEntrySuccess()
{
	logger.trace() << "\t...succesfully read entry";
}

char * createBufferAndReadLimeDataIntoIt(LimeReader * r, size_t nbytes)
{
	char * buffer = new char[nbytes + 1];
	int error = limeReaderReadData(buffer, &nbytes, r);
	if(error != 0) 
		throw Print_Error_Message("Error while reading data from LIME file...", __FILE__, __LINE__);
	buffer[nbytes] = 0;
	return buffer;
}

void LimeFileReader::handleLimeEntry_xlf(Sourcefileparameters & parameters, char * buffer, std::string lime_type)
{
	if ( LimeFileReader::limeFileProp.numberOfFermionicEntries > 1 )
	{
		logger.warn() << "Reading more than one fermion field is not implemented yet! Skip this entry...";
		return;
	}
	logger_readLimeEntry( lime_type);
	
	sourcefileParameters::setFromLimeEntry_xlf(parameters, buffer);
}

void LimeFileReader::handleLimeEntry_ildg(Sourcefileparameters & parameters, char * buffer, std::string lime_type, size_t numberOfBytes)
{
	logger_readLimeEntry( lime_type);
	sourcefileParameters::setFromLimeEntry_ildg(parameters, buffer, numberOfBytes);
}

void LimeFileReader::handleLimeEntry_scidacChecksum(char * buffer, std::string lime_type, size_t numberOfBytes)
{
	logger_readLimeEntry( lime_type );
	sourcefileParameters::setFromLimeEntry_scidacChecksum(parameters, buffer, numberOfBytes);
}

void LimeFileReader::handleLimeEntry_inverter(std::string lime_type) throw(std::logic_error)
{
	if ( LimeFileReader::limeFileProp.numberOfFermionicEntries > 1 )
	{
		logger.warn() << "Reading more than one fermion field is not implemented yet! Skip this entry...";
		return;
	}
	else
	{
		logger_readLimeEntry( lime_type );
		throw std::logic_error("parsing of inverter infos is not implemented yet. Aborting...");
		//todo: implement similar to xlf infos parsing
		//parameters: "solver", " epssq", " noiter", " kappa", "mu", " time", " hmcversion", " date"};
		
		//todo: this should be moved elsewhere!
		this->parameters.numberOfFermionFieldsRead++;
	}
}

void LimeFileReader::handleLimeEntry_etmcPropagator(std::string lime_type) throw(std::logic_error)
{
	if ( LimeFileReader::limeFileProp.numberOfFermionicEntries > 1 )
	{
		logger.warn() << "Reading more than one fermion field is not implemented yet! Skip this entry...";
		return;
	}
	throw std::logic_error("Reading of etmc propagator not yet implemented. Aborting...");
	logger_readLimeEntry( lime_type );
}

LimeFileProperties LimeFileReader::extractMetaDataFromLimeEntry(LimeHeaderData limeHeaderData)
{
	logger.trace() << "Extracting meta data from LIME entry...";
	LimeFileProperties props;

	if ( checkLimeEntryForBinaryData(limeHeaderData.limeEntryType, limeEntryTypes) == 1)
	{
		props.numberOfBinaryDataEntries = 1;		
	}
	//todo: create class for the different cases
	else
	{
		props.numberOfFermionicEntries += checkLimeEntryForFermionInformations(limeHeaderData.limeEntryType, limeEntryTypes);

		char * buffer = createBufferAndReadLimeDataIntoIt(limeReader,  limeHeaderData.numberOfBytes);
		std::string lime_type = limeHeaderData.limeEntryType;

		if(limeEntryTypes["scidac checksum"] == lime_type) 
		{
			handleLimeEntry_scidacChecksum(buffer, lime_type, limeHeaderData.numberOfBytes); 
		}
		else if( limeEntryTypes["inverter"] == lime_type )
		{
			handleLimeEntry_inverter(lime_type);
		}
		else if(limeEntryTypes["xlf"]  == lime_type ) 
		{
			handleLimeEntry_xlf(this->parameters, buffer, lime_type);
		}
		else if ( limeEntryTypes["ildg"] == lime_type )
		{
			handleLimeEntry_ildg(this->parameters, buffer, lime_type, limeHeaderData.numberOfBytes);
		}
		else if( limeEntryTypes["etmc propagator"] == lime_type )
		{
			handleLimeEntry_etmcPropagator(lime_type);
		}
		else
		{
			logger.warn() << "Do not know LIME entry type \"" + lime_type + "\"";
		}
		
		delete[] buffer;
		logger_readLimeEntrySuccess();
	}

	props.numberOfEntries = 1;

	return props;
}

int LimeFileReader::getReadTrajectoryNumber() noexcept
{
	return this->parameters.trajectorynr;
}

double LimeFileReader::getReadPlaquetteValue() noexcept
{
	return this->parameters.plaquettevalue;
}