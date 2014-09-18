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

LimeFileReader::LimeFileReader(std::string sourceFilenameIn, int precision, char ** data) : 
	sourceFilename(sourceFilenameIn), desiredPrecision(precision)
{
	checkIfFileExists(sourceFilename);
	
	extractMetadataFromLimeFile();
	
	extractDataFromLimeFile(data);
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

void LimeFileReader::extractDataFromLimeFile(char ** destination)
{
	readDataFromLimeFile(destination);
	//todo: put conversion to numbers in here...
}

void LimeFileReader::readDataFromLimeFile(char ** destination)
{
	logger.trace() << "Reading data from LIME file \"" << sourceFilename << "\"...";
	readLimeFile(destination);
	logger.trace() << "\tsuccesfully read data from LIME file " << sourceFilename;
}


size_t LimeFileReader::sizeOfGaugefieldBuffer()
{
	return this->parameters.num_entries * sizeof(hmc_float);
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

void LimeFileReader::extractBinaryDataFromLimeEntry( LimeHeaderData limeHeaderData, char ** destination)
{
	if( checkLimeEntryForBinaryData(limeHeaderData.limeEntryType) )
	{
		if (limeFileProp.numberOfBinaryDataEntries == 1)
		{
			//todo: generalize for diff. field types..
			checkBufferSize(limeHeaderData.numberOfBytes, sizeOfGaugefieldBuffer());
			destination[limeFileProp.numberOfBinaryDataEntries-1] = createBuffer(sizeOfGaugefieldBuffer());
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
	
	this->parameters.printMetaDataToScreen(sourceFilename);
	
	//todo: this may be unified with a check against the inputparameters..
	checkPrecision(desiredPrecision, this->parameters.prec);
}

void LimeFileReader::readMetaDataFromLimeFile()
{
	logger.trace() << "Reading metadata from LIME file \"" << sourceFilename << "\"...";
	readLimeFile(NULL);
	limeFileProp.readMetaData = true;
	logger.trace() << "\tsuccesfully read metadata from LIME file " << sourceFilename;
}

void LimeFileReader::readLimeFile(char ** destination)
{
	//TODO: this construction is not nice, but currently necessary as the file is potentially read multiple times
	openFile();
	
	goThroughLimeRecords(destination);
	
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

void LimeFileReader::goThroughLimeRecords(char ** destination)
{
	int statusOfLimeReader = LIME_SUCCESS;
	while( (statusOfLimeReader = limeReaderNextRecord(limeReader)) != LIME_EOF ) 
	{
		checkLimeRecordReadForFailure(statusOfLimeReader);
		extractInformationFromLimeEntry(destination);
	}
}

void LimeFileReader::extractInformationFromLimeEntry(char ** destination)
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
			extractBinaryDataFromLimeEntry(limeHeaderData, destination);
		}
	}
	else
	{
		logger.fatal() << "Error while reading LIME entry. MB flag is \"" + boost::lexical_cast<std::string>(limeHeaderData.MB_flag) + "\"";
	}
}

extern "C" {
#include <lime_fixed_types.h>
#include <libxml/parser.h>
#include <libxml/xmlreader.h>
#include <stdio.h>
#include <string.h>
//this is for htons
#include <arpa/inet.h>
}

#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <boost/regex.hpp>
#include <vector>
#include <algorithm>
#include <boost/lexical_cast.hpp>

#define ENDIAN (htons(1) == 1)

static std::string removeNewlines(std::string in)
{
  in.erase(std::remove(in.begin(), in.end(), '\n'), in.end());
  return in;
}

static std::string trimStringBeforeEqual(std::string in)
{
	unsigned found = in.find_last_of("=");
	std::string withoutEqual = in.substr(found+1);
	found = ( withoutEqual ).find_first_of(" ");
	return  withoutEqual.substr(found+1);
}

static void logger_readLimeEntry(std::string type)
{
	logger.trace() << "\tfound \"" << type << "\" entry...";
}

static void logger_readLimeEntrySuccess()
{
	logger.trace() << "\t...succesfully read entry";
}

class helper
{
public:
	boost::regex re;
	std::string str;
	std::string value;
};

static double castStringToDouble(std::string in)
{
	return boost::lexical_cast<double>(in);
}

static int castStringToInt(std::string in)
{
	return boost::lexical_cast<int>(in);
}

static void fillHelperMap_xlf(std::map<std::string, helper> & helperMap)
{
	helperMap["plaquette"].re = boost::regex ("plaquette\\s+=\\s+[\\+\\-]*\\d+\\.\\d+");
	helperMap["trajectory_nr"].re = boost::regex ("trajectory nr\\s+=\\s+[\\+\\-]*\\d+");
	helperMap["beta"].re = boost::regex ("beta\\s+=\\s+[\\+\\-]*\\d+.\\d+");
	helperMap["kappa"].re = boost::regex ("kappa\\s+=\\s+[\\+\\-]*\\d+.\\d+");
	helperMap["mu"].re = boost::regex ("mu\\s+=\\s+[\\+\\-]*\\d+.\\d+");
	helperMap["c2_rec"].re= boost::regex ("c2_rec\\s+=\\s+[\\+\\-]*\\d+.\\d+");
 	helperMap["time"].re = boost::regex ("time\\s+=\\s+[\\+\\-]*\\d+");
 	helperMap["hmcversion"].re = boost::regex ("hmcversion\\s+=\\s+\\d.\\d+[a-z]*");
	helperMap["mubar"].re = boost::regex ("mubar\\s+=\\s+[\\+\\-]*\\d+.\\d+");
	helperMap["epsilonbar"].re = boost::regex ("epsilonbar\\s+=\\s+[\\+\\-]*\\d+.\\d+");
 	helperMap["date"].re = boost::regex ("date\\s+=\\s+[\\s\\.a-zA-Z\\d\\:]+");
}

static void setParametersToValues_xlf(Sourcefileparameters & parameters, std::map<std::string, helper>  helperMap)
{
	parameters.plaquettevalue = castStringToDouble(helperMap["plaquette"].value);
	parameters.kappa = castStringToDouble(helperMap["kappa"].value);
	parameters.mu = castStringToDouble(helperMap["mu"].value);
	parameters.beta = castStringToDouble(helperMap["beta"].value);
	parameters.c2_rec = castStringToDouble(helperMap["c2_rec"].value);
	parameters.mubar = castStringToDouble(helperMap["mubar"].value);
	parameters.epsilonbar = castStringToDouble(helperMap["epsilonbar"].value);

	parameters.trajectorynr = castStringToInt(helperMap["trajectory_nr"].value);
	parameters.time = castStringToInt(helperMap["time"].value);

	parameters.hmcversion = helperMap["hmcversion"].value;
	parameters.date = helperMap["date"].value;
}

static void mapStringToHelperMap(std::string str, std::map<std::string, helper> &  helperMap)
{
	logger.trace() << "Going through string:";
	logger.trace() << str;
	
	//todo: find out about ::iterator
	for (std::map<std::string, helper>::iterator it = helperMap.begin(); it != helperMap.end(); it++)
	{
		logger.trace() << "Found \"" + it->first + "\"";
		
		//todo: add check if re is found
		
		//todo: make this more beautiful
		//http://stackoverflow.com/questions/10058606/c-splitting-a-string-by-a-character
		//start/end points of tokens in str
		std::vector<std::string> tokens;
		boost::sregex_token_iterator begin(str.begin(), str.end(), it->second.re), end;
		std::copy(begin, end, std::back_inserter(tokens));

		it->second.str = removeNewlines( tokens[0] );
		it->second.value =  trimStringBeforeEqual(it->second.str);
	}
}

const int xmlNodeType_startElement = 1;

void extractXmlValuesBasedOnMap(xmlTextReaderPtr reader, std::map<std::string, std::string> * map)
{
	xmlChar *name, *value;
	name = xmlTextReaderName(reader);
	int type = xmlTextReaderNodeType(reader);
	std::string tmp(reinterpret_cast<char*> (name) );

	for (std::map<std::string, std::string>::iterator it = map->begin(); it != map->end(); it++)
	{
		if(it->first == tmp && type == xmlNodeType_startElement)
		{
			xmlTextReaderRead(reader);
			value = xmlTextReaderValue(reader);
			if(value != nullptr) {
				std::string tmp2 (reinterpret_cast<char*> (value) );
				it->second = tmp2;
				xmlFree(value);
				logger.trace() << "Found xml value:\t" + tmp2;
			}
		}
	}

	xmlFree(name);
}

void goThroughBufferWithXmlReaderAndExtractInformationBasedOnMap(const char * buffer, int size, std::map<std::string, std::string> & helperMap)
{
	xmlTextReaderPtr reader;
	int returnValue;

	//This fct. takes a URL as 3rd argument, but this seems to be unimportant here, so it is left out.
	//See http://xmlsoft.org/html/libxml-xmlreader.html#xmlReaderForMemory
	reader = xmlReaderForMemory(buffer, size, NULL, NULL, 0);

	if (reader != NULL) {
		returnValue = xmlTextReaderRead(reader);
		while (returnValue == 1) 
		{
			extractXmlValuesBasedOnMap(reader, &helperMap);
			returnValue = xmlTextReaderRead(reader);
		}
		if (returnValue == -1) 
		{
			throw Print_Error_Message( "There was an error in the XML parser...", __FILE__, __LINE__);
		}
		xmlFreeTextReader(reader);
	} 
	else 
		throw Print_Error_Message( "There was an error in the XML parser...", __FILE__, __LINE__);
}

//NOTE: these two functions are similar to some in the meta package,
//      but I would rather not include the latter here.
int calcNumberOfEntriesForDiracFermionfield(const Sourcefileparameters params)
{
  //latSize sites, 4 dirac indices, Nc colour indices, 2 complex indices
  return (int) (params.lx) * (params.ly) * (params.lz) * (params.lt) * NC * NSPIN * 2;
}
int calcNumberOfEntriesForGaugefield(const Sourcefileparameters params)
{
  // latSize sites, 4 links, 2 complex indices -> 9 complex numbers per link
  return (int) (params.lx) * (params.ly) * (params.lz) * (params.lt) * 2 * 4 * 9;
}

int LimeFileReader::calcNumberOfEntriesBasedOnFieldType(std::string fieldType)
{
	if(fieldType == "diracFermion") 
	{
		return calcNumberOfEntriesForDiracFermionfield( this->parameters );
	} 
	else if( fieldType == "su3gauge") 
	{
		return calcNumberOfEntriesForGaugefield( this->parameters );
	} 
	else 
	{
		throw Print_Error_Message("Unknown ildg field type \"" + fieldType + "\"", __FILE__, __LINE__);
	}
	return 0; //to get rid of a warning
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

static void setParametersToValues_ildg(Sourcefileparameters & parameters, std::map <std::string, std::string> helperMap)
{
	parameters.prec = castStringToInt(helperMap["precision"]);
	parameters.lx = castStringToInt(helperMap["lx"]);
	parameters.ly = castStringToInt(helperMap["ly"]);
	parameters.lz = castStringToInt(helperMap["lz"]);
	parameters.lt = castStringToInt(helperMap["lt"]);

	parameters.field = helperMap["field"];
}

static void fillHelperMap_ildg(std::map<std::string, std::string> & helperMap)
{
	helperMap["field"] = "";
	helperMap["precision"] = "";
	helperMap["lx"] = "";
	helperMap["ly"] = "";
	helperMap["lz"] = "";
	helperMap["lt"] = "";
}

static void fillHelperMap_scidacChecksum(std::map<std::string, std::string> & helperMap)
{
	helperMap["suma"] = "";
	helperMap["sumb"] = "";
}

static void setParametersToValues_scidacChecksum(Sourcefileparameters & parameters, std::map <std::string, std::string> helperMap)
{
	uint32_t suma, sumb;
	
	std::stringstream tmp, tmp2;
	tmp << helperMap["suma"];
	tmp >> std::hex >> suma;
	tmp2 << helperMap["sumb"];
	tmp2 >> std::hex >> sumb;
	logger.fatal() << suma;
	logger.fatal() << sumb;
	
	parameters.checksum =  Checksum(suma, sumb);
}

int LimeFileReader::checkLimeEntryForFermionInformations(std::string lime_type)
{
	return ( 
		limeEntryTypes["propagator"] == lime_type || 
		limeEntryTypes["inverter"] == lime_type || 
		limeEntryTypes["etmc-propagator"] == lime_type
		) ? 1 : 0;
}

bool LimeFileReader::checkLimeEntryForBinaryData(std::string lime_type)
{
	return	( 
		limeEntryTypes["scidac binary data"] == lime_type || 
		limeEntryTypes["ildg binary data"] == lime_type
		) ? true : false;
}

LimeFileProperties LimeFileReader::extractMetaDataFromLimeEntry(LimeHeaderData limeHeaderData)
{
	logger.trace() << "Extracting meta data from LIME entry...";
	LimeFileProperties props;

	if ( checkLimeEntryForBinaryData(limeHeaderData.limeEntryType) == 1)
	{
		props.numberOfBinaryDataEntries = 1;		
	}
	//todo: create class for the different cases
	else
	{
		//todo: is this meaningful?
		props.numberOfFermionicEntries += checkLimeEntryForFermionInformations(limeHeaderData.limeEntryType);
		if ( LimeFileReader::limeFileProp.numberOfFermionicEntries > 1 )
		{
			logger.warn() << "Reading more than one fermion field is not implemented yet!";
		}

		char * buffer = createBufferAndReadLimeDataIntoIt(limeReader,  limeHeaderData.numberOfBytes);
		std::string lime_type = limeHeaderData.limeEntryType;

		if(limeEntryTypes["scidac checksum"] == lime_type) 
		{
			logger_readLimeEntry( lime_type );

			std::map<std::string, std::string> helperMap;
			fillHelperMap_scidacChecksum(helperMap);
			goThroughBufferWithXmlReaderAndExtractInformationBasedOnMap(buffer, limeHeaderData.numberOfBytes, helperMap);
			setParametersToValues_scidacChecksum(this->parameters, helperMap);		  
		}
		
		else if( limeEntryTypes["inverter"] == lime_type && LimeFileReader::limeFileProp.numberOfFermionicEntries < 2)
		{
			if ( limeFileProp.numberOfFermionicEntries > 1 )
			{
				logger.fatal() << "Reading more than one fermion field is not implemented yet. Aborting...";
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
		
		else if(limeEntryTypes["xlf"]  == lime_type && LimeFileReader::limeFileProp.numberOfFermionicEntries < 2) 
		{
			logger_readLimeEntry( lime_type);
			std::string str(buffer);
			std::map<std::string, helper> helperMap;
			fillHelperMap_xlf(helperMap);

			mapStringToHelperMap(str, helperMap);
			setParametersToValues_xlf(this->parameters, helperMap);
		}

		else if ( limeEntryTypes["ildg"] == lime_type )
		{
			logger_readLimeEntry( lime_type);
			std::map<std::string, std::string> helperMap;
			fillHelperMap_ildg(helperMap);
			
			goThroughBufferWithXmlReaderAndExtractInformationBasedOnMap(buffer, limeHeaderData.numberOfBytes, helperMap);
			
			setParametersToValues_ildg(this->parameters, helperMap);
			this->parameters.num_entries = calcNumberOfEntriesBasedOnFieldType(this->parameters.field);
		}	

		else if( limeEntryTypes["etmc propagator"] == lime_type && LimeFileReader::limeFileProp.numberOfFermionicEntries < 2 )
		{
			throw std::logic_error("Reading of etmc propagator not yet implemented. Aborting...");
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