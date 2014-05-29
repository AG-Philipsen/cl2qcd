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

#include "ildg_read_gaugefield.h"

#include "limeUtilities.hpp"

#include "../host_functionality/logger.hpp"
#include <sstream>

#include <iostream>
#include <vector>

extern "C" {
#include <lime_fixed_types.h>
#include <libxml/parser.h>
#include <libxml/xmlreader.h>
#include <stdio.h>
#include <string.h>
//this is for htons
#include <arpa/inet.h>
}

//#include "parser_utils.h"

#include <string>
#include <algorithm>    // copy
#include <iterator>     // back_inserter
#include <boost/regex.hpp>        // regex, sregex_token_iterator
#include <vector>
#include <algorithm>
#include <boost/lexical_cast.hpp>

#define ENDIAN (htons(1) == 1)

class LimeEntryTypes
{
public:
	typedef std::map<std::string, std::string> Mapper;
	static Mapper mapper;
	std::string operator[](std::string key)
	{
		return mapper[key];
	}
};

//todo: move to some .cpp file
LimeEntryTypes::Mapper LimeEntryTypes::mapper = { {"propagator", "propagator-info"}, {"xlf", "xlf-info"} , {"inverter", "inverter-info"}, {"gauge-checksum-copy", "gauge-scidac-checksum-copy"}, {"etmc-propagator", "etmc-propagator-format"},  { "scidac binary data", "scidac-binary-data"}, {"scidac checksum", "scidac-checksum"}, {"ildg", "ildg-format"}, {"ildg binary data", "ildg-binary-data"}  };

LimeEntryTypes limeEntryTypes;

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
	helperMap["plaquette"].re = boost::regex ("plaquette\\s+=\\s+\\d\\.\\d+");
	helperMap["trajectory_nr"].re = boost::regex ("trajectory nr\\s+=\\s+\\d+");
	helperMap["beta"].re = boost::regex ("beta\\s+=\\s+\\d.\\d+");
	helperMap["kappa"].re = boost::regex ("kappa\\s+=\\s+\\d.\\d+");
	helperMap["mu"].re = boost::regex ("mu\\s+=\\s+\\d.\\d+");
	helperMap["c2_rec"].re= boost::regex ("c2_rec\\s+=\\s+\\d.\\d+");
 	helperMap["time"].re = boost::regex ("time\\s+=\\s+[\\+\\-]\\d+");
 	helperMap["hmcversion"].re = boost::regex ("hmcversion\\s+=\\s+\\d.\\d+[a-z]*");
	helperMap["mubar"].re = boost::regex ("mubar\\s+=\\s+\\d.\\d+");
	helperMap["epsilonbar"].re = boost::regex ("epsilonbar\\s+=\\s+\\d.\\d+");
 	helperMap["date"].re = boost::regex ("date\\s+=\\s+[\\s\\.a-zA-Z\\d\\:]+");
}

static void setParametersToValues_xlf(sourcefileparameters & parameters, std::map<std::string, helper>  helperMap)
{
	parameters.plaquettevalue_source = castStringToDouble(helperMap["plaquette"].value);
	parameters.kappa_source = castStringToDouble(helperMap["kappa"].value);
	parameters.mu_source = castStringToDouble(helperMap["mu"].value);
	parameters.beta_source = castStringToDouble(helperMap["beta"].value);
	parameters.c2_rec_source = castStringToDouble(helperMap["c2_rec"].value);
	parameters.mubar_source = castStringToDouble(helperMap["mubar"].value);
	parameters.epsilonbar_source = castStringToDouble(helperMap["epsilonbar"].value);

	parameters.trajectorynr_source = castStringToInt(helperMap["trajectory_nr"].value);
	parameters.time_source = castStringToInt(helperMap["time"].value);

	parameters.hmcversion_source = helperMap["hmcversion"].value;
	parameters.date_source = helperMap["date"].value;
}

static void mapStringToHelperMap(std::string str, std::map<std::string, helper> &  helperMap)
{
	//todo: find out about ::iterator
	for (std::map<std::string, helper>::iterator it = helperMap.begin(); it != helperMap.end(); it++)
	{
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

void get_XLF_infos(const char * buffer, sourcefileparameters & parameters)
{
	std::string str(buffer);
	//todo: make class out of this
	std::map<std::string, helper> helperMap;

	fillHelperMap_xlf(helperMap);
	mapStringToHelperMap(str, helperMap);
	setParametersToValues_xlf(parameters, helperMap);
}

void get_inverter_infos(const char * buffer, sourcefileparameters & parameters)
{
	std::string str(buffer);
	logger.fatal() << str;
	throw std::logic_error("parsing of inverter infos is not implemented yet. Aborting...");
	//todo: implement similar to xlf infos parsing
	//parameters: "solver", " epssq", " noiter", " kappa", "mu", " time", " hmcversion", " date"};
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
			std::string tmp2 (reinterpret_cast<char*> (value) );
			it->second = tmp2;
			xmlFree(value);
		}

	}

	xmlFree(name);
}


void get_XML_infos(const char * buffer, int size, std::map<std::string, std::string> & helperMap)
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

//todo: merge with xml fcts. above!!
void get_checksum(const char * buffer, int size, sourcefileparameters & parameters)
{
	uint32_t suma, sumb;

	// see comment in get_XML_infos(..)
	xmlTextReaderPtr reader = xmlReaderForMemory(buffer, size, NULL, NULL, 0);
	if(reader == NULL) {
	  throw Print_Error_Message( "There was an error the XML parser...", __FILE__, __LINE__);
	}

	int xml_state = xmlTextReaderRead(reader);
	while(xml_state == 1) {
		std::stringstream name_buf;
		xmlChar* name_xml = xmlTextReaderName(reader);
		name_buf << name_xml;
		free(name_xml);
		std::string name = name_buf.str();


		if("suma" == name) {
			xmlTextReaderRead(reader);
			xmlChar * value = xmlTextReaderValue(reader);
			if(value != nullptr) {
				std::stringstream tmp;
				tmp << value;
				tmp >> std::hex >> suma;
				free(value);
			}
			xmlTextReaderRead(reader);
		}
		if("sumb" == name) {
			xmlTextReaderRead(reader);
			xmlChar * value = xmlTextReaderValue(reader);
			if(value != nullptr) {
				std::stringstream tmp;
				tmp << value;
				tmp >> std::hex >> sumb;
				free(value);
			}
			xmlTextReaderRead(reader);
		}

		xml_state = xmlTextReaderRead(reader);
	}
	if(xml_state == -1) {
	  throw Print_Error_Message( "There was an error the XML parser...", __FILE__, __LINE__);
	}
	xmlFreeTextReader(reader);

	parameters.checksum =  Checksum(suma, sumb);
}

//NOTE: these two functions are similar to some in the meta package,
//      but I would rather not include the latter here.
int sourcefileparameters::calcNumberOfEntriesForDiracFermionfield()
{
  //latSize sites, 4 dirac indices, Nc colour indices, 2 complex indices
  return (int) (lx_source) * (ly_source) * (lz_source) * (lt_source) * NC * NSPIN * 2;
}
int sourcefileparameters::calcNumberOfEntriesForGaugefield()
{
  // latSize sites, 4 links, 2 complex indices -> 9 complex numbers per link
  return (int) (lx_source) * (ly_source) * (lz_source) * (lt_source) * 2 * 4 * 9;
}

int sourcefileparameters::calcNumberOfEntriesBasedOnFieldType(std::string fieldType)
{
  if(fieldType == "diracFermion") {
    return calcNumberOfEntriesForDiracFermionfield();
  } else if( fieldType == "su3gauge") {
    return calcNumberOfEntriesForGaugefield();
  } else {
	  throw Print_Error_Message("Unknown ildg field type \"" + fieldType + "\"", __FILE__, __LINE__);
    return 0; //to get rid of a warning
  }
}

void checkLimeRecordReadForFailure(int returnValueFromLimeRecordRead)
{
	if( returnValueFromLimeRecordRead != LIME_SUCCESS ) {
		std::ostringstream errorMessage;
		errorMessage << "\t\tlimeReaderNextRecord returned status = "  << returnValueFromLimeRecordRead;
		throw Print_Error_Message( errorMessage.str(), __FILE__, __LINE__);
	}
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

void sourcefileparameters::checkLimeEntryForInverterInfos(std::string lime_type, LimeReader *r, size_t nbytes)
{
	if( limeEntryTypes["inverter"] == lime_type)
	{
		if ( limeFileProp.numberOfFermionicEntries > 1 )
		{
			logger.fatal() << "Reading more than one fermion field is not implemented yet. Aborting...";
			return;
		}
	
		logger.trace() << "\tfound inverter-infos as lime_type " << lime_type ;
		
		char * buffer = createBufferAndReadLimeDataIntoIt(r, nbytes);
		get_inverter_infos(buffer, *this);
		delete[] buffer;
		logger_readLimeEntrySuccess();

		//todo: this should be moved elsewhere!
		numberOfFermionFieldsRead++;
	}
}

void sourcefileparameters::checkLimeEntryForXlfInfos(std::string lime_type, LimeReader *r, size_t nbytes)
{
	//!!read XLF info, only FIRST fermion is read!!
	if(limeEntryTypes["xlf"]  == lime_type && limeFileProp.numberOfFermionicEntries < 2) 
	{
		logger.trace() << "\tfound XLF-infos as lime_type " << lime_type;
		char * buffer = createBufferAndReadLimeDataIntoIt(r, nbytes);
		get_XLF_infos(buffer, *this);
		delete[] buffer;
		logger_readLimeEntrySuccess();
	}
}

static void setParametersToValues_xlm(sourcefileparameters & parameters, std::map <std::string, std::string> helperMap)
{
	parameters.prec_source = castStringToInt(helperMap["precision"]);
	parameters.lx_source = castStringToInt(helperMap["lx"]);
	parameters.ly_source = castStringToInt(helperMap["ly"]);
	parameters.lz_source = castStringToInt(helperMap["lz"]);
	parameters.lt_source = castStringToInt(helperMap["lt"]);
	parameters.flavours_source = castStringToInt(helperMap["flavours"]);

	parameters.field_source = helperMap["field"];
}

static void fillHelperMap_xml(std::map<std::string, std::string> & helperMap)
{
	helperMap["field"] = "";
	helperMap["precision"] = "";
	//todo: this is not in su3gauge entries! Distinct between cases
	helperMap["flavours"] = "0";
	helperMap["lx"] = "";
	helperMap["ly"] = "";
	helperMap["lz"] = "";
	helperMap["lt"] = "";
}

void sourcefileparameters::checkLimeEntryForXlmInfos(std::string lime_type, LimeReader *r, size_t nbytes)
{
	//!!read ildg format (gauge fields) or etmc-propagator-format (fermions), only FIRST fermion is read!!
	if
		(
			(limeEntryTypes["etmc propagator"] == lime_type || limeEntryTypes["ildg"] == lime_type) &&
			limeFileProp.numberOfFermionicEntries < 2 
			)
	{

		std::map<std::string, std::string> helperMap;
		fillHelperMap_xml(helperMap);

		logger.trace() << "\tfound XML-infos as lime_type";
		char * buffer = createBufferAndReadLimeDataIntoIt(r, nbytes);
		get_XML_infos(buffer, nbytes, helperMap);
		delete[] buffer;
		
		setParametersToValues_xlm(*this, helperMap);
		num_entries_source = calcNumberOfEntriesBasedOnFieldType(field_source);

		logger_readLimeEntrySuccess();
	}
}

void sourcefileparameters::checkLimeEntryForScidacChecksum(std::string lime_type, LimeReader *r, size_t nbytes)
{
	if(limeEntryTypes["scidac checksum"] == lime_type) 
	{
		logger_readLimeEntry( lime_type );
		char * buffer = createBufferAndReadLimeDataIntoIt(r, nbytes);
		get_checksum(buffer, nbytes, *this);
		delete[] buffer;
		logger_readLimeEntrySuccess();
	}
}

int checkLimeEntryForFermionInformations(std::string lime_type)
{
	return limeEntryTypes["propagator"] == lime_type ? 1 : 0;
}

int checkLimeEntryForBinaryData(std::string lime_type)
{
	return	( 
		limeEntryTypes["scidac binary data"] == lime_type || 
		limeEntryTypes["ildg binary data"] == lime_type
		) ? 1 : 0;
}

LimeFileProperties sourcefileparameters::extractMetaDataFromLimeEntry(LimeReader * r, LimeHeaderData limeHeaderData)
{
	LimeFileProperties props;
 
	props.numberOfFermionicEntries += checkLimeEntryForFermionInformations(limeHeaderData.limeEntryType);
    
	checkLimeEntryForInverterInfos(limeHeaderData.limeEntryType, r, limeHeaderData.numberOfBytes);
	
	checkLimeEntryForXlfInfos(limeHeaderData.limeEntryType, r, limeHeaderData.numberOfBytes);
	
	checkLimeEntryForXlmInfos(limeHeaderData.limeEntryType, r, limeHeaderData.numberOfBytes);
	
	checkLimeEntryForScidacChecksum(limeHeaderData.limeEntryType, r, limeHeaderData.numberOfBytes);
	
	props.numberOfBinaryDataEntries = checkLimeEntryForBinaryData(limeHeaderData.limeEntryType);
	props.numberOfEntries = 1;
	
	return props;
}

size_t sourcefileparameters::sizeOfGaugefieldBuffer()
{
  return num_entries_source * sizeof(hmc_float);
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

void sourcefileparameters::extractBinaryDataFromLimeEntry_NeedsDifferentName(LimeReader * r, LimeHeaderData limeHeaderData, char ** destination)
{
  if (limeFileProp.numberOfBinaryDataEntries == 1)
    {
      //todo: generalize for diff. field types..
      checkBufferSize(limeHeaderData.numberOfBytes, sizeOfGaugefieldBuffer());
      destination[limeFileProp.numberOfBinaryDataEntries-1] = createBuffer(sizeOfGaugefieldBuffer());
      limeReaderReadData(destination[limeFileProp.numberOfBinaryDataEntries-1], &limeHeaderData.numberOfBytes, r);
    }
  else
    {
      logger.fatal() << "Reading of more than one binary data entry from LIME file not yet implemented...";
    }
}

static void checkIfMetaDataHasBeenRead(bool readMetaData)
{
	if(! readMetaData)
	{
		throw std::logic_error("Did not read metadata yet. Will not read in binary data. Aborting...");
	}
}

void sourcefileparameters::extractBinaryDataFromLimeEntry(LimeReader * r, char ** destination, LimeHeaderData limeHeaderData)
{
	checkIfMetaDataHasBeenRead(limeFileProp.readMetaData);
	if( checkLimeEntryForBinaryData(limeHeaderData.limeEntryType) )
	{
		logger.fatal() << limeHeaderData.limeEntryType;
		extractBinaryDataFromLimeEntry_NeedsDifferentName(r, limeHeaderData, destination);
	}
}

void sourcefileparameters::extractInformationFromLimeEntry(LimeReader * r, char ** destination)
{
	LimeHeaderData limeHeaderData(r);
	if (limeHeaderData.MB_flag == 1) 
	{
		if( ! limeFileProp.readMetaData )
			{
				this->limeFileProp += extractMetaDataFromLimeEntry(r, limeHeaderData);
			}
		else 
			{
				extractBinaryDataFromLimeEntry(r, destination, limeHeaderData);
			}
	}
}

void sourcefileparameters::goThroughLimeRecords(LimeReader * r, char ** destination)
{
	int statusOfLimeReader = LIME_SUCCESS;
	while( (statusOfLimeReader = limeReaderNextRecord(r)) != LIME_EOF ) 
	{
		checkLimeRecordReadForFailure(statusOfLimeReader);
		extractInformationFromLimeEntry(r, destination);
	}
}

void sourcefileparameters::readLimeFile(char ** destination)
{
  FILE *limeFileOpenedForReading;
  LimeReader *limeReader;

  limeFileOpenedForReading = fopen (sourceFilename.c_str(), "r");
  limeReader = limeCreateReader(limeFileOpenedForReading);

  goThroughLimeRecords(limeReader, destination);

  limeDestroyReader(limeReader);
  fclose(limeFileOpenedForReading); 
}

void sourcefileparameters::readDataFromLimeFile(char ** destination)
{
  logger.trace() << "Reading data from LIME file \"" << sourceFilename << "\"...";
  readLimeFile(destination);
  logger.trace() << "\tsuccesfully read data from LIME file " << sourceFilename;
}

void sourcefileparameters::readMetaDataFromLimeFile()
{
  logger.trace() << "Reading metadata from LIME file \"" << sourceFilename << "\"...";
  readLimeFile(NULL);
	limeFileProp.readMetaData = true;
  logger.trace() << "\tsuccesfully read metadata from LIME file " << sourceFilename;
}

void checkIfFileExists(std::string file)
{
  FILE * checker;
  checker = fopen(file.c_str(), "r");
  if(checker == 0) {
    throw File_Exception(file);
  }
  fclose(checker);
  return;
}

void checkPrecision(int desiredPrecision, int actualPrecision)
{
  if(desiredPrecision != actualPrecision) 
    throw Print_Error_Message("\nThe desired precision and the one from the sourcefile do not match. Aborting", __FILE__, __LINE__);
}

void sourcefileparameters::extractMetadataFromLimeFile()
{
  readMetaDataFromLimeFile();

  printMetaDataToScreen(sourceFilename);

  //todo: this may be unified with a check against the inputparameters..
  checkPrecision(desiredPrecision, prec_source);  
}

void sourcefileparameters::extractDataFromLimeFile(char ** destination)
{
  readDataFromLimeFile(destination);
  //todo: put conversion to numbers in here...
}

//todo: make char ** std::vector<char*>
void sourcefileparameters::readsourcefile(std::string sourceFilenameIn, int desiredPrecisionIn, char ** destination)
{
	//todo: move all this to constructor!
  checkIfFileExists(sourceFilenameIn);
	
	sourceFilename = sourceFilenameIn;
	desiredPrecision = desiredPrecisionIn;

  extractMetadataFromLimeFile();
	
  extractDataFromLimeFile(destination);
}
