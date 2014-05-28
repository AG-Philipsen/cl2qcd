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

#include "parser_utils.h"

#define ENDIAN (htons(1) == 1)

//todo: use these instead of hardcoded strings.
//todo: move to .h ?
//possible limeEntryTypes
const std::vector<std::string> limeEntryTypes = {
  "propagator-type", "xlf-info", "inverter-info", "gauge-scidac-checksum-copy", "etmc-propagator-format",
  "scidac-binary-data", "scidac-checksum", "ildg-format", "ildg-binary-data"
};

void sourcefileparameters::get_XLF_infos(const char * filename)
{
	char hmcversion[50];
	char date[50];
	
	FILE * reader;
	reader = fopen(filename, "r");
	if (reader != NULL) {
		//there are " " inf front of most of the labels, the last one was added for a different style
		const char * tmparray [] = {"plaquette", " trajectory nr", " beta", "kappa", "mu", "c2_rec", " time", " hmcversion", " mubar", " epsilonbar", " date", " plaquette"};
		char tmp1[512];

		while ( fgets (tmp1, 512, reader) != NULL) {
			trim2(tmp1);
			if(strncmp(tmparray[0], tmp1, strlen(tmparray[0])) == 0) extrInfo_hmc_float(tmp1, strlen(tmparray[0]), strlen(tmp1), &plaquettevalue_source);
			if(strncmp(tmparray[1], tmp1, strlen(tmparray[1])) == 0) extrInfo_int(tmp1, strlen(tmparray[1]), strlen(tmp1), &trajectorynr_source);
			if(strncmp(tmparray[2], tmp1, strlen(tmparray[2])) == 0) extrInfo_beta(tmp1, strlen(tmparray[2]), strlen(tmp1), &beta_source, &kappa_source, &mu_source, &c2_rec_source);
			if(strncmp(tmparray[6], tmp1, strlen(tmparray[6])) == 0) extrInfo_int(tmp1, strlen(tmparray[6]), strlen(tmp1), &time_source);
			if(strncmp(tmparray[7], tmp1, strlen(tmparray[7])) == 0) extrInfo_char(tmp1, strlen(tmparray[7]), strlen(tmp1), hmcversion);
			if(strncmp(tmparray[8], tmp1, strlen(tmparray[8])) == 0) extrInfo_hmc_float(tmp1, strlen(tmparray[8]), strlen(tmp1), &mubar_source);
			if(strncmp(tmparray[9], tmp1, strlen(tmparray[9])) == 0) extrInfo_hmc_float(tmp1, strlen(tmparray[9]), strlen(tmp1), &epsilonbar_source);
			if(strncmp(tmparray[10], tmp1, strlen(tmparray[10])) == 0) extrInfo_char(tmp1, strlen(tmparray[10]), strlen(tmp1), date);
			if(strncmp(tmparray[11], tmp1, strlen(tmparray[11])) == 0) extrInfo_hmc_float(tmp1,  strlen(tmparray[11]), strlen(tmp1), &plaquettevalue_source);
		}
	} else throw File_Exception(filename);

	hmcversion_source = hmcversion;
	date_source = date;
	
	logger.trace() << "\tsuccesfully read XLFInfos";
	return;
}

void sourcefileparameters::get_inverter_infos(const char * filename, char * solver, char * hmcversion, char * date )
{
	FILE * reader;
	reader = fopen(filename, "r");
	if (reader != NULL) {
		//there are " " inf front of most of the labels
		const char * tmparray [] = {"solver", " epssq", " noiter", " kappa", "mu", " time", " hmcversion", " date"};
		char tmp1[512];
		while ( fgets (tmp1, 512, reader) != NULL ) {
			if(strncmp(tmparray[0], tmp1, strlen(tmparray[0])) == 0) extrInfo_char(tmp1, strlen(tmparray[0]), strlen(tmp1), solver);
			if(strncmp(tmparray[1], tmp1, strlen(tmparray[1])) == 0) extrInfo_hmc_float(tmp1,  strlen(tmparray[1]), strlen(tmp1), &epssq_source);
			if(strncmp(tmparray[2], tmp1, strlen(tmparray[2])) == 0) extrInfo_int (tmp1, strlen(tmparray[2]), strlen(tmp1), &noiter_source);
			if(strncmp(tmparray[3], tmp1, strlen(tmparray[3])) == 0) extrInfo_kappa(tmp1, strlen(tmparray[3]), strlen(tmp1), &kappa_solver_source, &mu_solver_source);
			if(strncmp(tmparray[5], tmp1, strlen(tmparray[5])) == 0) extrInfo_int (tmp1, strlen(tmparray[5]), strlen(tmp1), &time_solver_source);
			if(strncmp(tmparray[6], tmp1, strlen(tmparray[6])) == 0) extrInfo_char(tmp1, strlen(tmparray[6]), strlen(tmp1), hmcversion);
			if(strncmp(tmparray[7], tmp1, strlen(tmparray[7])) == 0) extrInfo_char(tmp1, strlen(tmparray[7]), strlen(tmp1), date);
		}
	} else throw File_Exception(filename);

	logger.trace() << "\tsuccesfully read InverterInfos" ;
	return;
}

void sourcefileparameters::get_XML_infos(const char * buffer, int size)
{
	xmlTextReaderPtr reader;
	int returnValue;
	int tmpArray[6];
	char field[100];

	//This fct. takes a URL as 3rd argument, but this seems to be unimportant here, so it is left out.
	//See http://xmlsoft.org/html/libxml-xmlreader.html#xmlReaderForMemory
	reader = xmlReaderForMemory(buffer, size, NULL, NULL, 0);

	if (reader != NULL) {
		returnValue = xmlTextReaderRead(reader);
		while (returnValue == 1) {
			get_XML_info_simple(reader, tmpArray, field);
			returnValue = xmlTextReaderRead(reader);
		}
		if (returnValue == -1) {
		  throw Print_Error_Message( "There was an error the XML parser...", __FILE__, __LINE__);
		}
		xmlFreeTextReader(reader);
	} else throw Print_Error_Message( "There was an error the XML parser...", __FILE__, __LINE__);

	prec_source = tmpArray[0];
	flavours_source = tmpArray[1];
	lx_source = tmpArray[2];
	ly_source = tmpArray[3];
	lz_source = tmpArray[4];
	lt_source = tmpArray[5];
	field_source = field;
}

//todo: merge with xml fcts. above!!
Checksum sourcefileparameters::get_checksum(const char * buffer, int size)
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

	return Checksum(suma, sumb);
}

// read in binary file and save it as readable file
// since tmLQCD always saves data with BigEndian one has to be careful

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

int sourcefileparameters::calcNumberOfEntriesBasedOnFieldType(char * fieldType)
{
  if(strcmp(fieldType, "diracFermion") == 0) {
    return calcNumberOfEntriesForDiracFermionfield();
  } else if(strcmp(fieldType, "su3gauge") == 0) {
    return calcNumberOfEntriesForGaugefield();
  } else {
    throw Print_Error_Message("Unknown ildg field type...", __FILE__, __LINE__);
    return 0; //to get rid of a warning
  }
}

int sourcefileparameters::calcNumberOfEntriesBasedOnFieldType(std::string fieldType)
{
  if(fieldType == "diracFermion") {
    return calcNumberOfEntriesForDiracFermionfield();
  } else if( fieldType == "su3gauge") {
    return calcNumberOfEntriesForGaugefield();
  } else {
    throw Print_Error_Message("Unknown ildg field type...", __FILE__, __LINE__);
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

char *createBufferAndReadLimeDataIntoIt(LimeReader * r, size_t nbytes)
{
  char * buffer = new char[nbytes + 1];
  int error = limeReaderReadData(buffer, &nbytes, r);
  if(error != 0) 
    throw Print_Error_Message("Something went wrong...", __FILE__, __LINE__);
  return buffer;
}

void createTemporaryFileToStoreStreamFromLimeReader(char * tmp_file_name, LimeReader *r, size_t nbytes)
{
  char * buffer = 0;
  FILE * tmp;

  buffer = createBufferAndReadLimeDataIntoIt(r, nbytes);
  tmp = fopen(tmp_file_name, "w");
  if(tmp == NULL) 
    throw Print_Error_Message("\t\terror in creating tmp file\n", __FILE__, __LINE__);
  fwrite(buffer, 1, nbytes, tmp);
  fclose(tmp);
  delete [] buffer;
}

void sourcefileparameters::checkLimeEntryForInverterInfos(std::string lime_type, int switcher, LimeReader *r, size_t nbytes)
{
  if("inverter-info" == lime_type)
    {
      if ( switcher > 1 )
	{
	  logger.fatal() << "Reading more than one fermion field is not implemented yet. Aborting...";
	  return;
	}
      
      char solvertype[50];
      char hmcversion_solver[50];
      char date_solver[50];
	  
      logger.trace() << "\tfound inverter-infos as lime_type " << lime_type ;
      numberOfFermionFieldsRead++;
      
      char tmp_file_name[] = "tmpfilenameone";
      createTemporaryFileToStoreStreamFromLimeReader(tmp_file_name, r, nbytes);
      
      get_inverter_infos(tmp_file_name, solvertype, hmcversion_solver, date_solver);
      
      remove(tmp_file_name);
      
      solvertype_source =  solvertype;
      hmcversion_solver_source = hmcversion_solver;
      date_solver_source = date_solver;
    }
}

void sourcefileparameters::checkLimeEntryForXlfInfos(std::string lime_type, int switcher, LimeReader *r, size_t nbytes)
{
  //!!read XLF info, only FIRST fermion is read!!
  if("xlf-info" == lime_type && switcher < 2) {
    FILE * tmp;

    logger.trace() << "\tfound XLF-infos as lime_type " << lime_type;
    char tmp_file_name[] = "tmpfilenametwo";
    createTemporaryFileToStoreStreamFromLimeReader(tmp_file_name, r, nbytes);
    
    get_XLF_infos(tmp_file_name);
    
    remove(tmp_file_name);
  }
}

void sourcefileparameters::checkLimeEntryForXlmInfos(std::string lime_type, int switcher, LimeReader *r, size_t nbytes)
{
  //!!read ildg format (gauge fields) or etmc-propagator-format (fermions), only FIRST fermion is read!!
  if(("etmc-propagator-format" == lime_type || "ildg-format" == lime_type) && switcher < 2 ) {
    logger.trace() << "\tfound XML-infos as lime_type \"" << lime_type << "\"";

		char * buffer = createBufferAndReadLimeDataIntoIt(r, nbytes);
		get_XML_infos(buffer, nbytes);
		delete[] buffer;
    logger.trace() << "\tsuccesfully read XMLInfos";
    
    num_entries_source = calcNumberOfEntriesBasedOnFieldType(field_source);
  }
}

void sourcefileparameters::checkLimeEntryForScidacChecksum(std::string lime_type, LimeReader *r, size_t nbytes)
{
  if("scidac-checksum" == lime_type) {
    logger.trace() << "\tfound scidac-checksum as lime_type" << lime_type;
    char * buffer = createBufferAndReadLimeDataIntoIt(r, nbytes);
    buffer[nbytes] = 0;

    checksum = get_checksum(buffer, nbytes);
    delete[] buffer;
  }
}

int checkLimeEntryForFermionInformations(std::string lime_type, int switcher)
{
  if("propagator-type" == lime_type) {
    if (switcher == 0) {
      logger.info() << "\tfile contains fermion informations" ;
    } else {
      logger.info() << "\tfile contains informations about more than one fermion: " << switcher;
    }
    return 1;
  }
  else 
    return 0;
}

void sourcefileparameters::checkLimeEntry(int *  numberOfFermionEntries, LimeReader * r, LimeHeaderData limeHeaderData)
{
  *numberOfFermionEntries += checkLimeEntryForFermionInformations(limeHeaderData.limeEntryType, *numberOfFermionEntries);
    
  checkLimeEntryForInverterInfos(limeHeaderData.limeEntryType, *numberOfFermionEntries, r, limeHeaderData.numberOfBytes);
  
  checkLimeEntryForXlfInfos(limeHeaderData.limeEntryType, *numberOfFermionEntries, r, limeHeaderData.numberOfBytes);
  
  checkLimeEntryForXlmInfos(limeHeaderData.limeEntryType, *numberOfFermionEntries, r, limeHeaderData.numberOfBytes);
  
  checkLimeEntryForScidacChecksum(limeHeaderData.limeEntryType, r, limeHeaderData.numberOfBytes);
}

void sourcefileparameters::checkPrecision(int desiredPrecision)
{
  if(desiredPrecision != prec_source) 
    throw Print_Error_Message("\nThe desired precision and the one from the sourcefile do not match, will not read data!!!", __FILE__, __LINE__);
}

size_t sourcefileparameters::sizeOfGaugefieldBuffer()
{
  return num_entries_source * sizeof(hmc_float);
}

char* sourcefileparameters::createBufferForGaugefield(int num_entries)
{
  size_t datasize = sizeOfGaugefieldBuffer();
  char * buffer = new char[datasize];
  return buffer;
}

void sourcefileparameters::checkSizeOfBinaryDataForGaugefield(size_t actualSize)
{
  size_t expectedSize = sizeOfGaugefieldBuffer();
  if(actualSize != expectedSize) {
    throw Invalid_Parameters("Binary data does not have expected size.", expectedSize, actualSize);
  }
}

int sourcefileparameters::extractBinaryDataFromLimeEntry_NeedsDifferentName(LimeReader * r, LimeHeaderData limeHeaderData, char ** destination, int numberOfBinaryDataEntries)
{
  if (numberOfBinaryDataEntries == 0)
    {
      //todo: generalize for diff. field types..
      checkSizeOfBinaryDataForGaugefield(limeHeaderData.numberOfBytes);
      destination[numberOfBinaryDataEntries] = createBufferForGaugefield(num_entries_source);
      limeReaderReadData(destination[numberOfBinaryDataEntries], &limeHeaderData.numberOfBytes, r);
    }
  else
    {
      logger.fatal() << "Reading of more than one binary data entry from LIME file not yet implemented...";
    }
  return (numberOfBinaryDataEntries)+1;
}

LimeFileProperties sourcefileparameters::extractBinaryDataFromLimeEntry(LimeReader * r, char ** destination, int numberOfBinaryDataEntries)
{
	LimeHeaderData limeHeaderData(r);
	LimeFileProperties limeFileProp;
	
	if (limeHeaderData.MB_flag == 1) 
	{
		if( 
				limeEntryTypes[5] == limeHeaderData.limeEntryType || 
				limeEntryTypes[8] == limeHeaderData.limeEntryType
			)
		{
			logger.fatal() << limeHeaderData.limeEntryType;
			limeFileProp.numberOfBinaryDataEntries += extractBinaryDataFromLimeEntry_NeedsDifferentName(r, limeHeaderData, destination, numberOfBinaryDataEntries);
		}
		limeFileProp.numberOfEntries += 1;
	}
	
	return limeFileProp;
}

LimeFileProperties sourcefileparameters::extractMetaDataFromLimeEntry(LimeReader * r)
{
	LimeFileProperties limeFileProp;
	LimeHeaderData limeHeaderData(r);

	if (limeHeaderData.MB_flag == 1) 
	{
		//todo: fix this!
		int numberOfFermionEntries = 0;
		checkLimeEntry(&numberOfFermionEntries, r, limeHeaderData);
		limeFileProp.numberOfEntries += 1;
	}

	return limeFileProp;
}

LimeFileProperties sourcefileparameters::extractInformationFromLimeEntry(LimeReader * r, char ** destination, bool readMetaData, int numberOfBinaryDataEntries)
{
  if( readMetaData)
    {
      return extractMetaDataFromLimeEntry(r);
    }
  else 
    {
      return extractBinaryDataFromLimeEntry(r, destination, numberOfBinaryDataEntries);
    }
}

void sourcefileparameters::goThroughLimeRecords(LimeReader * r, char ** destination, bool readMetaData)
{
  LimeFilePropertiesCollector limeFileProp;
  int statusOfLimeReader = 0;

  while( (statusOfLimeReader = limeReaderNextRecord(r)) != LIME_EOF ) {
    checkLimeRecordReadForFailure(statusOfLimeReader);
    limeFileProp += extractInformationFromLimeEntry(r, destination, readMetaData, limeFileProp.numberOfBinaryDataEntries);
  }
}

void sourcefileparameters::readLimeFile(std::string sourceFilename, char ** destination, bool readMetaData)
{
  FILE *limeFileOpenedForReading;
  LimeReader *limeReader;

  limeFileOpenedForReading = fopen (sourceFilename.c_str(), "r");
  limeReader = limeCreateReader(limeFileOpenedForReading);

  goThroughLimeRecords(limeReader, destination, readMetaData);

  limeDestroyReader(limeReader);
  fclose(limeFileOpenedForReading); 
}

void sourcefileparameters::readDataFromLimeFile(std::string sourceFilename, char ** destination)
{
  logger.trace() << "Reading data from LIME file \"" << sourceFilename << "\"...";
  readLimeFile(sourceFilename, destination, false);
  logger.trace() << "\tsuccesfully read data from LIME file " << sourceFilename;
}

void sourcefileparameters::readMetaDataFromLimeFile(std::string sourceFilename)
{
  logger.trace() << "Reading metadata from LIME file \"" << sourceFilename << "\"...";
  readLimeFile(sourceFilename, NULL, true);
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

void sourcefileparameters::printMetaDataToScreen(std::string file)
{
  logger.info() << "*************************************************************" ;
  logger.info() << "*************************************************************" ;
  logger.info() << "Metadata from file " << file << ":";
  logger.info() << "\treading XML-data gave:";
  logger.info() << "\t\tfield type:\t" << field_source ;
  logger.info() << "\t\tprecision:\t" << prec_source ;
  logger.info() << "\t\tlx:\t\t" << lx_source ;
  logger.info() << "\t\tly:\t\t" << ly_source ;
  logger.info() << "\t\tlz:\t\t" << lz_source ;
  logger.info() << "\t\tlt:\t\t" << lt_source ;
  logger.info() << "\t\tflavours:\t" << flavours_source ;
  logger.info() << "\treading XLF-data gave:";
  logger.info() << "\t\tplaquette:\t" << plaquettevalue_source;
  logger.info() << "\t\ttrajectorynr:\t" << trajectorynr_source;
  logger.info() << "\t\tbeta:\t\t" << beta_source;
  logger.info() << "\t\tkappa:\t\t" << kappa_source;
  logger.info() << "\t\tmu:\t\t" << mu_source;
  logger.info() << "\t\tc2_rec:\t\t" << c2_rec_source;
  logger.info() << "\t\ttime:\t\t" << time_source;
  logger.info() << "\t\thmc-version:\t" << hmcversion_source;
  logger.info() << "\t\tmubar:\t\t" << mubar_source;
  logger.info() << "\t\tepsilonbar:\t" << epsilonbar_source;
  logger.info() << "\t\tdate:\t\t" << date_source;
  if(numberOfFermionFieldsRead != 0) {
    logger.info() << "\treading inverter-data gave:";
    logger.info() << "\t\tsolvertype:\t" << solvertype_source;
    logger.info() << "\t\tepssq:\t\t" << std::setprecision(30) << epssq_source;
    logger.info() << "\t\tnoiter:\t\t" << noiter_source;
    logger.info() << "\t\tkappa_solver:\t" << kappa_solver_source;
    logger.info() << "\t\tmu_solver:\t" << mu_solver_source;
    logger.info() << "\t\ttime_solver:\t" << time_solver_source;
    logger.info() << "\t\thmc-ver_solver:\t" << hmcversion_solver_source;
    logger.info() << "\t\tdate_solver:\t" << date_solver_source;
  }
  logger.info() << "\tfile-checksum:\t" << checksum;
  logger.info() << "*************************************************************" ;
}

void sourcefileparameters::set_defaults()
{
	lx_source = 0;
	ly_source = 0;
	lz_source = 0;
	lt_source = 0;
	prec_source = 0;
	num_entries_source = 0;
	flavours_source = 0;
	trajectorynr_source = 0;
	time_source = 0;
	time_solver_source = 0;
	noiter_source = 0;
	plaquettevalue_source = 0;
	beta_source = 0;
	kappa_source = 0;
	mu_source = 0;
	c2_rec_source = 0;
	mubar_source = 0;
	epsilonbar_source = 0;
	epssq_source = 0;
	kappa_solver_source = 0;
	mu_solver_source = 0;
	
	return;
}

void sourcefileparameters::extractMetadataFromLimeFile(std::string sourceFilename, int desiredPrecision)
{
  readMetaDataFromLimeFile(sourceFilename);

  printMetaDataToScreen(sourceFilename);

  //todo: this may be unified with a check against the inputparameters..
  checkPrecision(desiredPrecision);  
}

void sourcefileparameters::extractDataFromLimeFile(std::string sourceFilename, char ** destination)
{
  readDataFromLimeFile(sourceFilename, destination);
  //todo: put conversion to numbers in here...
}

//todo: make char ** std::vector<char*>
void sourcefileparameters::readsourcefile(std::string sourceFilename, int desiredPrecision, char ** destination)
{
  checkIfFileExists(sourceFilename);

  extractMetadataFromLimeFile(sourceFilename, desiredPrecision);
   
  extractDataFromLimeFile(sourceFilename, destination);
}

LimeFilePropertiesCollector:: ~LimeFilePropertiesCollector()
{
  logger.trace() << "Found " << numberOfEntries << " LIME records.";
  logger.trace() << "Found " << numberOfBinaryDataEntries << " binary entries in LIME file";
}

void LimeFileProperties::operator+=(LimeFileProperties other)
{
  this->numberOfEntries += other.numberOfEntries;
  this->numberOfBinaryDataEntries += other.numberOfBinaryDataEntries;
}
