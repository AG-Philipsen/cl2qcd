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

#include "ildgIo_gaugefield.hpp"

#include "limeUtilities.hpp"

#include "../host_functionality/logger.hpp"
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
	helperMap["plaquette"].re = boost::regex ("plaquette\\s+=\\s+[\\+\\-]*\\d\\.\\d+");
	helperMap["trajectory_nr"].re = boost::regex ("trajectory nr\\s+=\\s+[\\+\\-]*\\d+");
	helperMap["beta"].re = boost::regex ("beta\\s+=\\s+[\\+\\-]*\\d.\\d+");
	helperMap["kappa"].re = boost::regex ("kappa\\s+=\\s+[\\+\\-]*\\d.\\d+");
	helperMap["mu"].re = boost::regex ("mu\\s+=\\s+[\\+\\-]*\\d.\\d+");
	helperMap["c2_rec"].re= boost::regex ("c2_rec\\s+=\\s+[\\+\\-]*\\d.\\d+");
 	helperMap["time"].re = boost::regex ("time\\s+=\\s+[\\+\\-]*\\d+");
 	helperMap["hmcversion"].re = boost::regex ("hmcversion\\s+=\\s+\\d.\\d+[a-z]*");
	helperMap["mubar"].re = boost::regex ("mubar\\s+=\\s+[\\+\\-]*\\d.\\d+");
	helperMap["epsilonbar"].re = boost::regex ("epsilonbar\\s+=\\s+[\\+\\-]*\\d.\\d+");
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
int calcNumberOfEntriesForDiracFermionfield(const sourcefileparameters_values params)
{
  //latSize sites, 4 dirac indices, Nc colour indices, 2 complex indices
  return (int) (params.lx_source) * (params.ly_source) * (params.lz_source) * (params.lt_source) * NC * NSPIN * 2;
}
int calcNumberOfEntriesForGaugefield(const sourcefileparameters_values params)
{
  // latSize sites, 4 links, 2 complex indices -> 9 complex numbers per link
  return (int) (params.lx_source) * (params.ly_source) * (params.lz_source) * (params.lt_source) * 2 * 4 * 9;
}

int sourcefileparameters::calcNumberOfEntriesBasedOnFieldType(std::string fieldType)
{
	if(fieldType == "diracFermion") 
	{
		return calcNumberOfEntriesForDiracFermionfield( *this );
	} 
	else if( fieldType == "su3gauge") 
	{
		return calcNumberOfEntriesForGaugefield( *this);
	} 
	else 
	{
		throw Print_Error_Message("Unknown ildg field type \"" + fieldType + "\"", __FILE__, __LINE__);
	}
	return 0; //to get rid of a warning
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

static void setParametersToValues_ildg(sourcefileparameters & parameters, std::map <std::string, std::string> helperMap)
{
	parameters.prec_source = castStringToInt(helperMap["precision"]);
	parameters.lx_source = castStringToInt(helperMap["lx"]);
	parameters.ly_source = castStringToInt(helperMap["ly"]);
	parameters.lz_source = castStringToInt(helperMap["lz"]);
	parameters.lt_source = castStringToInt(helperMap["lt"]);

	parameters.field_source = helperMap["field"];
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

static void setParametersToValues_scidacChecksum(sourcefileparameters & parameters, std::map <std::string, std::string> helperMap)
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

int sourcefileparameters::checkLimeEntryForFermionInformations(std::string lime_type)
{
	return ( 
		limeEntryTypes["propagator"] == lime_type || 
		limeEntryTypes["inverter"] == lime_type || 
		limeEntryTypes["etmc-propagator"] == lime_type
		) ? 1 : 0;
}

bool sourcefileparameters::checkLimeEntryForBinaryData(std::string lime_type)
{
	return	( 
		limeEntryTypes["scidac binary data"] == lime_type || 
		limeEntryTypes["ildg binary data"] == lime_type
		) ? true : false;
}

LimeFileProperties sourcefileparameters::extractMetaDataFromLimeEntry(LimeReader * r, LimeHeaderData limeHeaderData)
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
		if ( sourcefileparameters::limeFileProp.numberOfFermionicEntries > 1 )
		{
			logger.warn() << "Reading more than one fermion field is not implemented yet!";
		}

		char * buffer = createBufferAndReadLimeDataIntoIt(r,  limeHeaderData.numberOfBytes);
		std::string lime_type = limeHeaderData.limeEntryType;

		if(limeEntryTypes["scidac checksum"] == lime_type) 
		{
			logger_readLimeEntry( lime_type );

			std::map<std::string, std::string> helperMap;
			fillHelperMap_scidacChecksum(helperMap);
			goThroughBufferWithXmlReaderAndExtractInformationBasedOnMap(buffer, limeHeaderData.numberOfBytes, helperMap);
			setParametersToValues_scidacChecksum(*this, helperMap);		  
		}
		
		else if( limeEntryTypes["inverter"] == lime_type && sourcefileparameters::limeFileProp.numberOfFermionicEntries < 2)
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
				numberOfFermionFieldsRead++;
			}
		}
		
		else if(limeEntryTypes["xlf"]  == lime_type && sourcefileparameters::limeFileProp.numberOfFermionicEntries < 2) 
		{
			logger_readLimeEntry( lime_type);
			std::string str(buffer);
			std::map<std::string, helper> helperMap;			
			fillHelperMap_xlf(helperMap);

			mapStringToHelperMap(str, helperMap);
			setParametersToValues_xlf(*this, helperMap);
		}

		else if ( limeEntryTypes["ildg"] == lime_type )
		{
			logger_readLimeEntry( lime_type);
			std::map<std::string, std::string> helperMap;
			fillHelperMap_ildg(helperMap);
			
			goThroughBufferWithXmlReaderAndExtractInformationBasedOnMap(buffer, limeHeaderData.numberOfBytes, helperMap);
			
			setParametersToValues_ildg(*this, helperMap);
			num_entries_source = calcNumberOfEntriesBasedOnFieldType(field_source);
		}	

		else if( limeEntryTypes["etmc propagator"] == lime_type && sourcefileparameters::limeFileProp.numberOfFermionicEntries < 2 )
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

void sourcefileparameters::extractBinaryDataFromLimeEntry(LimeReader * r, LimeHeaderData limeHeaderData, char ** destination)
{
	if( checkLimeEntryForBinaryData(limeHeaderData.limeEntryType) )
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
}

void sourcefileparameters::extractInformationFromLimeEntry(LimeReader * r, char ** destination)
{
	LimeHeaderData limeHeaderData(r);
	logger.trace() << "Found entry in LIME file of type \"" + limeHeaderData.limeEntryType + "\"";
	if (limeHeaderData.MB_flag != LIME_ERR_MBME) 
	{
		if( ! limeFileProp.readMetaData )
		{
			this->limeFileProp += extractMetaDataFromLimeEntry(r, limeHeaderData);
		}
		else 
		{
			extractBinaryDataFromLimeEntry(r, limeHeaderData, destination);
		}
	}
	else
	{
		logger.fatal() << "Error while reading LIME entry. MB flag is \"" + boost::lexical_cast<std::string>(limeHeaderData.MB_flag) + "\"";
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

//todo: this must be readsourcefile_gaugefield or so, and then one has to check if the entry is in fact "su3gauge"
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

//todo: refactor!

#include <assert.h>

#include "../host_functionality/logger.hpp"
#include <sstream>

#include <time.h>

void write_gaugefield (
  char * binary_data, n_uint64_t num_bytes, Checksum checksum,
  int lx, int ly, int lz, int lt, int prec, int trajectorynr, hmc_float plaquettevalue, hmc_float beta, hmc_float kappa, hmc_float mu, hmc_float c2_rec, hmc_float epsilonbar, hmc_float mubar,
  const char * hmc_version, const char * filename)
{

	logger.info() << "writing gaugefield to lime-file...";

	time_t current_time;
	FILE *outputfile;
	outputfile = fopen(filename, "w");
	int MB_flag;
	int ME_flag;
	n_uint64_t length_xlf_info = 0, length_ildg_format = 0, length_scidac_checksum = 0;
	LimeRecordHeader * header_ildg_format, *header_scidac_checksum, * header_ildg_binary_data, * header_xlf_info;

	//set values
	const char * field_out = "su3gauge";
	time_t rawtime;
	time ( &current_time );
	const char * date = ctime (&current_time);

	// TODO replace this whole block by something templated
	//get binary data
	//here it must not be assumed that the argument prec and sizeof(hmc_float) are the same!!
	logger.debug() << "  num_bytes = " << num_bytes;

	if(sizeof(hmc_float) * 8 != prec) {
		throw Invalid_Parameters("Precision does not match executables.", sizeof(hmc_float) * 8, prec);
	}

	//write xlf-info to string, should look like this
	/*
	char xlf_info [] = "plaquette = 6.225960e-01\n trajectory nr = 67336\n beta = 6.000000, kappa = 0.177000, mu = 0.500000, c2_rec = 0.000000\n time = 1278542490\n hmcversion = 5.1.5\n mubar = 0.000000\n epsilonbar = 0.000000\n date = Thu Jul  8 00:41:30 2010\n";
	*/

	//note: it is not clear what "time" is supposed to be
	char dummystring[1000];
	char xlf_info[1000];
	sprintf(xlf_info, "%s", "plaquette = ");
	sprintf(dummystring, "%f", plaquettevalue);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, "\n trajectory nr = ");
	sprintf(dummystring, "%i", trajectorynr);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, "\n beta = ");
	sprintf(dummystring, "%f", beta);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, ", kappa = ");
	sprintf(dummystring, "%f", kappa);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, ", mu = ");
	sprintf(dummystring, "%f", mu);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, ", c2_rec = ");
	sprintf(dummystring, "%f", c2_rec);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, "\n time = ");
	sprintf(dummystring, "%i", current_time);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, "\n hmcversion = ");
	sprintf(dummystring, "%s", hmc_version);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, "\n mubar = ");
	sprintf(dummystring, "%f", mubar);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, "\n epsilonbar = ");
	sprintf(dummystring, "%f", epsilonbar);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, "\n date = ");
	sprintf(dummystring, "%s", date);
	strcat(xlf_info, dummystring);

	length_xlf_info = strlen(xlf_info);

	//write scidac checksum, this is stubb
	std::string scidac_checksum;
	{
		std::ostringstream tmp;
		tmp << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<scidacChecksum>\n<version>1.0</version>\n";
		tmp << "<suma>" << std::hex << checksum.get_suma() << "</suma>\n";
		tmp << "<sumb>" << std::hex << checksum.get_sumb() << "</sumb>\n";
		tmp << "</scidacChecksum>";
		scidac_checksum = tmp.str();
	}
	length_scidac_checksum = scidac_checksum.length();

	//write ildg_format to string, should look like this:
	/*
	char ildg_format [] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n  <version>1.0</version>\n  <field>su3gauge</field>\n  <precision>64</precision>\n  <lx>4</lx>\n  <ly>4</ly>\n  <lz>4</lz>\n  <lt>4</lt>\n</ildgFormat>";
	*/
	char dummystring2[1000];
	char ildg_format[1000];
	sprintf(ildg_format, "%s", "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n  <version>1.0</version>\n  <field>");
	strcat(ildg_format, field_out);
	strcat(ildg_format, "</field>\n  <precision>");
	sprintf(dummystring2, "%i", prec);
	strcat(ildg_format, dummystring2);
	strcat(ildg_format, "</precision>\n  <lx>");
	sprintf(dummystring2, "%i", lx);
	strcat(ildg_format, dummystring2);
	strcat(ildg_format, "</lx>\n  <ly>");
	sprintf(dummystring2, "%i", ly);
	strcat(ildg_format, dummystring2);
	strcat(ildg_format, "</ly>\n  <lz>");
	sprintf(dummystring2, "%i", lz);
	strcat(ildg_format, dummystring2);
	strcat(ildg_format, "</lz>\n  <lt>");
	sprintf(dummystring2, "%i", lt);
	strcat(ildg_format, dummystring2);
	strcat(ildg_format, "</lt>\n</ildgFormat>");

	length_ildg_format = strlen(ildg_format);


	//write the lime file
	MB_flag = 1;

	LimeWriter *writer;
	writer = limeCreateWriter(outputfile);

	const char * types[] = {"xlf-info", "ildg-format", "ildg-binary-data", "scidac-checksum"};

	//xlf-info
	ME_flag = 1;
	header_xlf_info = limeCreateHeader(MB_flag, ME_flag, (char*) types[0], length_xlf_info);
	limeWriteRecordHeader(header_xlf_info, writer);
	limeDestroyHeader(header_xlf_info);
	limeWriteRecordData( xlf_info, &length_xlf_info, writer);
	logger.debug() << "  xlf-info written";

	//ildg-format
	ME_flag = 2;
	header_ildg_format = limeCreateHeader(MB_flag, ME_flag, (char*) types[1], length_ildg_format);
	limeWriteRecordHeader(header_ildg_format, writer);
	limeDestroyHeader(header_ildg_format);
	limeWriteRecordData( ildg_format, &length_ildg_format, writer);
	logger.debug() << "  ildg-format written";

	//binary data
	ME_flag = 3;
	header_ildg_binary_data = limeCreateHeader(MB_flag, ME_flag, (char*) types[2], num_bytes);
	limeWriteRecordHeader(header_ildg_binary_data, writer);
	limeDestroyHeader(header_ildg_binary_data);
	limeWriteRecordData(binary_data, &num_bytes, writer);
	logger.debug() << "  ildg_binary_data written";

	//scidac-checksum
	ME_flag = 4;
	header_scidac_checksum = limeCreateHeader(MB_flag, ME_flag, (char*) types[3], length_scidac_checksum);
	limeWriteRecordHeader(header_scidac_checksum, writer);
	limeDestroyHeader(header_scidac_checksum);
	limeWriteRecordData(const_cast<char*>(scidac_checksum.c_str()), &length_scidac_checksum, writer);
	logger.debug() << "  scidac-checksum written";

	//closing
	fclose(outputfile);
	limeDestroyWriter(writer);
	logger.info() << "  " << (float) ( (float) (length_xlf_info + length_ildg_format + num_bytes + length_scidac_checksum) / 1024 / 1024 ) << " MBytes were written to the lime file " << filename;
}

