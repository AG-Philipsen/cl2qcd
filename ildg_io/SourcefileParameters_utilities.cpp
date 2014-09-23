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

#include "SourcefileParameters_utilities.hpp"

#include "../executables/exceptions.h"

#include <iterator>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

extern "C" {
#include <lime_fixed_types.h>
#include <libxml/parser.h>
#include <libxml/xmlreader.h>
#include <stdio.h>
#include <string.h>
//this is for htons
#include <arpa/inet.h>
}

const int xmlNodeType_startElement = 1;

static double castStringToDouble(std::string in)
{
	try
	{
		return boost::lexical_cast<double>(in);
	}
	catch (std::bad_cast)
	{
		throw Print_Error_Message( "Could not cast string \"" + in + "\" to double!", __FILE__, __LINE__);
	}
}

static int castStringToInt(std::string in)
{
	try
	{
		return boost::lexical_cast<int>(in);
	}
	catch (std::bad_cast)
	{
		throw Print_Error_Message( "Could not cast string \"" + in + "\" to int!", __FILE__, __LINE__);
	}
}

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

class ParserMap_xlf
{
public:
	boost::regex re;
	std::string str;
	std::string value = "-1";
};

static std::map<std::string, ParserMap_xlf> createParserMap_xlf()
{
	std::map<std::string, ParserMap_xlf> parserMap;
	
	parserMap["plaquette"].re = boost::regex ("plaquette\\s*=\\s*[\\+\\-]*\\d+[.]*\\d*");
	parserMap["trajectory_nr"].re = boost::regex ("trajectory nr\\s+=\\s+[\\+\\-]*\\d+");
	parserMap["beta"].re = boost::regex ("beta\\s*=\\s*[\\+\\-]*\\d+[.]*\\d*");
	parserMap["kappa"].re = boost::regex ("kappa\\s+=\\s+[\\+\\-]*\\d+[.]*\\d*");
	parserMap["mu"].re = boost::regex ("mu\\s+=\\s+[\\+\\-]*\\d+[.]*\\d*");
	parserMap["c2_rec"].re= boost::regex ("c2_rec\\s+=\\s+[\\+\\-]*\\d+[.]*\\d*");
 	parserMap["time"].re = boost::regex ("time\\s+=\\s+[\\+\\-]*\\d+");
 	parserMap["hmcversion"].re = boost::regex ("hmcversion\\s+=\\s+\\d.+\\d+[a-z]*");
	parserMap["mubar"].re = boost::regex ("mubar\\s+=\\s+[\\+\\-]*\\d+[.]*\\d*");
	parserMap["epsilonbar"].re = boost::regex ("epsilonbar\\s+=\\s+[\\+\\-]*\\d+[.]*\\d*");
 	parserMap["date"].re = boost::regex ("date\\s+=\\s+[\\s\\.a-zA-Z\\d\\:]+");
	
	return parserMap;
}

static void setParametersToValues_xlf(Sourcefileparameters & parameters, std::map<std::string, ParserMap_xlf>  parserMap)
{
	logger.trace() << "setting sourcefile parameters from xlf entry...";
	parameters.plaquettevalue = castStringToDouble(parserMap["plaquette"].value);
	parameters.kappa = castStringToDouble(parserMap["kappa"].value);
	parameters.mu = castStringToDouble(parserMap["mu"].value);
	parameters.beta = castStringToDouble(parserMap["beta"].value);
	parameters.c2_rec = castStringToDouble(parserMap["c2_rec"].value);
	parameters.mubar = castStringToDouble(parserMap["mubar"].value);
	parameters.epsilonbar = castStringToDouble(parserMap["epsilonbar"].value);

	parameters.trajectorynr = castStringToInt(parserMap["trajectory_nr"].value);
	parameters.time = castStringToInt(parserMap["time"].value);

	parameters.hmcversion = parserMap["hmcversion"].value;
	parameters.date = parserMap["date"].value;
}

static void mapStringToParserMap(std::string str, std::map<std::string, ParserMap_xlf> &  parserMap)
{
	logger.trace() << "Going through string:";
	logger.trace() << str;
	
	//todo: find out about ::iterator
	for (std::map<std::string, ParserMap_xlf>::iterator it = parserMap.begin(); it != parserMap.end(); it++)
	{
		logger.trace() << "Check for \"" + it->first + "\"";
		
		//todo: make this more beautiful
		//http://stackoverflow.com/questions/10058606/c-splitting-a-string-by-a-character
		//start/end points of tokens in str
		std::vector<std::string> tokens;
		boost::sregex_token_iterator begin(str.begin(), str.end(), it->second.re), end;
		if(begin == end) //re is not found!
		{
			logger.fatal() << "Did not find match for \"" << it->first << "\"!";
		}
		else
		{
			std::copy(begin, end, std::back_inserter(tokens));
			logger.trace() << "Found match for \"" << it->first << "\":";
			for (int i = 0; i<(int) tokens.size(); i++)
				logger.trace() << tokens[i];
			
			it->second.str = removeNewlines( tokens[0] );
			it->second.value =  trimStringBeforeEqual(it->second.str);
			logger.trace() << "Match for \"" + it->first << "\": " + it->second.str;
			logger.trace() << "Extracted value for \"" + it->first << "\": " + it->second.value;
		}
	}
	logger.trace() << "done...";
}

void sourcefileParameters::setFromLimeEntry_xlf(Sourcefileparameters & parameters, char * buffer)
{
	std::string str(buffer);
	std::map<std::string, ParserMap_xlf> parserMap = createParserMap_xlf();

	mapStringToParserMap(str, parserMap);
	setParametersToValues_xlf(parameters, parserMap);
}

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

void goThroughBufferWithXmlReaderAndExtractInformationBasedOnMap(const char * buffer, int size, std::map<std::string, std::string> & parserMap)
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
			extractXmlValuesBasedOnMap(reader, &parserMap);
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

static void setParametersToValues_ildg(Sourcefileparameters & parameters, std::map <std::string, std::string> parserMap)
{
	parameters.prec = castStringToInt(parserMap["precision"]);
	parameters.lx = castStringToInt(parserMap["lx"]);
	parameters.ly = castStringToInt(parserMap["ly"]);
	parameters.lz = castStringToInt(parserMap["lz"]);
	parameters.lt = castStringToInt(parserMap["lt"]);

	parameters.field = parserMap["field"];
}

static void fillParserMap_ildg(std::map<std::string, std::string> & parserMap)
{
	parserMap["field"] = "";
	parserMap["precision"] = "";
	parserMap["lx"] = "";
	parserMap["ly"] = "";
	parserMap["lz"] = "";
	parserMap["lt"] = "";
}

static void fillParserMap_scidacChecksum(std::map<std::string, std::string> & parserMap)
{
	parserMap["suma"] = "";
	parserMap["sumb"] = "";
}

static void setParametersToValues_scidacChecksum(Sourcefileparameters & parameters, std::map <std::string, std::string> parserMap)
{
	uint32_t suma, sumb;
	
	std::stringstream tmp, tmp2;
	tmp << parserMap["suma"];
	tmp >> std::hex >> suma;
	tmp2 << parserMap["sumb"];
	tmp2 >> std::hex >> sumb;
	
	parameters.checksum =  Checksum(suma, sumb);
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

int calcNumberOfEntriesBasedOnFieldType(const Sourcefileparameters params) throw(Print_Error_Message)
{
	std::string fieldType = params.field;
	
	if(fieldType == "diracFermion") 
	{
		return calcNumberOfEntriesForDiracFermionfield( params );
	} 
	else if( fieldType == "su3gauge") 
	{
		return calcNumberOfEntriesForGaugefield( params );
	} 
	else 
	{
		throw Print_Error_Message("Unknown ildg field type \"" + fieldType + "\"", __FILE__, __LINE__);
	}
	return 0; //to get rid of a warning
}

void sourcefileParameters::setFromLimeEntry_ildg(Sourcefileparameters & parameters, char * buffer, size_t numberOfBytes)
{
	std::map<std::string, std::string> parserMap;
	fillParserMap_ildg(parserMap);
			
	goThroughBufferWithXmlReaderAndExtractInformationBasedOnMap(buffer, numberOfBytes, parserMap);
			
	setParametersToValues_ildg(parameters, parserMap);
	parameters.num_entries = calcNumberOfEntriesBasedOnFieldType(parameters);
}

void sourcefileParameters::setFromLimeEntry_scidacChecksum(Sourcefileparameters & parameters, char * buffer, size_t numberOfBytes)
{
	std::map<std::string, std::string> parserMap;
	fillParserMap_scidacChecksum(parserMap);
	goThroughBufferWithXmlReaderAndExtractInformationBasedOnMap(buffer, numberOfBytes, parserMap);
	setParametersToValues_scidacChecksum(parameters, parserMap);
}
	