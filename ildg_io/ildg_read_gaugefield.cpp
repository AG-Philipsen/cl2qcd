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
#include "checksum.h"


extern "C" {
#include <lime.h>
#include <lime_fixed_types.h>
#include <libxml/parser.h>
#include <libxml/xmlreader.h>
#include <stdio.h>
#include <string.h>
//this is for htons
#include <arpa/inet.h>
}
#include <cassert>

#define ENDIAN (htons(1) == 1)

//todo: use these instead of hardcoded strings.
//todo: move to .h ?
//possible limeEntryTypes
const char * limeEntryTypes[] = {
  "propagator-type", "xlf-info", "inverter-info", "gauge-scidac-checksum-copy", "etmc-propagator-format",
  "scidac-binary-data", "scidac-checksum", "ildg-format", "ildg-binary-data"
};

void extrInfo_hmc_float(const char * in1, int len1, int len2, hmc_float * dest);

// two strings in xlf-info and inverter-info are complicated because there are several vars saved in them
// this is not a beautiful implementation!!
void extrInfo_beta(const char * in1, int len1, int len2, hmc_float * dest1, hmc_float * dest2, hmc_float * dest3, hmc_float * dest4);

void extrInfo_kappa(const char * in1, int len1, int len2, hmc_float * dest1, hmc_float * dest2);

void extrInfo_int(const char * in1, int len1, int len2, int * dest);
// the \n at the end is overwritten by \0
void extrInfo_char(const char * in1, int len1, int len2, char * dest);

void get_XLF_infos(char * filename, hmc_float * plaquettevalue, int * trajectorynr, hmc_float * beta, hmc_float * kappa, hmc_float * mu,
                   hmc_float * c2_rec, int * time, char * hmcversion, hmc_float * mubar, hmc_float * epsilonbar, char * date );

void get_inverter_infos(char * filename, char * solver, hmc_float * epssq, int * noiter, hmc_float * kappa_solver, hmc_float * mu_solver,
                        int * time, char * hmcversion, char * date );

// http://xmlsoft.org/xmlreader.html
// compile with gcc ReadXML.c $(xml2-config --cflags) -Wall $(xml2-config --libs)

void trim(char * buff);

void get_XML_info_simple(xmlTextReaderPtr reader, int numbers[6], char * field);

void get_XML_infos(const char * buffer, int size, const char* filename, int * prec, int * lx, int * ly, int * lz, int *lt, int * flavours, char * field_out );

// get XML Infos: file to be read + parameters
void read_meta_data(const char * file, int * lx, int * ly, int * lz, int * lt, int * prec, char * field_out, int * num_entries,
                    int * flavours, hmc_float * plaquettevalue, int * trajectorynr, hmc_float * beta, hmc_float * kappa, hmc_float * mu, hmc_float * c2_rec, int * time, char * hmcversion, hmc_float * mubar, hmc_float * epsilonbar, char * date,
                    char * solvertype, hmc_float * epssq, int * noiter, hmc_float * kappa_solver, hmc_float * mu_solver,  int * time_solver, char * hmcversion_solver, char * date_solver, int * fermion, Checksum * checksum);

void read_data(const char * file, char * data, size_t bytes);

void read_tmlqcd_file(char * file,
                      int * lx, int * ly, int * lz, int * lt, int * prec, char * field_out, int * num_entries, int * flavours,
                      hmc_float * plaquettevalue, int * trajectorynr, hmc_float * beta, hmc_float * kappa, hmc_float * mu, hmc_float * c2_rec, int * time, char * hmcversion, hmc_float * mubar, hmc_float * epsilonbar, char * date,
                      char * solvertype, hmc_float * epssq, int * noiter, hmc_float * kappa_solver, hmc_float * mu_solver, int * time_solver, char * hmcversion_solver, char * date_solver,
                      hmc_float * array, int * hmc_prec, Checksum * checksum);

Checksum get_checksum(const char * buffer, int size, const char * filename);

//CP:
//changed some things for C++
//   char tmp[len2-len1-3] to "tmp = new char[len2-len1-3; .... free(tmp);" and so forth
//   char arrays are now const char arrays, except where this produces problems with fwrite
//   nasty workaroung by Lars:
//     std::string buffer;
//     limeReaderReadData ((void*) buffer.c_str(),(n_uint64_t *) &nbytes, r);
//     char * buffer2 = new char[nbytes+1];
//     strcpy(buffer2, buffer.c_str());
//     fwrite(buffer2, 1, sizeof(char)*nbytes, tmp);
//   changed the use of tmpnam to some fixed filename

void extrInfo_hmc_float(const char * in1, int len1, int len2, hmc_float * dest)
{
	char * tmp = new char[len2 - len1 - 3 + 1];
	strncpy(tmp, &in1[len1 + 3], len2 - len1 - 3);
	tmp[len2 - len1 - 3] = '\0';
	*dest = atof(tmp);
	delete[] tmp;
}

// two strings in xlf-info and inverter-info are complicated because there are several vars saved in them
// this is not a beautiful implementation!!
void extrInfo_beta(const char * in1, int len1, int len2, hmc_float * dest1, hmc_float * dest2, hmc_float * dest3, hmc_float * dest4)
{
	char * tmp = new char[len2 - len1 - 3 + 1];
	int cutoff;
	strncpy(tmp, &in1[len1 + 3], len2 - len1 - 3);
	tmp[len2 - len1 - 3] = '\0';
	//every number is saved with 6 digits after the "."
	//find the "."
	cutoff = strchr(tmp, '.') - tmp + 1;
	char * beta = new char [cutoff + 6 + 1];
	strncpy(beta, tmp, cutoff + 6);
	beta[cutoff + 6] = '\0';
	*dest1 = atof(beta);
	//cut of the part ", kappa = "
	strcpy(tmp, &tmp[cutoff + 6 + 10]);
	cutoff = strchr(tmp, '.') - tmp + 1;
	char * kappa = new char [cutoff + 6 + 1];
	strncpy(kappa, tmp, cutoff + 6);
	kappa[cutoff + 6] = '\0';
	*dest2 = atof(kappa);
	//cut of the part ", mu = "
	strcpy(tmp, &tmp[cutoff + 6 + 7]);
	cutoff = strchr(tmp, '.') - tmp + 1;
	char * mu = new char [cutoff + 6 + 1];
	strncpy(mu, tmp, cutoff + 6);
	mu[cutoff + 6] = '\0';
	*dest3 = atof(mu);
	//cut of the part ", c2_rec = "
	strcpy(tmp, &tmp[cutoff + 6 + 11]);
	cutoff = strchr(tmp, '.') - tmp + 1;
	char * c2_rec = new char [cutoff + 6 + 1];
	strncpy(c2_rec, tmp, cutoff + 6);
	c2_rec[cutoff + 6] = '\0';
	*dest4 = atof(c2_rec);
	delete [] beta;
	delete [] kappa;
	delete [] c2_rec;
	delete [] mu;
	delete [] tmp;
}

void extrInfo_kappa(const char * in1, int len1, int len2, hmc_float * dest1, hmc_float * dest2)
{
	char * tmp = new char[len2 - len1 - 3 + 1];
	int cutoff;
	strncpy(tmp, &in1[len1 + 3], len2 - len1 - 3);
	tmp[len2 - len1 - 3] = '\0';
	//every number is saved with 6 digits after the "."
	//find the "."
	cutoff = strchr(tmp, '.') - tmp + 1;
	char * kappa = new char [cutoff + 6 + 1];
	strncpy(kappa, tmp, cutoff + 6);
	kappa[cutoff + 6] = '\0';
	*dest1 = atof(kappa);
	//cut of the part ", mu = "
	strcpy(tmp, &tmp[cutoff + 6 + 7]);
	cutoff = strchr(tmp, '.') - tmp + 1;
	char * mu = new char [cutoff + 6 + 1];
	strncpy(mu, tmp, cutoff + 6);
	mu[cutoff + 6] = '\0';
	*dest2 = atof(mu);
	delete [] tmp;
	delete [] kappa;
	delete [] mu;
	delete [] tmp;
}

void extrInfo_int(const char * in1, int len1, int len2, int * dest)
{
	char * tmp = new char[len2 - len1 - 3 + 1];
	strncpy(tmp, &in1[len1 + 3], len2 - len1 - 3);
	tmp[len2 - len1 - 3] = '\0';
	*dest = (int) atoi(tmp);
	delete [] tmp;
}

// the \n at the end is overwritten by \0
void extrInfo_char(const char * in1, int len1, int len2, char * dest)
{
	char * tmp = new char[len2 - len1 - 4 + 1];
	strncpy(tmp, &in1[len1 + 3], len2 - len1 - 3);
	tmp[len2 - len1 - 4] = '\0';
	strcpy(dest, tmp);
	delete [] tmp;
}

void trim2(char * buff)
{
	int i = 0, j = 0;
	int len = (int)strlen(buff);
	while (i != len) {
		if (buff[i] != '\n' || buff[i] != ' ')
			buff[j++] = buff[i];
		i++;
	}
	buff[j] = 0;
}


void sourcefileparameters::get_XLF_infos(const char * filename, char * hmcversion, char * date )
{
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

// from http://www.codecodex.com/wiki/Remove_blanks_from_a_string#C
void trim(char * buff)
{
	int i = 0, j = 0;
	int len = (int)strlen(buff);
	while (i != len) {
		if (buff[i] != ' ')
			buff[j++] = buff[i];
		i++;
	}
	buff[j] = 0;
}

// http://xmlsoft.org/xmlreader.html
// compile with gcc ReadXML.c $(xml2-config --cflags) -Wall $(xml2-config --libs)

void get_XML_info_simple(xmlTextReaderPtr reader, int numbers[6], char * field)
{
	xmlChar *name, *value;
	name = xmlTextReaderName(reader);
	int type = xmlTextReaderNodeType(reader);
	/*unsigned */
	const char * cmpr[] = {"field", "precision", "flavours", "lx", "ly", "lz", "lt"};
	//check if the desired info follows
	//sometimes there are additional " " that have to be removed with trim(string)
	if (strcmp((char*)name, cmpr[0]) == 0 && type == 1) {
		xmlTextReaderRead(reader);
		value = xmlTextReaderValue(reader);
		trim((char*) value);
		strcpy(field, (char*)value);
		xmlFree(value);
	}
	if (strcmp((char*)name, cmpr[1]) == 0 && type == 1) {
		xmlTextReaderRead(reader);
		value = xmlTextReaderValue(reader);
		numbers[0] = atoi((char*)value);
		xmlFree(value);
	}
	if (strcmp((char*)name, cmpr[2]) == 0 && type == 1) {
		xmlTextReaderRead(reader);
		value = xmlTextReaderValue(reader);
		numbers[1] = atoi((char*)value);
		xmlFree(value);
	}
	if (strcmp((char*)name, cmpr[3]) == 0 && type == 1) {
		xmlTextReaderRead(reader);
		value = xmlTextReaderValue(reader);
		numbers[2] = atoi((char*)value);
		xmlFree(value);
	}
	if (strcmp((char*)name, cmpr[4]) == 0 && type == 1) {
		xmlTextReaderRead(reader);
		value = xmlTextReaderValue(reader);
		numbers[3] = atoi((char*)value);
		xmlFree(value);
	}
	if (strcmp((char*)name, cmpr[5]) == 0 && type == 1) {
		xmlTextReaderRead(reader);
		value = xmlTextReaderValue(reader);
		numbers[4] = atoi((char*)value);
		xmlFree(value);
	}
	if (strcmp((char*)name, cmpr[6]) == 0 && type == 1) {
		xmlTextReaderRead(reader);
		value = xmlTextReaderValue(reader);
		numbers[5] = atoi((char*)value);
		xmlFree(value);
	}
	xmlFree(name);
}

void sourcefileparameters::get_XML_infos(const char * buffer, int size, char * field_out )
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
	strcpy(field_out, field);
}

//todo: merge with xml fcts. above!!
Checksum get_checksum(const char * buffer, int size)
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
      
      strcpy(solvertype_source, solvertype);
      strcpy(hmcversion_solver_source, hmcversion_solver);
      strcpy(date_solver_source, date_solver);
    }
}

void sourcefileparameters::checkLimeEntryForXlfInfos(std::string lime_type, int switcher, LimeReader *r, size_t nbytes)
{
  //!!read XLF info, only FIRST fermion is read!!
  if("xlf-info" == lime_type && switcher < 2) {
    char hmcversion[50];
    char date[50];
    FILE * tmp;

    logger.trace() << "\tfound XLF-infos as lime_type " << lime_type;
    char tmp_file_name[] = "tmpfilenametwo";
    createTemporaryFileToStoreStreamFromLimeReader(tmp_file_name, r, nbytes);
    
    get_XLF_infos(tmp_file_name, hmcversion, date);
    
    remove(tmp_file_name);

    strcpy(hmcversion_source, hmcversion);
    strcpy(date_source, date);
  }
}

void sourcefileparameters::checkLimeEntryForXlmInfos(std::string lime_type, int switcher, LimeReader *r, size_t nbytes)
{
  char field_out[100];
  //!!read ildg format (gauge fields) or etmc-propagator-format (fermions), only FIRST fermion is read!!
  if(("etmc-propagator-format" == lime_type || "ildg-format" == lime_type) && switcher < 2 ) {
    logger.trace() << "\tfound XML-infos as lime_type" << lime_type;
    char * buffer = createBufferAndReadLimeDataIntoIt(r, nbytes);

    get_XML_infos(buffer, nbytes, field_out );
    delete[] buffer;
    logger.trace() << "\tsuccesfully read XMLInfos";
    
    num_entries_source = calcNumberOfEntriesBasedOnFieldType(field_out);
  }
  strcpy(field_source, field_out);
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

int sourcefileparameters::extractInformationFromLimeEntry(LimeReader * r)
{
  int numberOfFermionEntries = 0;
  LimeHeaderData limeHeaderData(r);
  if (limeHeaderData.MB_flag == 1) {
    checkLimeEntry(&numberOfFermionEntries, r, limeHeaderData);
    return 1;
  }
  return 0;
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

int sourcefileparameters::extractBinaryDataFromLimeEntry(LimeReader * r, char ** destination, int * numberOfBinaryDataEntries)
{
    LimeHeaderData limeHeaderData(r);
    if (limeHeaderData.MB_flag == 1) {
      if( strcmp (limeEntryTypes[5], limeHeaderData.limeEntryType.c_str()) == 0 || strcmp (limeEntryTypes[8], limeHeaderData.limeEntryType.c_str()) == 0  )
	{
	  *numberOfBinaryDataEntries= extractBinaryDataFromLimeEntry_NeedsDifferentName(r, limeHeaderData, destination, *numberOfBinaryDataEntries);
	}
      return 1;
    }
    return 0;
}

void sourcefileparameters::goThroughLimeRecordForMetaData(LimeReader * r)
{
  int numberOfLimeEntries = 0;
  int statusOfLimeReader = 0;
  while( (statusOfLimeReader = limeReaderNextRecord(r)) != LIME_EOF ) {
    checkLimeRecordReadForFailure(statusOfLimeReader);
    numberOfLimeEntries += extractInformationFromLimeEntry(r);
  }
  logger.trace() << "Found " << numberOfLimeEntries << " LIME records.";
}

void sourcefileparameters::goThroughLimeRecordForData(LimeReader * r, char ** destination)
{
  int numberOfLimeEntries = 0;
  int statusOfLimeReader = 0;
  int numberOfBinaryDataEntries = 0;
  while( (statusOfLimeReader = limeReaderNextRecord(r)) != LIME_EOF ) {
    checkLimeRecordReadForFailure(statusOfLimeReader);
    numberOfLimeEntries += extractBinaryDataFromLimeEntry(r, destination, &numberOfBinaryDataEntries);
  }
  logger.trace() << "Found " << numberOfLimeEntries << " LIME records.";
  logger.trace() << "Found " << numberOfBinaryDataEntries << " binary entries in LIME file";
}

void sourcefileparameters::readLimeFile(std::string sourceFilename, char ** destination, bool readMetaData)
{
  FILE *limeFileOpenedForReading;
  LimeReader *limeReader;

  limeFileOpenedForReading = fopen (sourceFilename.c_str(), "r");
  limeReader = limeCreateReader(limeFileOpenedForReading);

  if( readMetaData)
    {
      goThroughLimeRecordForMetaData(limeReader);
    }
  else 
    {
      goThroughLimeRecordForData(limeReader, destination);
    }

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
  logger.trace() << "\treading XML-data gave:";
  logger.info() << "\t\tfield type:\t" << field_source ;
  logger.info() << "\t\tprecision:\t" << prec_source ;
  logger.info() << "\t\tlx:\t\t" << lx_source ;
  logger.info() << "\t\tly:\t\t" << ly_source ;
  logger.info() << "\t\tlz:\t\t" << lz_source ;
  logger.info() << "\t\tlt:\t\t" << lt_source ;
  logger.info() << "\t\tflavours:\t" << flavours_source ;
  logger.trace() << "\treading XLF-data gave:";
  logger.info() << "\t\tplaquette:\t" << plaquettevalue_source;
  logger.debug() << "\t\ttrajectorynr:\t" << trajectorynr_source;
  logger.info() << "\t\tbeta:\t\t" << beta_source;
  logger.info() << "\t\tkappa:\t\t" << kappa_source;
  logger.info() << "\t\tmu:\t\t" << mu_source;
  logger.debug() << "\t\tc2_rec:\t\t" << c2_rec_source;
  logger.debug() << "\t\ttime:\t\t" << time_source;
  logger.info() << "\t\thmc-version:\t" << hmcversion_source;
  logger.debug() << "\t\tmubar:\t\t" << mubar_source;
  logger.debug() << "\t\tepsilonbar:\t" << epsilonbar_source;
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

void sourcefileparameters::readsourcefile(std::string sourceFilename, int desiredPrecision, char ** destination)
{
  checkIfFileExists(sourceFilename);

  extractMetadataFromLimeFile(sourceFilename, desiredPrecision);
   
  extractDataFromLimeFile(sourceFilename, destination);
}
