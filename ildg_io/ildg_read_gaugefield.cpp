/*
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


void get_XLF_infos(const char * filename, hmc_float * plaquettevalue, int * trajectorynr, hmc_float * beta, hmc_float * kappa, hmc_float * mu,
                   hmc_float * c2_rec, int * time, char * hmcversion, hmc_float * mubar, hmc_float * epsilonbar, char * date )
{
	FILE * reader;
	reader = fopen(filename, "r");
	if (reader != NULL) {
		//there are " " inf front of most of the labels, the last one was added for a different style
		const char * tmparray [] = {"plaquette", " trajectory nr", " beta", "kappa", "mu", "c2_rec", " time", " hmcversion", " mubar", " epsilonbar", " date", " plaquette"};
		char tmp1[512];

		while ( fgets (tmp1, 512, reader) != NULL) {
			trim2(tmp1);
			if(strncmp(tmparray[0], tmp1, strlen(tmparray[0])) == 0) extrInfo_hmc_float(tmp1, strlen(tmparray[0]), strlen(tmp1), plaquettevalue);
			if(strncmp(tmparray[1], tmp1, strlen(tmparray[1])) == 0) extrInfo_int(tmp1, strlen(tmparray[1]), strlen(tmp1), trajectorynr);
			if(strncmp(tmparray[2], tmp1, strlen(tmparray[2])) == 0) extrInfo_beta(tmp1, strlen(tmparray[2]), strlen(tmp1), beta, kappa, mu, c2_rec);
			if(strncmp(tmparray[6], tmp1, strlen(tmparray[6])) == 0) extrInfo_int(tmp1, strlen(tmparray[6]), strlen(tmp1), time);
			if(strncmp(tmparray[7], tmp1, strlen(tmparray[7])) == 0) extrInfo_char(tmp1, strlen(tmparray[7]), strlen(tmp1), hmcversion);
			if(strncmp(tmparray[8], tmp1, strlen(tmparray[8])) == 0) extrInfo_hmc_float(tmp1, strlen(tmparray[8]), strlen(tmp1), mubar);
			if(strncmp(tmparray[9], tmp1, strlen(tmparray[9])) == 0) extrInfo_hmc_float(tmp1, strlen(tmparray[9]), strlen(tmp1), epsilonbar);
			if(strncmp(tmparray[10], tmp1, strlen(tmparray[10])) == 0) extrInfo_char(tmp1, strlen(tmparray[10]), strlen(tmp1), date);
			if(strncmp(tmparray[11], tmp1, strlen(tmparray[11])) == 0) extrInfo_hmc_float(tmp1,  strlen(tmparray[11]), strlen(tmp1), plaquettevalue);
		}
	} else throw File_Exception(filename);

	return;
}

void get_inverter_infos(const char * filename, char * solver, hmc_float * epssq, int * noiter, hmc_float * kappa_solver, hmc_float * mu_solver,
                        int * time, char * hmcversion, char * date )
{
	FILE * reader;
	reader = fopen(filename, "r");
	if (reader != NULL) {
		//there are " " inf front of most of the labels
		const char * tmparray [] = {"solver", " epssq", " noiter", " kappa", "mu", " time", " hmcversion", " date"};
		char tmp1[512];
		while ( fgets (tmp1, 512, reader) != NULL ) {
			if(strncmp(tmparray[0], tmp1, strlen(tmparray[0])) == 0) extrInfo_char(tmp1, strlen(tmparray[0]), strlen(tmp1), solver);
			if(strncmp(tmparray[1], tmp1, strlen(tmparray[1])) == 0) extrInfo_hmc_float(tmp1,  strlen(tmparray[1]), strlen(tmp1), epssq);
			if(strncmp(tmparray[2], tmp1, strlen(tmparray[2])) == 0) extrInfo_int (tmp1, strlen(tmparray[2]), strlen(tmp1), noiter);
			if(strncmp(tmparray[3], tmp1, strlen(tmparray[3])) == 0) extrInfo_kappa(tmp1, strlen(tmparray[3]), strlen(tmp1), kappa_solver, mu_solver);
			if(strncmp(tmparray[5], tmp1, strlen(tmparray[5])) == 0) extrInfo_int (tmp1, strlen(tmparray[5]), strlen(tmp1), time);
			if(strncmp(tmparray[6], tmp1, strlen(tmparray[6])) == 0) extrInfo_char(tmp1, strlen(tmparray[6]), strlen(tmp1), hmcversion);
			if(strncmp(tmparray[7], tmp1, strlen(tmparray[7])) == 0) extrInfo_char(tmp1, strlen(tmparray[7]), strlen(tmp1), date);
		}
	} else throw File_Exception(filename);

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

void sourcefileparameters::get_XML_infos(const char * buffer, int size, const char * filename, int * prec, int * lx, int * ly, int * lz, int *lt, int * flavours, char * field_out )
{
	xmlTextReaderPtr reader;
	int ret;
	int tmpArray[6];
	char field[100];

	reader = xmlReaderForMemory(buffer, size, filename, nullptr, 0);
	if (reader != NULL) {
		ret = xmlTextReaderRead(reader);
		while (ret == 1) {
			get_XML_info_simple(reader, tmpArray, field);
			ret = xmlTextReaderRead(reader);
		}
		xmlFreeTextReader(reader);
		if (ret != 0) {
			logger.warn() << filename << ": failed to parse";
		}
	} else throw File_Exception(filename);

	*prec = tmpArray[0], *flavours = tmpArray[1];
	*lx = tmpArray[2], *ly = tmpArray[3], *lz = tmpArray[4], *lt = tmpArray[5],
	 strcpy(field_out, field);
	return;
}

Checksum get_checksum(const char * buffer, int size, const char * filename)
{
	uint32_t suma, sumb;

	xmlTextReaderPtr reader = xmlReaderForMemory(buffer, size, filename, nullptr, 0);
	if(reader == nullptr) {
		throw File_Exception(filename);
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
	if(xml_state) {
		logger.warn() << filename << ": failed to parse";
	}
	xmlFreeTextReader(reader);

	return Checksum(suma, sumb);
}

// read in binary file and save it as readable file
// since tmLQCD always saves data with BigEndian one has to be careful

// get XML Infos: file to be read + parameters
void sourcefileparameters::read_meta_data(const char * file, int * lx, int * ly, int * lz, int * lt, int * prec, char * field_out, int * num_entries,
                    int * flavours, hmc_float * plaquettevalue, int * trajectorynr, hmc_float * beta, hmc_float * kappa, hmc_float * mu, hmc_float * c2_rec, int * time, char * hmcversion, hmc_float * mubar, hmc_float * epsilonbar, char * date,
                    char * solvertype, hmc_float * epssq, int * noiter, hmc_float * kappa_solver, hmc_float * mu_solver,  int * time_solver, char * hmcversion_solver, char * date_solver, int * fermion, Checksum * checksum)
{
	FILE *fp;
	int MB_flag, ME_flag, msg, rec, status, first, switcher = 0;
	size_t bytes_pad;
	n_uint64_t nbytes;

	//possible lime_types
	const char * lime_types[] = {
		"propagator-type", "xlf-info", "inverter-info", "gauge-scidac-checksum-copy", "etmc-propagator-format",
		"scidac-binary-data", "scidac-checksum", "ildg-format", "ildg-binary-data"
	};

	//read lime file
	fp = fopen (file, "r");
	LimeReader *r;
	r = limeCreateReader(fp);
	first = 1;
	msg = 0;
	//go through the lime-entries
	while( (status = limeReaderNextRecord(r)) != LIME_EOF ) {
		if( status != LIME_SUCCESS ) {
			char errmsg[256];
			sprintf(errmsg, "\t\tlimeReaderNextRecord returned status = %d\n", status);
			throw Print_Error_Message(errmsg);
		}
		if (MB_flag == 1 || first) {
			first = 0;
			rec = 0;
			msg++;
		}
		rec++;
		//read header data
		nbytes    = limeReaderBytes(r);
		std::string lime_type(limeReaderType(r));
		bytes_pad = limeReaderPadBytes(r);
		MB_flag   = limeReaderMBFlag(r);
		ME_flag   = limeReaderMEFlag(r);
		if("propagator-type" == lime_type) {
			if (switcher == 0) {
				logger.info() << "\tfile contains fermion informations" ;
				switcher ++;
			} else {
				logger.info() << "\tfile contains informations about more than one fermion: " << switcher;
				switcher ++;
			}
		}
		//!!read the inverter-infos for FIRST fermion infos only!!
		if("inverter-info" == lime_type && switcher == 1) {
			logger.trace() << "\tfound inverter-infos as lime_type " << lime_type ;
			*fermion = *fermion + 1;
			//!!create tmporary file to read in data, this can be done better
			FILE * tmp;
			const char tmp_file_name[] = "tmpfilenameone";
			tmp = fopen(tmp_file_name, "w");
			if(tmp == NULL) throw Print_Error_Message("\t\terror in creating tmp file\n");

			char * buffer = new char[nbytes + 1];
			int error = limeReaderReadData(buffer, &nbytes, r);
			if(error != 0) throw Print_Error_Message("Something went wrong...", __FILE__, __LINE__);
			fwrite(buffer, 1, nbytes, tmp);
			fclose(tmp);
			delete [] buffer;
			buffer = 0;

			get_inverter_infos(tmp_file_name, solvertype, epssq, noiter, kappa_solver, mu_solver, time_solver, hmcversion_solver, date_solver);
			logger.trace() << "\tsuccesfully read InverterInfos" ;

			remove(tmp_file_name);
		}
		//!!read XLF info, only FIRST fermion is read!!
		if("xlf-info" == lime_type && switcher < 2) {
			logger.trace() << "\tfound XLF-infos as lime_type " << lime_type;
			//!!create tmporary file to read in data, this can be done better
			FILE * tmp;
			const char tmp_file_name[] = "tmpfilenametwo";
			tmp = fopen(tmp_file_name, "w");
			if(tmp == NULL) {
				throw Print_Error_Message("\t\terror in creating tmp file\n");
			}
			char * buffer = new char[nbytes + 1];
			int error = limeReaderReadData (buffer, &nbytes, r);
			if(error != 0) throw Print_Error_Message("Something went wrong...", __FILE__, __LINE__);
			fwrite(buffer, 1, nbytes, tmp);
			fclose(tmp);
			delete [] buffer;
			buffer = 0;

			get_XLF_infos(tmp_file_name, plaquettevalue, trajectorynr, beta, kappa, mu, c2_rec, time, hmcversion, mubar, epsilonbar, date);
			logger.trace() << "\tsuccesfully read XLFInfos";

			remove(tmp_file_name);
		}
		//!!read ildg format (gauge fields) or etmc-propagator-format (fermions), only FIRST fermion is read!!
		if(("etmc-propagator-format" == lime_type || "ildg-format" == lime_type) && switcher < 2 ) {
			logger.trace() << "\tfound XML-infos as lime_type" << lime_type;
			//!!create tmporary file to read in data, this can be done better
			char * buffer = new char[nbytes + 1];
			int error = limeReaderReadData(buffer, &nbytes, r);
			if(error != 0) throw Print_Error_Message("Something went wrong...", __FILE__, __LINE__);

			get_XML_infos(buffer, nbytes, file, prec, lx, ly, lz, lt, flavours, field_out );
			delete[] buffer;
			buffer = 0;
			logger.trace() << "\tsuccesfully read XMLInfos";

			// different sizes for fermions or gauge fields
			if(strcmp(field_out, "diracFermion") == 0) {
				//latSize sites, 4 dirac indices, Nc colour indices, 2 complex indices
				*num_entries = (int) (*lx) * (*ly) * (*lz) * (*lt) * NC * NSPIN * 2;
			} else if(strcmp(field_out, "su3gauge") == 0) {
				// latSize sites, 4 links, 2 complex indices -> 9 complex numbers per link
				*num_entries = (int) (*lx) * (*ly) * (*lz) * (*lt) * 2 * 4 * 9;
			} else {
				throw Print_Error_Message("\tError in read_meta_infos()");
			}
		}
		if("scidac-checksum" == lime_type) {
			char * buffer = new char[nbytes + 1];
			int error = limeReaderReadData (buffer, &nbytes, r);
			if(error != 0) throw Print_Error_Message("Something went wrong...", __FILE__, __LINE__);
			buffer[nbytes] = 0;
			*checksum = get_checksum(buffer, nbytes, file);
			delete[] buffer;
		}
	}
	limeDestroyReader(r);
	fclose(fp);
	return;
}

//LZ removed last, unused parameter char* field_out
void read_data(const char * file, char * data, size_t bytes)
{
	FILE *fp;
	int MB_flag, ME_flag, msg, rec, status, first, cter = 0;
	char *lime_type;
	size_t bytes_pad;
	n_uint64_t nbytes;

	//possible lime_types
	const char * lime_types[] = {
		"propagator-type", "xlf-info", " inverter-info", "gauge-scidac-checksum-copy", "etmc-propagator-format",
		"scidac-binary-data", "scidac-checksum", "ildg-format", "ildg-binary-data"
	};

	//read lime file
	fp = fopen (file, "r");
	LimeReader *r;
	r = limeCreateReader(fp);
	first = 1;
	msg = 0;
	while( (status = limeReaderNextRecord(r)) != LIME_EOF ) {
		if( status != LIME_SUCCESS ) {
			char errmsg[256];
			sprintf(errmsg, "\tlimeReaderNextRecord returned status = %d\n",
			        status);
			throw Print_Error_Message(errmsg);
		}
		if (MB_flag == 1 || first) {
			first = 0;
			rec = 0;
			msg++;
		}
		rec++;
		//read header data
		nbytes    = limeReaderBytes(r);
		lime_type = limeReaderType(r);
		bytes_pad = limeReaderPadBytes(r);
		MB_flag   = limeReaderMBFlag(r);
		ME_flag   = limeReaderMEFlag(r);
		//!! read data only for the FIRST entry!!
		if( (strcmp (lime_types[5], lime_type) == 0 || strcmp (lime_types[8], lime_type) == 0 ) && cter < 1) {
			if(nbytes != bytes) {
				throw Invalid_Parameters("Binary data does not have expected size.", bytes, nbytes);
			}
			cter ++;
			limeReaderReadData(data, &nbytes, r);
		}
	}

	limeDestroyReader(r);
	fclose(fp);
}

void sourcefileparameters::read_tmlqcd_file(const char * file, char ** array, int * hmc_prec, Checksum * checksum)
{
	int lx, ly, lz, lt, prec, num_entries, flavours, trajectorynr, time, time_solver, noiter;
	hmc_float plaquettevalue, beta, kappa, mu, c2_rec, mubar, epsilonbar, epssq, kappa_solver, mu_solver;
	char field_out[100];
	char hmcversion[50];
	char date[50];
	char solvertype[50];
	char hmcversion_solver[50];
	char date_solver[50];

	int fermion = 0;
	FILE * checker;
	checker = fopen(file, "r");
	if(checker == 0) {
	  throw File_Exception(file);
	} else {
		logger.info() << " " ;
		logger.info() << "*************************************************************" ;
		logger.info() << "*************************************************************" ;
		logger.info() << "Reading gaugefield configuration from file " << file << "...";
	}
	fclose(checker);

	logger.info() << "\tMetadata:" ;
	read_meta_data(file, &lx, &ly, &lz, &lt, &prec, field_out, &num_entries, &flavours, &plaquettevalue, &trajectorynr,
	               &beta, &kappa, &mu, &c2_rec, &time, hmcversion, &mubar, &epsilonbar, date,
	               solvertype, &epssq, &noiter, &kappa_solver, &mu_solver, &time_solver, hmcversion_solver, date_solver, &fermion, checksum);

	lx_source = lx;
	ly_source = ly;
	lz_source = lz;
	lt_source = lt;
	prec_source = prec;
	strcpy(field_source, field_out);
	num_entries_source = num_entries;
	flavours_source = flavours;
	plaquettevalue_source = plaquettevalue;
	trajectorynr_source = trajectorynr;
	beta_source = beta;
	kappa_source = kappa;
	mu_source = mu;
	c2_rec_source = c2_rec;
	time_source = time;
	strcpy(hmcversion_source, hmcversion);
	mubar_source = mubar;
	epsilonbar_source = epsilonbar;
	strcpy(date_source, date);
	strcpy(solvertype_source, solvertype);
	epssq_source = epssq;
	noiter_source = noiter;
	kappa_solver_source = kappa_solver;
	mu_solver_source = mu_solver;
	time_solver_source = time_solver;
	strcpy(hmcversion_solver_source, hmcversion_solver);
	strcpy(date_solver_source, date_solver);

	
	logger.trace() << "\treading XML-file gave:";
	logger.info() << "\t\tfield type:\t" << field_out ;
	logger.info() << "\t\tprecision:\t" << prec ;
	logger.info() << "\t\tlx:\t\t" << lx_source ;
	logger.info() << "\t\tlx:\t\t" << lx ;
	logger.info() << "\t\tly:\t\t" << ly ;
	logger.info() << "\t\tlz:\t\t" << lz ;
	logger.info() << "\t\tlt:\t\t" << lt ;
	logger.debug() << "\t\tflavours:\t" << flavours ;
	logger.trace() << "\treading XLF-data gave:";
	logger.info() << "\t\tplaquette:\t" << plaquettevalue;
	logger.debug() << "\t\ttrajectorynr:\t" << trajectorynr;
	logger.info() << "\t\tbeta:\t\t" << beta;
	logger.info() << "\t\tkappa:\t\t" << kappa;
	logger.info() << "\t\tmu:\t\t" << mu;
	logger.debug() << "\t\tc2_rec:\t\t" << c2_rec;
	logger.debug() << "\t\ttime:\t\t" << time;
	logger.info() << "\t\thmc-version:\t" << hmcversion;
	logger.debug() << "\t\tmubar:\t\t" << mubar;
	logger.debug() << "\t\tepsilonbar:\t" << epsilonbar;
	logger.info() << "\t\tdate:\t\t" << date;
	if(fermion != 0) {
		logger.info() << "\treading inverter-data gave:";
		logger.info() << "\t\tsolvertype:\t" << solvertype;
		logger.info() << "\t\tepssq:\t\t" << std::setprecision(30) << epssq;
		logger.info() << "\t\tnoiter:\t\t" << noiter;
		logger.info() << "\t\tkappa_solver:\t" << kappa_solver;
		logger.info() << "\t\tmu_solver:\t" << mu_solver;
		logger.info() << "\t\ttime_solver:\t" << time_solver;
		logger.info() << "\t\thmc-ver_solver:\t" << hmcversion_solver;
		logger.info() << "\t\tdate_solver:\t" << date_solver;
	}
	logger.debug() << "\tfile-checksum:\t" << *checksum;
	
	if(*hmc_prec != prec) {
		throw Print_Error_Message("\nthe precision of hmc and sourcefile do not match, will not read data!!!", __FILE__, __LINE__);
	} else {
		logger.trace() << "reading data..";
		//!!note: the read-routines were not changed, the array is just set to the values of the num_array`s
		size_t datasize = num_entries * sizeof(hmc_float);
		*array = new char[datasize];
		read_data(file, *array, datasize);
		logger.trace() << "\tsuccesfully read in data";
	}
	logger.trace() << "\nsuccesfully read tmlqcd-file " << file;
	logger.info() << "*************************************************************" ;
	logger.info() << " ";
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


void sourcefileparameters::readsourcefile(const char * file, int precision, char ** array)
{

	int  prec_tmp;

	//CP
	//this was done because i am lazy
	prec_tmp = precision;
	read_tmlqcd_file( file, array, &prec_tmp, &checksum);

	return;
}

