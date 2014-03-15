#ifndef _PARSER_UTILS_H_
#define _PARSER_UTILS_H_

/**
 * Functions used to extract informations from lime entries stored in buffers.
 */ 

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


#endif
