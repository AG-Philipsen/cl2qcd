#ifndef _READGAUGEH_
#define _READGAUGEH_

#include "hmcerrs.h"
#include "types.h"

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>


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

#define ENDIAN (htons(1) == 1)


void extrInfo_hmc_float(const char * in1, const char * in2, int len1, int len2, hmc_float * dest);

// two strings in xlf-info and inverter-info are complicated because there are several vars saved in them
// this is not a beautiful implementation!!
void extrInfo_beta(const char * in1, const char * in2, int len1, int len2, hmc_float * dest1, hmc_float * dest2, hmc_float * dest3, hmc_float * dest4);

void extrInfo_int(const char * in1, const char * in2, int len1, int len2, int * dest);
// the \n at the end is overwritten by \0
void extrInfo_char(const char * in1, const char * in2, int len1, int len2, char * dest);

hmc_error get_XLF_infos(char * filename, hmc_float * plaquettevalue, int * trajectorynr, hmc_float * beta, hmc_float * kappa, hmc_float * mu, 
		      hmc_float * c2_rec, int * time, char * hmcversion, hmc_float * mubar, hmc_float * epsilonbar, char * date );

hmc_error get_inverter_infos(char * filename, char * solver, hmc_float * epssq, int * noiter, hmc_float * kappa_solver, hmc_float * mu_solver, 
		      int * time, char * hmcversion, char * date );

// http://xmlsoft.org/xmlreader.html
// compile with gcc ReadXML.c $(xml2-config --cflags) -Wall $(xml2-config --libs)

void trim(char * buff);

void get_XML_info_simple(xmlTextReaderPtr reader, int numbers[6], char * field);

hmc_error get_XML_infos(const char * filename, int * prec, int * lx, int * ly, int * lz, int *lt, int * flavours, char * field_out );

// get XML Infos: file to be read + parameters
hmc_error read_meta_data(char * file, int * lx, int * ly, int * lz, int * lt, int * prec, char * field_out, int * num_entries, 
		 int * flavours, hmc_float * plaquettevalue, int * trajectorynr, hmc_float * beta, hmc_float * kappa, hmc_float * mu, hmc_float * c2_rec, int * time, char * hmcversion, hmc_float * mubar, hmc_float * epsilonbar, char * date, 
		 char * solvertype, hmc_float * epssq, int * noiter, hmc_float * kappa_solver, hmc_float * mu_solver,  int * time_solver, char * hmcversion_solver, char * date_solver, int * fermion);
    
hmc_error read_binary_data_single(char * file, float * numArray, int num_entries, int filelength );
    
int read_data_single(char * file, float * num_array_single, int num_entries, char * field_out);
      
int read_binary_data_double(char * file, double * numArray, int num_entries, int filelength );
    
int read_data_double(char * file, double * num_array_double, int num_entries, char * field_out);

hmc_error read_tmlqcd_file(char * file, 
			   int * lx, int * ly, int * lz, int * lt, int * prec, char * field_out, int * num_entries, int * flavours, 
			   hmc_float * plaquettevalue, int * trajectorynr, hmc_float * beta, hmc_float * kappa, hmc_float * mu, hmc_float * c2_rec, int * time, char * hmcversion, hmc_float * mubar, hmc_float * epsilonbar, char * date, 
			   char * solvertype, hmc_float * epssq, int * noiter, hmc_float * kappa_solver, hmc_float * mu_solver, int * time_solver, char * hmcversion_solver, char * date_solver,
			   hmc_float * array, int * hmc_prec);

			   
class sourcefileparameters {
 public:
  sourcefileparameters(){set_defaults();};
  hmc_error readsourcefile(char * file, int precision, hmc_float ** array);
  hmc_error set_defaults();
  void val_assign_source(int * out, int in);
  void val_assign_source(hmc_float * out, hmc_float in);  
  int lx_source, ly_source, lz_source, lt_source, prec_source, num_entries_source, flavours_source, 
      trajectorynr_source, time_source, time_solver_source, noiter_source;
  hmc_float plaquettevalue_source, beta_source, kappa_source, mu_source, c2_rec_source, mubar_source, epsilonbar_source,
      epssq_source, kappa_solver_source, mu_solver_source;
  char field_source[100]; char hmcversion_source[50]; char date_source[50]; char solvertype_source[50];
  char hmcversion_solver_source[50]; char date_solver_source[50];
};

#endif
