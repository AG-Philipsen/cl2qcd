#ifndef _WRITEGAUGEFIELDH
#define _WRITEGAUGEFIELDH

#include "hmcerrs.h"
#include "types.h"

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <time.h>

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

hmc_error make_binary_data_single(hmc_float * array, char * out, const int array_size, const int num_entries);

hmc_error make_binary_data_double(hmc_float * array, char * out, const int array_size, const int num_entries);

hmc_error write_gaugefield (
		ildg_gaugefield * array, int array_size, 
		int lx, int ly, int lz, int lt, int prec, int trajectorynr, hmc_float plaquettevalue, hmc_float beta, hmc_float kappa, hmc_float mu, hmc_float c2_rec, hmc_float epsilonbar, hmc_float mubar, 
		const char * hmc_version, const char * filename);

#endif