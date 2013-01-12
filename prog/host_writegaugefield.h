/** @file
 * Storing of gaugefield to files
 */

#ifndef _WRITEGAUGEFIELDH_
#define _WRITEGAUGEFIELDH_

#include "types.h"
#include "exceptions.h"

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

#include "checksum.h"


#define ENDIAN (htons(1) == 1)

/**
 * Write the gaugefield to a file
 *
 * \param array The float array representing the gaugefield.
 * \param array_size The number of floats in the gaugefield array.
 *
 * \todo complete documentation
 */
void write_gaugefield (
  char * binary_data, n_uint64_t num_bytes, Checksum checksum,
  int lx, int ly, int lz, int lt, int prec, int trajectorynr, hmc_float plaquettevalue, hmc_float beta, hmc_float kappa, hmc_float mu, hmc_float c2_rec, hmc_float epsilonbar, hmc_float mubar,
  const char * hmc_version, const char * filename);

#endif /* _WRITEGAUGEFIELDH_ */
