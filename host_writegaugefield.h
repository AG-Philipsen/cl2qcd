/** @file
 * Storing of gaugefield to files
 *
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

#ifndef _WRITEGAUGEFIELDH_
#define _WRITEGAUGEFIELDH_

#include "common_header_files/types.h"
#include "meta/exceptions.h"

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
