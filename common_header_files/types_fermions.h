/** @file
 * Common fermion types used by HMC, both on host and OpenCL device.
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

#ifndef _TYPES_FERMIONSH_
#define _TYPES_FERMIONSH_

// relies on "types.h"

typedef struct {
	hmc_complex e0;
	hmc_complex e1;
	hmc_complex e2;
} su3vec
#ifdef _USEDOUBLEPREC_
__attribute__((aligned(16)));
#else
__attribute__((aligned(8)));
#endif

typedef struct {
	su3vec e0;
	su3vec e1;
	su3vec e2;
	su3vec e3;
} spinor
#ifdef _USEDOUBLEPREC_
__attribute__((aligned(32)));
#else
__attribute__((aligned(16)));
#endif

typedef struct {
	su3vec e0;
	su3vec e1;
} halfspinor;

/**
 * The type used for storing spinors on the device.
 */
#ifdef _USE_SOA_
typedef hmc_complex spinorStorageType;
typedef hmc_complex staggeredStorageType;
#else
typedef spinor spinorStorageType;
typedef su3vec staggeredStorageType;
#endif

#endif /* _TYPES_FERMIONSH_ */

