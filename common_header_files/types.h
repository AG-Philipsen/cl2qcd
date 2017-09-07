/** @file
 * Common types used by all program parts, both on host and OpenCL device.
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

#ifndef _TYPESH_
#define _TYPESH_

/// this definitions rely on "globaldefs.h"

#ifndef _INKERNEL_
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif
#endif

#ifdef _INKERNEL_
#define CONST __constant
#else /* _INKERNEL_ */
#define CONST const
#endif /* _INKERNEL_ */

/** The floating precision type used by hmc, can be 32 or 64 bit. */
#ifdef _USEDOUBLEPREC_
typedef double hmc_float __attribute__((aligned(8)));
#else
typedef float hmc_float __attribute__((aligned(4)));
#endif

/** Complex number type, precision is the same as for hmc_float */
#ifdef _INKERNEL_
typedef struct {
	hmc_float re;
	hmc_float im;
} hmc_complex
#ifdef _USEDOUBLEPREC_
__attribute__((aligned(16)));
#else
__attribute__((aligned(8)));
#endif
typedef struct {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e10;
	hmc_complex e11;
} Matrixsu2;
typedef struct {
	hmc_float e00;
	hmc_float e01;
	hmc_float e10;
	hmc_float e11;
} Matrixsu2_pauli;
#else
struct hmc_complex {
	hmc_float re;
	hmc_float im;
};
struct Matrixsu2 {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e10;
	hmc_complex e11;
} ;
struct Matrixsu2_pauli {
	hmc_float e00;
	hmc_float e01;
	hmc_float e10;
	hmc_float e11;
} ;
#endif

#ifdef _INKERNEL_
__constant hmc_complex hmc_complex_one = {1., 0.};
__constant hmc_complex hmc_complex_zero = {0., 0.};
__constant hmc_complex hmc_complex_minusone = { -1., 0.};
__constant hmc_complex hmc_complex_i = {0., 1.};
__constant hmc_complex hmc_complex_minusi = {0., -1.};
#else
/** A complex 1 */
hmc_complex const hmc_complex_one = {1., 0.};
/** A complex -1 */
hmc_complex const hmc_complex_minusone = { -1., 0.};
/** A complex 0 */
hmc_complex const hmc_complex_zero = {0., 0.};
/** A complex i */
hmc_complex const hmc_complex_i = {0., 1.};
/** A complex -i */
hmc_complex const hmc_complex_minusi = {0., -1.};
#endif

//matrix definitions
#ifdef _INKERNEL_
//a generic 3x3 matrix
typedef struct {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e02;
	hmc_complex e10;
	hmc_complex e11;
	hmc_complex e12;
	hmc_complex e20;
	hmc_complex e21;
	hmc_complex e22;
} Matrix3x3;
//an su3 matrix
typedef struct {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e02;
	hmc_complex e10;
	hmc_complex e11;
	hmc_complex e12;
	hmc_complex e20;
	hmc_complex e21;
	hmc_complex e22;
} Matrixsu3;
#else
//an su3 matrix_
struct Matrixsu3 {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e02;
	hmc_complex e10;
	hmc_complex e11;
	hmc_complex e12;
	hmc_complex e20;
	hmc_complex e21;
	hmc_complex e22;
};
struct Matrix3x3 {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e02;
	hmc_complex e10;
	hmc_complex e11;
	hmc_complex e12;
	hmc_complex e20;
	hmc_complex e21;
	hmc_complex e22;
};
#endif // ifdef _INKERNEL_


#ifdef _USE_SOA_
typedef hmc_complex Matrix3x3StorageType;
#else
typedef Matrix3x3 Matrix3x3StorageType;
#endif

#ifdef _USE_SOA_
typedef hmc_complex Matrixsu3StorageType;
#else
typedef Matrixsu3 Matrixsu3StorageType;
#endif

//CP: this an algebraelement
typedef struct {
	hmc_float e0;
	hmc_float e1;
	hmc_float e2;
	hmc_float e3;
	hmc_float e4;
	hmc_float e5;
	hmc_float e6;
	hmc_float e7;
} ae;

#ifdef _USE_SOA_
typedef hmc_float aeStorageType;
#else
typedef ae aeStorageType;
#endif

#ifndef _INKERNEL_ //Kernels will not take namespaces etc.
namespace common {
	enum startcondition { cold_start = 1, hot_start, start_from_source };
	enum action { wilson = 1, clover, twistedmass, tlsym, iwasaki, dbw2, rooted_stagg };
	enum integrator { leapfrog = 1, twomn };
	enum pbp_version {std = 1, tm_one_end_trick};
	enum solver { cg = 1, bicgstab, bicgstab_save };
	enum sourcetypes {point = 1, volume, timeslice, zslice};
	enum sourcecontents {one = 1, z4, gaussian, z2};
}
#endif

#endif /* _TYPESH_ */
