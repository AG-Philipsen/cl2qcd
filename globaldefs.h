/** @file
 * Global parameters, specified as macros
 *
 * @todo Some, or even most of these, would be better represented as constants
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

#ifndef _GLOBALSH_
#define _GLOBALSH_

/** Number of colors */
#define NC 3
#define NSPIN 4

/** Number of dimensions of the lattice */
#define NDIM 4

//EVEN ODD
#define EVEN 0
#define ODD 1
#define OE 0
#define EO 1

#define ARG_DEF -1

/**
 * PI
 * @todo Rather use PI from stdlib
 */
#define PI  3.14159265358979

#define su2_entries 4

// Definition of numeric constants for the symmetric structure constants d_ijk of su(3) suited for OpenCL
#ifdef _INKERNEL_
/** 1/2 */
#define F_1_2  0.5
/** 1/(2*sqrt(3)) */
#define F_1_2S3 0.288675134594813
/** 1/sqrt(3) */
#define F_1_S3  0.577350269189626
#else
/** 1/2 */
#define F_1_2   (static_cast<hmc_float>(0.5))
/** 1/(2*sqrt(3)) */
#define F_1_2S3 (static_cast<hmc_float>(0.288675134594813))
/** 1/sqrt(3) */
#define F_1_S3  (static_cast<hmc_float>(0.577350269189626))
#endif //_INKERNEL_

#endif


