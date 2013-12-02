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

/** @file
 * Device code implementing complex numbers
 */

//opencl_operations_complex.cl

#include "operations_complex.h"

/** This type can be used to store hmc_complex as it has the same size in bytes */
#ifdef _USEDOUBLEPREC_
typedef float4 hmc_complex_store_type;
#else
typedef float2 hmc_complex_store_type;
#endif

/**
 * Workaround for complex constatns not being loaded properly on APP 2.5
 *
 * For a further discussion of this bug see
 * http://forums.amd.com/devforum/messageview.cfm?catid=390&threadid=156057&enterthread=y
 *
 * @todo make this only be used on APP 2.5
 */
inline hmc_complex complexLoadHack(__global const hmc_complex * p)
{
#ifdef USE_CPLX_LOAD_HACK
	union {
		hmc_complex_store_type v;
		hmc_complex c;
	} tmp;
	tmp.v = *((__global const hmc_complex_store_type*) p);
	return tmp.c;
#else
	return p[0];
#endif
}
