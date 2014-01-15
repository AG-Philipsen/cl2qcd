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
 * Definitions required in the Device code.
 */

//opencl_header.cl

#ifdef ENABLE_PRINTF
// we def out anything using printf if unsupported,
// therefore we also don't need to enable printf
#pragma OPENCL EXTENSION cl_amd_printf : enable
#endif

#ifdef _USEDOUBLEPREC_
#ifdef _DEVICE_DOUBLE_EXTENSION_AMD_
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif
#ifdef _DEVICE_DOUBLE_EXTENSION_KHR_
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif
#endif

/// this requires "globaldefs.h" and  "types.h"
#ifdef _USE_BLOCKED_LOOPS_
#define PARALLEL_FOR(VAR, LIMIT) \
size_t _block_size = (LIMIT + get_global_size(0) - 1) / get_global_size(0); \
for(size_t VAR = get_global_id(0) * _block_size; VAR < (get_global_id(0) + 1) * _block_size && VAR < LIMIT; ++VAR)
#else /* _USE_BLOCKED_LOOPS_ */
#if (NSPACE / 16 ) * 16 == NSPACE /* On Cypress with APP 2.6 we need different strategies for NSPACE % 16 == 0 and NSPACE % 16 == 0 */
#define PARALLEL_FOR(VAR, LIMIT) \
size_t _global_size = get_global_size(0); /* required to avoid performance regression an Cypress with APP 2.6 */ \
for(size_t VAR = get_global_id(0); VAR < LIMIT; VAR += _global_size)
#else /* NSPACE % 16 == 0 */
#define PARALLEL_FOR(VAR, LIMIT) \
for(size_t VAR = get_global_id(0); VAR < LIMIT; VAR += get_global_size(0))
#endif /* NSPACE % 16 == 0 */
#endif /* _USE_BLOCKED_LOOPS_ */
