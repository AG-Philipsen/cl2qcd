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

#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define DOUBLE_ENABLED
#else /* cl_khr_fp64 */
#ifdef cl_amd_fp64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#define DOUBLE_ENABLED
#endif /* cl_amd_fp64 */
#endif /* cl_khr_fp64 */

#include "types.h"
#include "types_fermions.h"

//typedef struct { hmc_float re; hmc_float im; } __attribute__((aligned (16))) alignedDpComplex;

#define alignedDpComplex hmc_complex

alignedDpComplex make_alignedDpComplex(const hmc_float re, const hmc_float im)
{
	return (alignedDpComplex) {
		re, im
	};
}

//typedef struct {
//	alignedDpComplex e0;
//	alignedDpComplex e1;
//	alignedDpComplex e2;
//} __attribute((aligned(16))) aligned16DpSu3vec;

#define aligned16DpSu3vec su3vec

aligned16DpSu3vec make_aligned16DpSu3vec(const alignedDpComplex e0, const alignedDpComplex e1, const alignedDpComplex e2)
{
	return (aligned16DpSu3vec) {
		e0, e1, e2
	};
}

//typedef struct {
//	aligned16DpSu3vec e0;
//	aligned16DpSu3vec e1;
//	aligned16DpSu3vec e2;
//	aligned16DpSu3vec e3;
//} dpSpinor;

#define dpSpinor spinor

dpSpinor make_dpSpinor(const aligned16DpSu3vec e0, const aligned16DpSu3vec e1, const aligned16DpSu3vec e2, const aligned16DpSu3vec e3)
{
	return (dpSpinor) {
		e0, e1, e2, e3
	};
}

dpSpinor getDpSpinorFullestSOA(__global const hmc_float * const restrict in, const size_t i, const size_t stride)
{
	return make_dpSpinor(make_aligned16DpSu3vec(make_alignedDpComplex(in[ 0 * stride + i], in[ 1 * stride + i]),
	                     make_alignedDpComplex(in[ 2 * stride + i], in[ 3 * stride + i]),
	                     make_alignedDpComplex(in[ 4 * stride + i], in[ 5 * stride + i])),
	                     make_aligned16DpSu3vec(make_alignedDpComplex(in[ 6 * stride + i], in[ 7 * stride + i]),
	                         make_alignedDpComplex(in[ 8 * stride + i], in[ 9 * stride + i]),
	                         make_alignedDpComplex(in[10 * stride + i], in[11 * stride + i])),
	                     make_aligned16DpSu3vec(make_alignedDpComplex(in[12 * stride + i], in[13 * stride + i]),
	                         make_alignedDpComplex(in[14 * stride + i], in[15 * stride + i]),
	                         make_alignedDpComplex(in[16 * stride + i], in[17 * stride + i])),
	                     make_aligned16DpSu3vec(make_alignedDpComplex(in[18 * stride + i], in[19 * stride + i]),
	                         make_alignedDpComplex(in[20 * stride + i], in[21 * stride + i]),
	                         make_alignedDpComplex(in[22 * stride + i], in[23 * stride + i])));
}
void putDpSpinorFullestSOA(__global hmc_float * const restrict out, const size_t i, const dpSpinor val, const size_t stride)
{
	out[ 0 * stride + i] = val.e0.e0.re;
	out[ 1 * stride + i] = val.e0.e0.im;
	out[ 2 * stride + i] = val.e0.e1.re;
	out[ 3 * stride + i] = val.e0.e1.im;
	out[ 4 * stride + i] = val.e0.e2.re;
	out[ 5 * stride + i] = val.e0.e2.im;
	out[ 6 * stride + i] = val.e1.e0.re;
	out[ 7 * stride + i] = val.e1.e0.im;
	out[ 8 * stride + i] = val.e1.e1.re;
	out[ 9 * stride + i] = val.e1.e1.im;
	out[10 * stride + i] = val.e1.e2.re;
	out[11 * stride + i] = val.e1.e2.im;
	out[12 * stride + i] = val.e2.e0.re;
	out[13 * stride + i] = val.e2.e0.im;
	out[14 * stride + i] = val.e2.e1.re;
	out[15 * stride + i] = val.e2.e1.im;
	out[16 * stride + i] = val.e2.e2.re;
	out[17 * stride + i] = val.e2.e2.im;
	out[18 * stride + i] = val.e3.e0.re;
	out[19 * stride + i] = val.e3.e0.im;
	out[20 * stride + i] = val.e3.e1.re;
	out[21 * stride + i] = val.e3.e1.im;
	out[22 * stride + i] = val.e3.e2.re;
	out[23 * stride + i] = val.e3.e2.im;
}

__kernel void copyDpSpinorFullestSOARestricted(__global hmc_float * const restrict out, __global const hmc_float * const restrict in, const ulong elems, const ulong dummy)
{
	for(size_t i = get_global_id(0); i < elems; i += get_global_size(0)) {
		dpSpinor tmp = getDpSpinorFullestSOA(in, i, elems);
		putDpSpinorFullestSOA(out, i, tmp, elems);
	}
}
