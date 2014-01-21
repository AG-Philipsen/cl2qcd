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
 * Some utility kernels for buffer managment
 */

__kernel void copy_16_bytes(__global float4 * const restrict to, __global const float4 * const restrict from)
{
	// one float4 is exactly 16 bytes, so there is only work for one thread.

	if(get_global_id(0) == 0) {
		to[0] = from[0];
	}
}

#ifdef PARALLEL_FOR
__kernel void clear_bytes(__global char * const restrict dest, const ulong bytes)
{
	PARALLEL_FOR(i, bytes) {
		dest[i] = 0;
	}
}

__kernel void clear_float4(__global float4 * const restrict dest, const ulong elems)
{
	PARALLEL_FOR(i, elems) {
		dest[i] = (float4) {
			0.f, 0.f, 0.f, 0.f
		};
	}
}
#endif
