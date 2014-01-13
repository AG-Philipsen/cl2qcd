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

//eoprec operations
__kernel void convert_to_eoprec(__global spinor * const restrict even, __global spinor * const restrict odd, __global spinor const * const restrict in)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	for(int n = id; n < EOPREC_SPINORFIELDSIZE_MEM; n += global_size) {
		st_index pos = get_even_site(n);
		even[n] = in[get_site_idx(pos)];
		pos = get_odd_site(n);
		odd[n] = in[get_site_idx(pos)];
	}
	return;
}

__kernel void convert_from_eoprec(__global spinor const * const restrict even, __global spinor const * const restrict odd, __global spinor * const restrict out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	for(int n = id; n < EOPREC_SPINORFIELDSIZE_MEM; n += global_size) {
		st_index pos = get_even_site(n);
		out[get_site_idx(pos)] = even[n];
		pos = get_odd_site(n);
		out[get_site_idx(pos)] = odd[n];
	}
	return;
}

__kernel void convertSpinorfieldToSOA_eo(__global spinorStorageType * const restrict out, __global const spinor * const restrict in)
{
	for(uint i = get_global_id(0); i < EOPREC_SPINORFIELDSIZE_MEM; i += get_global_size(0)) {
		putSpinor_eo(out, i, in[i]);
	}
}

__kernel void convertSpinorfieldFromSOA_eo(__global spinor * const restrict out, __global const spinorStorageType * const restrict in)
{
	for(uint i = get_global_id(0); i < EOPREC_SPINORFIELDSIZE_MEM; i += get_global_size(0)) {
		out[i] = getSpinor_eo(in, i);
	}
}
