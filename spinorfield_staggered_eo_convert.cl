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
__kernel void convert_to_eoprec_staggered(__global su3vec * const restrict even, __global su3vec * const restrict odd, __global su3vec const * const restrict in)
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

__kernel void convert_from_eoprec_staggered(__global su3vec const * const restrict even, __global su3vec const * const restrict odd, __global su3vec * const restrict out)
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

__kernel void convert_staggered_field_to_SoA_eo(__global staggeredStorageType * const restrict out, __global const su3vec * const restrict in)
{
	for(uint i = get_global_id(0); i < EOPREC_SPINORFIELDSIZE_MEM; i += get_global_size(0)) {
		put_su3vec_to_field_eo(out, i, in[i]);
	}
}

__kernel void convert_staggered_field_from_SoA_eo(__global su3vec * const restrict out, __global const staggeredStorageType * const restrict in)
{
	for(uint i = get_global_id(0); i < EOPREC_SPINORFIELDSIZE_MEM; i += get_global_size(0)) {
		out[i] = get_su3vec_from_field_eo(in, i);
	}
}
