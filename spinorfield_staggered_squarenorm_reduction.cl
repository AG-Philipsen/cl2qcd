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

// Description of variables of global_squarenorm_reduction kernel:
//  - dest: The final result, namely the sum of the squarenorms of all staggered fields.
//  - result_tmp: This is the vector result filled by the kernel global_squarenorm_staggered.
//  - elems: It is the number of components of the vector result_tmp, namely the variable "num_groups"
//            of the kernel global_squarenorm_staggered.

__kernel void global_squarenorm_reduction(__global hmc_float* dest, __global hmc_float* result_tmp, const uint elems)
{
	uint id = get_global_id(0);
	hmc_float tmp = 0;
	if(id == 0) {
		for (uint i = 0; i < elems; i++) {
			tmp += result_tmp[i];
		}
		*dest = tmp;
	}
}