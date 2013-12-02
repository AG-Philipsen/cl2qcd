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

//this calculates the sum of all elements in the matrix
hmc_float sum_up_matrix3x3(const Matrix3x3 in)
{
	return in.e00.re + in.e00.im +
	       in.e01.re + in.e01.im +
	       in.e02.re + in.e02.im +
	       in.e10.re + in.e10.im +
	       in.e11.re + in.e11.im +
	       in.e12.re + in.e12.im +
	       in.e20.re + in.e20.im +
	       in.e21.re + in.e21.im +
	       in.e22.re + in.e22.im ;
}


__kernel void staple_test(__global const Matrixsu3StorageType * const restrict field, __global hmc_float * const restrict out)
{
	PARALLEL_FOR(id, VOL4D_LOCAL) {
		st_index pos = (id % 2 == 0) ? get_even_st_idx_local(id / 2) : get_odd_st_idx_local(id / 2);
		Matrix3x3 V;
		hmc_float res = 0;
		hmc_float res2;
		for(int mu = 0; mu < NDIM; mu++) {
			V = calc_staple(field, pos.space, pos.time, mu);
			res2 = sum_up_matrix3x3(V);
			res += res2;
		}
		int global_pos = get_pos(pos.space, pos.time);
		out[global_pos] = res;
	}

}
