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

__kernel void localQ_test(__global Matrixsu3StorageType * field, __global hmc_float * out)
{
	//CP: this is essentially the plaquette-kernel. The result is however different since the normalization is missing
	PARALLEL_FOR(id_tmp, VOL4D_LOCAL) {
		st_index pos = (id_tmp % 2 == 0) ? get_even_st_idx_local(id_tmp / 2) : get_odd_st_idx_local(id_tmp / 2);
		hmc_float res = 0.;

		//NOTE: The kernel crashes on ATI GPUs if one start with mu=0 here, although this does not make any sense since the nu-loop then does not contribute!!!
		for(int mu = 1; mu < NDIM; mu++) {
			for(int nu = 0; nu < mu; nu++) {
				Matrix3x3 tmp;
				tmp = local_Q_plaquette(field, pos.space, pos.time, mu, nu);
				res += trace_matrix3x3(tmp).re;
			}
		}

		int global_pos = get_pos(pos.space, pos.time);
		out[global_pos] = res;
	}

}
