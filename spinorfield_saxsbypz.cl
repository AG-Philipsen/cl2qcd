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

//alpha*x + beta*y + z
__kernel void saxsbypz(__global const spinor * const restrict x, __global const spinor * const restrict y, __global const spinor * const restrict z, __global const hmc_complex * const restrict alpha, __global const hmc_complex * const restrict beta, __global spinor * const restrict out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	hmc_complex alpha_tmp = complexLoadHack(alpha);
	hmc_complex beta_tmp = complexLoadHack(beta);
	for(int id_mem = id; id_mem < SPINORFIELDSIZE_MEM; id_mem += global_size) {
		spinor x_tmp = x[id_mem];
		x_tmp = spinor_times_complex(x_tmp, alpha_tmp);
		spinor y_tmp = y[id_mem];
		y_tmp = spinor_times_complex(y_tmp, beta_tmp);
		spinor z_tmp = z[id_mem];

		out[id_mem] = spinor_acc_acc(y_tmp, x_tmp, z_tmp);
	}

	return;
}
