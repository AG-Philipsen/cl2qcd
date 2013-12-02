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

// Description of variables of saxsbypz:
//  - x: The first input staggered field (an su3vec per each site => vector of VOL4D
//       components that are su3vec varibles)
//  - y: The second input staggered field (an su3vec per each site => vector of VOL4D
//       components that are su3vec varibles)
//  - z: The third input staggered field (an su3vec per each site => vector of VOL4D
//       components that are su3vec varibles)
//  - alpha: The complex number by which x has to be multiplied
//  - beta:  The complex number by which y has to be multiplied
//  - out: The output staggered field: alpha*x+beta*y+z (site by site)

__kernel void saxpbypz_staggered(__global const su3vec * const restrict x, __global const su3vec * const restrict y, __global const su3vec * const restrict z, __global const hmc_complex * const restrict alpha, __global const hmc_complex * const restrict beta, __global su3vec * const restrict out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	hmc_complex alpha_tmp = complexLoadHack(alpha);
	hmc_complex beta_tmp = complexLoadHack(beta);
	for(int id_mem = id; id_mem < SPINORFIELDSIZE_MEM; id_mem += global_size) {
		su3vec x_tmp = x[id_mem];
		x_tmp = su3vec_times_complex(x_tmp, alpha_tmp);
		su3vec y_tmp = y[id_mem];
		y_tmp = su3vec_times_complex(y_tmp, beta_tmp);
		su3vec z_tmp = z[id_mem];

		out[id_mem] = su3vec_acc_acc(x_tmp, y_tmp, z_tmp);
	}

	return;
}
