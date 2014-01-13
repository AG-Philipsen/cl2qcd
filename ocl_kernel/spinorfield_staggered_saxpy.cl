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

// Description of variables of saxpy:
//  - x: The first input staggered field (an su3vec per each site => vector of VOL4D
//       components that are su3vec varibles)
//  - y: The second input staggered field (an su3vec per each site => vector of VOL4D
//       components that are su3vec varibles)
//  - alpha: The complex number by which x has to be multiplied
//  - out: The output staggered field: alpha*x+y (site by site)

__kernel void saxpy_staggered(__global su3vec * x, __global su3vec * y, __global const hmc_complex * alpha_p, __global su3vec * out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	const hmc_complex alpha = complexLoadHack(alpha_p);

	for(int id_mem = id; id_mem < SPINORFIELDSIZE_MEM; id_mem += global_size) {
		const su3vec x_tmp = x[id_mem];
		const su3vec y_tmp = y[id_mem];
		const su3vec tmp = su3vec_times_complex(x_tmp, alpha);
		out[id_mem] = su3vec_acc(tmp, y_tmp);
	}
}

//For the moment this kernel is not needed. 
//Uncomment out the region and adapt to staggered fermions (i.e. spinors <---> su3vec) if needed.
/*
// the arguments have been hacked to work on apple
__kernel void saxpy_staggered_arg(__global spinor* x, __global spinor* y, const hmc_float alpha_re, const hmc_float alpha_im, __global spinor* out)
{
	const int id = get_global_id(0);
	const int global_size = get_global_size(0);

	const hmc_complex alpha = (hmc_complex) {
		alpha_re, alpha_im
	};

	for(int id_mem = id; id_mem < SPINORFIELDSIZE_MEM; id_mem += global_size) {
		const spinor x_tmp = x[id_mem];
		const spinor y_tmp = y[id_mem];
		const spinor tmp = spinor_times_complex(x_tmp, alpha);
		out[id_mem] = spinor_dim(y_tmp, tmp);
	}
}
*/