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

// -alpha*x + y
//CP: defined with a minus!!!

__kernel void saxpy(__global spinor* x, __global spinor* y, __global const hmc_complex * alpha_p, __global spinor* out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	const hmc_complex alpha = complexLoadHack(alpha_p);

	for(int id_mem = id; id_mem < SPINORFIELDSIZE_MEM; id_mem += global_size) {
		const spinor x_tmp = x[id_mem];
		const spinor y_tmp = y[id_mem];
		const spinor tmp = spinor_times_complex(x_tmp, alpha);
		out[id_mem] = spinor_dim(y_tmp, tmp);
	}
}

// the arguments have been hacked to work on apple
__kernel void saxpy_arg(__global spinor* x, __global spinor* y, const hmc_float alpha_re, const hmc_float alpha_im, __global spinor* out)
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
