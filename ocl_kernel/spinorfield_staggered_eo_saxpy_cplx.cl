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
//  - x: The first input staggered field (an su3vec per each site => vector of VOL4D/2
//       components that are su3vec varibles)
//  - y: The second input staggered field (an su3vec per each site => vector of VOL4D/2
//       components that are su3vec varibles)
//  - alpha: The complex number by which x has to be multiplied
//  - out: The output staggered field: alpha*x+y (site by site)

__kernel void saxpy_cplx_staggered_eoprec(__global const staggeredStorageType * const x, __global const staggeredStorageType * const y, __global const hmc_complex * const alpha, __global staggeredStorageType * const out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	const hmc_complex alpha_tmp = complexLoadHack(alpha);
	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		su3vec x_tmp = get_su3vec_from_field_eo(x, id_mem);
		su3vec y_tmp = get_su3vec_from_field_eo(y, id_mem);
		x_tmp = su3vec_times_complex(x_tmp, alpha_tmp);
		x_tmp = su3vec_acc(y_tmp, x_tmp);
		put_su3vec_to_field_eo(out, id_mem, x_tmp);
	}
}

__kernel void saxpy_cplx_arg_staggered_eoprec(__global const spinorStorageType * const x, __global const spinorStorageType * const y, const hmc_float alpha_re, const hmc_float alpha_im, __global spinorStorageType * const out)
{
	const int id = get_global_id(0);
	const int global_size = get_global_size(0);

	const hmc_complex alpha = (hmc_complex) {alpha_re, alpha_im};

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		su3vec x_tmp = get_su3vec_from_field_eo(x, id_mem);
		su3vec y_tmp = get_su3vec_from_field_eo(y, id_mem);
		x_tmp = su3vec_times_complex(x_tmp, alpha);
		x_tmp = su3vec_acc(y_tmp, x_tmp);
		put_su3vec_to_field_eo(out, id_mem, x_tmp);
	}
}
