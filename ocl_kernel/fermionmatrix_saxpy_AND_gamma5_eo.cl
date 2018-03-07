/*
 * Copyright (c) 2014 Francesca Cuteri
 * Copyright (c) 2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

// -alpha*x + y
//CP: defined with a minus!!!
__kernel void saxpy_AND_gamma5_eo(__global const spinorStorageType * const x, __global const spinorStorageType * const y, const hmc_float alpha_re, const hmc_float alpha_im, __global spinorStorageType * const out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);
	const hmc_complex alpha = (hmc_complex) {
		alpha_re, alpha_im
	};
	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		spinor x_tmp = getSpinor_eo(x, id_mem);
		spinor y_tmp = getSpinor_eo(y, id_mem);
		x_tmp = spinor_times_complex(x_tmp, alpha);
		x_tmp = spinor_dim(y_tmp, x_tmp);
		x_tmp = gamma5_local(x_tmp);
		putSpinor_eo(out, id_mem, x_tmp);
	}
}
