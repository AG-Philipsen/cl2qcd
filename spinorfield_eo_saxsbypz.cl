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
__kernel void saxsbypz_eoprec(__global const spinorStorageType * const x, __global const spinorStorageType * const y, __global const spinorStorageType * const z, __global const hmc_complex * const alpha, __global hmc_complex * beta, __global spinorStorageType * const out)
{
	const int id = get_global_id(0);
	const int global_size = get_global_size(0);

	const hmc_complex alpha_tmp = complexLoadHack(alpha);
	const hmc_complex beta_tmp = complexLoadHack(beta);
	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		spinor x_tmp = getSpinor_eo(x, id_mem);
		x_tmp = spinor_times_complex(x_tmp, alpha_tmp);
		spinor y_tmp = getSpinor_eo(y, id_mem);
		y_tmp = spinor_times_complex(y_tmp, beta_tmp);
		spinor z_tmp = getSpinor_eo(z, id_mem);

		spinor out_tmp = spinor_acc_acc(y_tmp, x_tmp, z_tmp);
		putSpinor_eo(out, id_mem, out_tmp);
	}

	return;
}
