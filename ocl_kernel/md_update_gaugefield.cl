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

/**
 * @file operations used for the molecular dynamics update
 */
/**
 * @file operations used for the molecular dynamics update of the gaugefield
 * u_out = exp(i eps p_in) u_in
 */

__kernel void md_update_gaugefield(const hmc_float eps, __global const aeStorageType * const restrict p_in, __global Matrixsu3StorageType * const restrict u_inout)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	Matrixsu3 tmp;
	Matrixsu3 tmp2;

	//CP: it is GAUGEMOMENTASIZE = NDIM * SPINORFIELDSIZE
	PARALLEL_FOR(index, VOL4D_MEM * NDIM) {
		// an su3 algebra element has NC*NC-1 = 8 hmc_float entries
		// &(p_in[index*8]) should point to the right position for the pos-th element of the long gaugemomentum vector p_in
		tmp2 = build_su3matrix_by_exponentiation(getAe(p_in, index), eps);
		tmp2 = project_su3(tmp2);
		tmp = getSU3(u_inout, index);
		tmp2 = multiply_matrixsu3(tmp2, tmp);
		putSU3(u_inout, index, tmp2);
	}
}
