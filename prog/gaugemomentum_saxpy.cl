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
 * @file
 * Kernel for the gaugemomentum saxpy operation.
 *  @param x The first input gaugemomentum field (an ae per each site)
 *  @param y The second input gaugemomentum field (an su3vec per each site)
 *  @param alpha The REAL number which x has to be multiplied by
 *  @param out The output gaugemomentum field: alpha*x+y (site by site)
 */
__kernel void gaugemomentum_saxpy(__global const aeStorageType * const x, __global const aeStorageType * const y, __global const hmc_float * const alpha, __global aeStorageType * const out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	for(int id_mem = id; id_mem < GAUGEMOMENTASIZE_MEM; id_mem += global_size) {
		ae x_tmp = getAe(x, id_mem);
		ae y_tmp = getAe(y, id_mem);
		putAe(out, id_mem, acc_factor_times_algebraelement(y_tmp, *alpha, x_tmp));
	}
}



#if 0
/**
 * @file operations used for the molecular dynamics update of the gauge momenta
 * p_out = p_in - eps/2 force(u_in, phi)
 * It is assumed that the force term has already been computed. Then one only has real-vectors and this is essentially adding one vector to another...
 */
__kernel void md_update_gaugemomenta(const hmc_float eps, __global aeStorageType * const restrict p_inout, __global const aeStorageType * const restrict force_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	for(int id_mem = id; id_mem < GAUGEMOMENTASIZE_MEM; id_mem += global_size) {
		update_gaugemomentum(getAe(force_in, id_mem), eps, id_mem, p_inout);
	}
}
#endif