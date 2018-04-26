/*
 * Copyright (c) 2013,2018 Alessandro Sciarra
 * Copyright (c) 2013 Matthias Bach
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
