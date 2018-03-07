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

__kernel void set_zero_spinorfield_stagg( __global su3vec * const restrict x )
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);

	for(int id_mem = id; id_mem < SPINORFIELDSIZE_MEM; id_mem += global_size) {
		x[id_mem] = set_su3vec_zero();
	}
	return;
}
