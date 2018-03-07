/*
 * Copyright (c) 2013,2014,2018 Alessandro Sciarra
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

__kernel void set_cold_spinorfield_stagg_eoprec(__global staggeredStorageType * const restrict out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	//NOTE: the normalization is the same as in the non eo case, putting |phi|^2 of the WHOLE field to 1!!
	hmc_float norm = 1. / sqrt((hmc_float)(3. * VOL4D_GLOBAL));
	for(int n = id; n < EOPREC_SPINORFIELDSIZE_MEM; n += global_size) {
		put_su3vec_to_field_eo(out, n, su3vec_times_real(set_su3vec_cold(), norm));
	}
}
