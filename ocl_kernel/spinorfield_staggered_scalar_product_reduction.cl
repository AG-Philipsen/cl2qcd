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

// The variable num_values could not be passed to this function since it is
// equal to get_num_groups(0). Nevertheless, since in get_work_size ls, gs
// and num_groups are checked to be corrected and since it is useful
// to have num_groups variable to allocate buffers, we decided to use here
// the num_values variable.
__kernel void scalar_product_reduction(__global hmc_complex* result_tmp, __global hmc_complex* result, const uint num_values)
{
	//!!CP: complex_acc cannot handle __global
	int id = get_global_id(0);
	if(id == 0) {
		hmc_complex tmp1;
		hmc_complex tmp2;
		tmp2 = complexLoadHack(&result_tmp[0]);
		for (int i = 1; i < num_values; i++) {
			tmp1 = complexLoadHack(&result_tmp[i]);
			tmp2 = complexadd(tmp2, tmp1);
		}
		(*result) = tmp2;
	}
}

__kernel void scalar_product_real_reduction(__global hmc_float* result_tmp, __global hmc_float* result, const uint num_values)
{
	int id = get_global_id(0);
	if(id == 0) {
		hmc_float tmp = result_tmp[0];
		for (int i = 1; i < num_values; i++) {
			tmp += result_tmp[i];
		}
		(*result) = tmp;
	}
}
