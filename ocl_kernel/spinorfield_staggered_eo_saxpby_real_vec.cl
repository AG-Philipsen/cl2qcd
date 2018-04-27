/*
 * Copyright (c) 2014,2018 Alessandro Sciarra
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

// Description of variables of saxpby:
//  - x: The first input staggered field (an su3vec per each site => vector of VOL4D/2
//       components that are su3vec varibles)
//  - y: The second input staggered field (an su3vec per each site => vector of VOL4D/2
//       components that are su3vec varibles)
//  - alpha: Vector of constants
//  - beta: Vector of constants
//  - index_alpha: Constant alpha to be used
//  - index_beta: Constant beta to be used
//  - out: The output staggered field: alpha*x+beta*y (site by site)

__kernel void saxpby_real_vec_staggered_eoprec(__global const staggeredStorageType* const x,
                                               __global const staggeredStorageType* const y,
                                               __global const hmc_float* const alpha, __global hmc_float* beta,
                                               const int index_alpha, const int index_beta,
                                               __global staggeredStorageType* const out)
{
    const int id          = get_global_id(0);
    const int global_size = get_global_size(0);

    for (int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
        su3vec x_tmp = get_su3vec_from_field_eo(x, id_mem);
        x_tmp        = su3vec_times_real(x_tmp, alpha[index_alpha]);
        su3vec y_tmp = get_su3vec_from_field_eo(y, id_mem);
        y_tmp        = su3vec_times_real(y_tmp, beta[index_beta]);

        su3vec out_tmp = su3vec_acc(y_tmp, x_tmp);
        put_su3vec_to_field_eo(out, id_mem, out_tmp);
    }
}
