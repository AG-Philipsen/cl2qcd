/*
 * Copyright (c) 2011-2013 Matthias Bach
 * Copyright (c) 2011 Christopher Pinke
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

// out = alpha*x
__kernel void sax(__global const spinor* const restrict x, __global const hmc_complex* const restrict alpha,
                  __global spinor* const restrict out)
{
    int id          = get_global_id(0);
    int global_size = get_global_size(0);

    hmc_complex alpha_tmp = complexLoadHack(alpha);
    for (int id_mem = id; id_mem < SPINORFIELDSIZE_MEM; id_mem += global_size) {
        spinor x_tmp = x[id_mem];
        out[id_mem]  = spinor_times_complex(x_tmp, alpha_tmp);
    }
}
