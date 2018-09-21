/*
 * Copyright (c) 2011-2013 Matthias Bach
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

__kernel void gamma5(__global spinor* const restrict inout)
{
    spinor out_tmp;

    PARALLEL_FOR (id_local, SPINORFIELDSIZE_LOCAL) {
        st_idx pos = get_st_idx_from_site_idx(id_local);
        out_tmp    = getSpinor(inout, get_site_idx(pos));
        out_tmp    = gamma5_local(out_tmp);
        putSpinor(inout, get_site_idx(pos), out_tmp);
    }
}
