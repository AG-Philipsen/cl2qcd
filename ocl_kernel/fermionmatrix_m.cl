/*
 * Copyright (c) 2011,2012 Christopher Pinke
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

/**
 * @file M normal Wilson fermionmatrix
 */
void dslash_for_site(__global const spinor* const restrict in, __global spinor* const restrict out,
                     __global const Matrixsu3StorageType* const restrict field, hmc_float kappa_in, st_idx const pos)
{
    spinor out_tmp = set_spinor_zero();
    spinor out_tmp2;

    // Diagonalpart: (this is simple here)
    out_tmp = getSpinor(in, get_site_idx(pos));

    // calc dslash (this includes mutliplication with kappa)
    out_tmp2 = dslash_unified_local(in, field, pos, TDIR, kappa_in);
    out_tmp  = spinor_dim(out_tmp, out_tmp2);
    out_tmp2 = dslash_unified_local(in, field, pos, XDIR, kappa_in);
    out_tmp  = spinor_dim(out_tmp, out_tmp2);
    out_tmp2 = dslash_unified_local(in, field, pos, YDIR, kappa_in);
    out_tmp  = spinor_dim(out_tmp, out_tmp2);
    out_tmp2 = dslash_unified_local(in, field, pos, ZDIR, kappa_in);
    out_tmp  = spinor_dim(out_tmp, out_tmp2);

    putSpinor(out, get_site_idx(pos), out_tmp);
}

__kernel void M_wilson(__global const spinor* const restrict in,
                       __global const Matrixsu3StorageType* const restrict field, __global spinor* const restrict out,
                       hmc_float kappa_in)
{
    PARALLEL_FOR (id_local, SPINORFIELDSIZE_LOCAL) {
        // st_idx pos = (evenodd == ODD) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local);
        st_idx pos = get_st_idx_from_site_idx(id_local);
        dslash_for_site(in, out, field, kappa_in, pos);
    }
}
