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
 * @file M_tm_minus this is the "total" twisted-mass fermionmatrix (no evenodd) for the lower flavor
 */
__kernel void M_tm_minus(__global const spinor* const restrict in,
                         __global const Matrixsu3StorageType* const restrict field, __global spinor* const restrict out,
                         hmc_float kappa_in, hmc_float mubar_in)
{
    spinor out_tmp, out_tmp2, plus;
    hmc_complex twistfactor       = {1., mubar_in};
    hmc_complex twistfactor_minus = {1., -1. * mubar_in};

    PARALLEL_FOR (id_local, SPINORFIELDSIZE_LOCAL) {
        st_idx pos = get_st_idx_from_site_idx(id_local);
        // get input spinor
        plus = getSpinor(in, get_site_idx(pos));

        // Diagonalpart: this is normal tm-diagonal matrix with negative imaginary part
        out_tmp = M_diag_tm_local(plus, twistfactor_minus, twistfactor);

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
}
