/*
 * Copyright (c) 2011,2012 Christopher Pinke
 * Copyright (c) 2011-2013 Matthias Bach
 * Copyright (c) 2013,2018 Alessandro Sciarra
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
 @file fermionmatrix-functions that are called from within kernels
*/

// local twisted-mass Diagonalmatrix:
//    (1+i*mubar*gamma_5)psi = (1, mubar)psi.0,1 (1,-mubar)psi.2,3
spinor inline M_diag_tm_local(spinor in, hmc_complex factor1, hmc_complex factor2)
{
    spinor tmp;
    tmp.e0 = su3vec_times_complex(in.e0, factor1);
    tmp.e1 = su3vec_times_complex(in.e1, factor1);
    tmp.e2 = su3vec_times_complex(in.e2, factor2);
    tmp.e3 = su3vec_times_complex(in.e3, factor2);
    return tmp;
}

/** @todo this can be optimized... */
// local gamma5:
//    (gamma_5)psi = (1)psi.0,1 (-1)psi.2,3
spinor inline gamma5_local(spinor in)
{
    spinor tmp;
    tmp.e0 = in.e0;
    tmp.e1 = in.e1;
    tmp.e2 = su3vec_times_real(in.e2, -1.);
    tmp.e3 = su3vec_times_real(in.e3, -1.);
    return tmp;
}

spinor dslash_unified_local(__global const spinor* const restrict in,
                            __global Matrixsu3StorageType const* const restrict field, const st_idx idx_arg,
                            const dir_idx dir, hmc_float kappa_in)
{
    // this is used to save the idx of the neighbors
    st_idx idx_neigh;

    spinor out_tmp, plus;
    site_idx nn;
    su3vec psi, phi;
    Matrixsu3 U;
    // this is used to save the BC-conditions...
    hmc_complex bc_tmp = (dir == TDIR) ? (hmc_complex){kappa_in * TEMPORAL_RE, kappa_in * TEMPORAL_IM}
                                       : (hmc_complex){kappa_in * SPATIAL_RE, kappa_in * SPATIAL_IM};
    out_tmp = set_spinor_zero();

    ///////////////////////////////////
    // mu = +dir
    idx_neigh = get_neighbor_from_st_idx(idx_arg, dir);
    // transform st_idx to site_idx {space, time} -> uint
    nn   = get_site_idx(idx_neigh);
    plus = getSpinor(in, nn);
    U    = getSU3(field, get_link_idx(dir, idx_arg));
    if (dir == XDIR) {
        /////////////////////////////////
        // Calculate (1 - gamma_1) y
        // with 1 - gamma_1:
        //| 1  0  0  i |       |       psi.e0 + i*psi.e3  |
        //| 0  1  i  0 | psi = |       psi.e1 + i*psi.e2  |
        //| 0  i  1  0 |       |(-i)*( psi.e1 + i*psi.e2) |
        //| i  0  0  1 |       |(-i)*( psi.e0 + i*psi.e3) |
        /////////////////////////////////
        psi        = su3vec_acc_i(plus.e0, plus.e3);
        phi        = su3matrix_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
        out_tmp.e3 = su3vec_dim_i(out_tmp.e3, psi);

        psi        = su3vec_acc_i(plus.e1, plus.e2);
        phi        = su3matrix_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
        out_tmp.e2 = su3vec_dim_i(out_tmp.e2, psi);
    } else if (dir == YDIR) {
        ///////////////////////////////////
        // Calculate (1 - gamma_2) y
        // with 1 - gamma_2:
        // | 1  0  0  1 |       |       psi.e0 + psi.e3  |
        // | 0  1 -1  0 | psi = |       psi.e1 - psi.e2  |
        // | 0 -1  1  0 |       |(-1)*( psi.e1 + psi.e2) |
        // | 1  0  0  1 |       |     ( psi.e0 + psi.e3) |
        ///////////////////////////////////
        psi        = su3vec_acc(plus.e0, plus.e3);
        phi        = su3matrix_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
        out_tmp.e3 = su3vec_acc(out_tmp.e3, psi);

        psi        = su3vec_dim(plus.e1, plus.e2);
        phi        = su3matrix_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
        out_tmp.e2 = su3vec_dim(out_tmp.e2, psi);
    } else if (dir == ZDIR) {
        ///////////////////////////////////
        // Calculate (1 - gamma_3) y
        // with 1 - gamma_3:
        // | 1  0  i  0 |        |       psi.e0 + i*psi.e2  |
        // | 0  1  0 -i |  psi = |       psi.e1 - i*psi.e3  |
        // |-i  0  1  0 |        |   i *(psi.e0 + i*psi.e2) |
        // | 0  i  0  1 |        | (-i)*(psi.e1 - i*psi.e3) |
        ///////////////////////////////////
        psi        = su3vec_acc_i(plus.e0, plus.e2);
        phi        = su3matrix_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
        out_tmp.e2 = su3vec_dim_i(out_tmp.e2, psi);

        psi        = su3vec_dim_i(plus.e1, plus.e3);
        phi        = su3matrix_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
        out_tmp.e3 = su3vec_acc_i(out_tmp.e3, psi);
    } else {  // TDIR
              // if chemical potential is activated, U has to be multiplied by appropiate factor
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real(U, EXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp = {COSCPI, SINCPI};
        U                   = multiply_matrixsu3_by_complex(U, cpi_tmp);
#endif
        ///////////////////////////////////
        // Calculate psi/phi = (1 - gamma_0) plus/y
        // with 1 - gamma_0:
        // | 1  0  1  0 |        | psi.e0 + psi.e2 |
        // | 0  1  0  1 |  psi = | psi.e1 + psi.e3 |
        // | 1  0  1  0 |        | psi.e1 + psi.e3 |
        // | 0  1  0  1 |        | psi.e0 + psi.e2 |
        ///////////////////////////////////
        // psi = 0. component of (1-gamma_0)y
        psi = su3vec_acc(plus.e0, plus.e2);
        // phi = U*psi
        phi        = su3matrix_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
        out_tmp.e2 = su3vec_acc(out_tmp.e2, psi);
        // psi = 1. component of (1-gamma_0)y
        psi = su3vec_acc(plus.e1, plus.e3);
        // phi = U*psi
        phi        = su3matrix_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
        out_tmp.e3 = su3vec_acc(out_tmp.e3, psi);
    }

    ///////////////////////////////////
    // mu = -dir
    idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir);
    // transform st_idx to site_idx {space, time} -> uint
    nn   = get_site_idx(idx_neigh);
    plus = getSpinor(in, nn);
    U    = getSU3(field, get_link_idx(dir, idx_neigh));
    // in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
    bc_tmp = (dir == TDIR) ? (hmc_complex){kappa_in * TEMPORAL_RE, kappa_in * MTEMPORAL_IM}
                           : (hmc_complex){kappa_in * SPATIAL_RE, kappa_in * MSPATIAL_IM};
    if (dir == XDIR) {
        ///////////////////////////////////
        // Calculate (1 + gamma_1) y
        // with 1 + gamma_1:
        // | 1  0  0 -i |       |       psi.e0 - i*psi.e3  |
        // | 0  1 -i  0 | psi = |       psi.e1 - i*psi.e2  |
        // | 0  i  1  0 |       |(-i)*( psi.e1 - i*psi.e2) |
        // | i  0  0  1 |       |(-i)*( psi.e0 - i*psi.e3) |
        ///////////////////////////////////
        psi        = su3vec_dim_i(plus.e0, plus.e3);
        phi        = su3matrix_dagger_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
        out_tmp.e3 = su3vec_acc_i(out_tmp.e3, psi);

        psi        = su3vec_dim_i(plus.e1, plus.e2);
        phi        = su3matrix_dagger_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
        out_tmp.e2 = su3vec_acc_i(out_tmp.e2, psi);
    } else if (dir == YDIR) {
        ///////////////////////////////////
        // Calculate (1 + gamma_2) y
        // with 1 + gamma_2:
        // | 1  0  0 -1 |       |       psi.e0 - psi.e3  |
        // | 0  1  1  0 | psi = |       psi.e1 + psi.e2  |
        // | 0  1  1  0 |       |     ( psi.e1 + psi.e2) |
        // |-1  0  0  1 |       |(-1)*( psi.e0 - psi.e3) |
        ///////////////////////////////////
        psi        = su3vec_dim(plus.e0, plus.e3);
        phi        = su3matrix_dagger_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
        out_tmp.e3 = su3vec_dim(out_tmp.e3, psi);

        psi        = su3vec_acc(plus.e1, plus.e2);
        phi        = su3matrix_dagger_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
        out_tmp.e2 = su3vec_acc(out_tmp.e2, psi);
    } else if (dir == ZDIR) {
        ///////////////////////////////////
        // Calculate (1 + gamma_3) y
        // with 1 + gamma_3:
        // | 1  0 -i  0 |       |       psi.e0 - i*psi.e2  |
        // | 0  1  0  i | psi = |       psi.e1 + i*psi.e3  |
        // | i  0  1  0 |       | (-i)*(psi.e0 - i*psi.e2) |
        // | 0 -i  0  1 |       |   i *(psi.e1 + i*psi.e3) |
        ///////////////////////////////////
        psi        = su3vec_dim_i(plus.e0, plus.e2);
        phi        = su3matrix_dagger_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
        out_tmp.e2 = su3vec_acc_i(out_tmp.e2, psi);

        psi        = su3vec_acc_i(plus.e1, plus.e3);
        phi        = su3matrix_dagger_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
        out_tmp.e3 = su3vec_dim_i(out_tmp.e3, psi);
    } else {  // TDIR
              // if chemical potential is activated, U has to be multiplied by appropiate factor
              // this is the same as at mu=0 in the imag. case, since U is taken to be U^+ later:
              //  (exp(iq)U)^+ = exp(-iq)U^+
              // as it should be
              // in the real case, one has to take exp(q) -> exp(-q)
#ifdef _CP_REAL_
        U = multiply_matrixsu3_by_real(U, MEXPCPR);
#endif
#ifdef _CP_IMAG_
        hmc_complex cpi_tmp2 = {COSCPI, SINCPI};
        U                    = multiply_matrixsu3_by_complex(U, cpi_tmp2);
#endif
        ///////////////////////////////////
        // Calculate psi/phi = (1 + gamma_0) y
        // with 1 + gamma_0:
        // | 1  0 -1  0 |       | psi.e0 - psi.e2 |
        // | 0  1  0 -1 | psi = | psi.e1 - psi.e3 |
        // |-1  0  1  0 |       | psi.e1 - psi.e2 |
        // | 0 -1  0  1 |       | psi.e0 - psi.e3 |
        ///////////////////////////////////
        psi = su3vec_dim(plus.e0, plus.e2);
        // phi = U*psi
        phi        = su3matrix_dagger_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
        out_tmp.e2 = su3vec_dim(out_tmp.e2, psi);
        // psi = 1. component of (1+gamma_0)y
        psi = su3vec_dim(plus.e1, plus.e3);
        // phi = U*psi
        phi        = su3matrix_dagger_times_su3vec(U, psi);
        psi        = su3vec_times_complex(phi, bc_tmp);
        out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
        out_tmp.e3 = su3vec_dim(out_tmp.e3, psi);
    }

    return out_tmp;
}
