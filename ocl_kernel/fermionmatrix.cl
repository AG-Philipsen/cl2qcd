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

//"local" dslash working on a particular link (n,t) in a specific direction
// NOTE: each component is multiplied by +KAPPA, so the resulting spinor has to be mutliplied by -1 to obtain the
// correct dslash!!! spinor dslash_local_0(__global const spinorfield * const restrict in,__global const
// ocl_s_gaugefield
// * const restrict field, const int n, const int t){
spinor dslash_local_0(__global const spinor* const restrict in,
                      __global const Matrixsu3StorageType* const restrict field, int n, int t, hmc_float kappa_in)
{
    spinor out_tmp, plus;
    int dir, nn;
    su3vec psi, phi;
    Matrixsu3 U;
    // this is used to save the BC-conditions...
    hmc_complex bc_tmp;
    out_tmp = set_spinor_zero();

    // go through the different directions
    ///////////////////////////////////
    // mu = 0
    ///////////////////////////////////
    dir = 0;
    ///////////////////////////////////
    // mu = +0
    nn   = get_neighbor_temporal(t);
    plus = get_spinor_from_field(in, n, nn);
    U    = getSU3(field, get_link_pos(dir, n, t));
    // if chemical potential is activated, U has to be multiplied by appropiate factor
#ifdef _CP_REAL_
    U = multiply_matrixsu3_by_real(U, EXPCPR);
#endif
#ifdef _CP_IMAG_
    hmc_complex cpi_tmp = {COSCPI, SINCPI};
    U                   = multiply_matrixsu3_by_complex(U, cpi_tmp);
#endif
    bc_tmp.re = kappa_in * TEMPORAL_RE;
    bc_tmp.im = kappa_in * TEMPORAL_IM;
    ///////////////////////////////////
    // Calculate psi/phi = (1 - gamma_0) plus/y
    // with 1 - gamma_0:
    // | 1  0  1  0 |        | psi.e0 + psi.e2 |
    // | 0  1  0  1 |  psi = | psi.e1 + psi.e3 |
    // | 1  0  1  0 |        | psi.e0 + psi.e2 |
    // | 0  1  0  1 |        | psi.e1 + psi.e3 |
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

    /////////////////////////////////////
    // mu = -0
    nn   = get_lower_neighbor_temporal(t);
    plus = get_spinor_from_field(in, n, nn);
    U    = getSU3(field, get_link_pos(dir, n, nn));
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
    // in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
    bc_tmp.re = kappa_in * TEMPORAL_RE;
    bc_tmp.im = kappa_in * MTEMPORAL_IM;
    ///////////////////////////////////
    // Calculate psi/phi = (1 + gamma_0) y
    // with 1 + gamma_0:
    // | 1  0 -1  0 |       | psi.e0 - psi.e2 |
    // | 0  1  0 -1 | psi = | psi.e1 - psi.e3 |
    // |-1  0  1  0 |       | psi.e2 - psi.e0 |
    // | 0 -1  0  1 |       | psi.e3 - psi.e1 |
    ///////////////////////////////////
    // psi = 0. component of (1+gamma_0)y
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

    return out_tmp;
}

spinor dslash_local_1(__global const spinor* const restrict in,
                      __global const Matrixsu3StorageType* const restrict field, int n, int t, hmc_float kappa_in)
{
    spinor out_tmp, plus;
    int dir, nn;
    su3vec psi, phi;
    Matrixsu3 U;
    // this is used to save the BC-conditions...
    hmc_complex bc_tmp;
    out_tmp = set_spinor_zero();

    // CP: all actions correspond to the mu = 0 ones
    ///////////////////////////////////
    // mu = 1
    ///////////////////////////////////
    dir = 1;

    ///////////////////////////////////
    // mu = +1
    nn        = get_neighbor(n, dir);
    plus      = get_spinor_from_field(in, nn, t);
    U         = getSU3(field, get_link_pos(dir, n, t));
    bc_tmp.re = kappa_in * SPATIAL_RE;
    bc_tmp.im = kappa_in * SPATIAL_IM;
    /////////////////////////////////
    // Calculate (1 - gamma_1) y
    // with 1 - gamma_1:
    //|  1  0  0  i |       |       psi.e0 + i*psi.e3  |
    //|  0  1  i  0 | psi = |       psi.e1 + i*psi.e2  |
    //|  0 -i  1  0 |       |(-i)*( psi.e1 + i*psi.e2) |
    //| -i  0  0  1 |       |(-i)*( psi.e0 + i*psi.e3) |
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

    ///////////////////////////////////
    // mu = -1
    nn   = get_lower_neighbor(n, dir);
    plus = get_spinor_from_field(in, nn, t);
    U    = getSU3(field, get_link_pos(dir, nn, t));
    // in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
    bc_tmp.re = kappa_in * SPATIAL_RE;
    bc_tmp.im = kappa_in * MSPATIAL_IM;
    ///////////////////////////////////
    // Calculate (1 + gamma_1) y
    // with 1 + gamma_1:
    // | 1  0  0 -i |       |       psi.e0 - i*psi.e3  |
    // | 0  1 -i  0 | psi = |       psi.e1 - i*psi.e2  |
    // | 0  i  1  0 |       | (i)*( psi.e1 - i*psi.e2) |
    // | i  0  0  1 |       | (i)*( psi.e0 - i*psi.e3) |
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

    return out_tmp;
}

spinor dslash_local_2(__global const spinor* const restrict in,
                      __global const Matrixsu3StorageType* const restrict field, int n, int t, hmc_float kappa_in)
{
    spinor out_tmp, plus;
    int dir, nn;
    su3vec psi, phi;
    Matrixsu3 U;
    // this is used to save the BC-conditions...
    hmc_complex bc_tmp;
    out_tmp = set_spinor_zero();
    ;

    ///////////////////////////////////
    // mu = 2
    ///////////////////////////////////
    dir = 2;

    ///////////////////////////////////
    // mu = +2
    nn        = get_neighbor(n, dir);
    plus      = get_spinor_from_field(in, nn, t);
    U         = getSU3(field, get_link_pos(dir, n, t));
    bc_tmp.re = kappa_in * SPATIAL_RE;
    bc_tmp.im = kappa_in * SPATIAL_IM;
    ///////////////////////////////////
    // Calculate (1 - gamma_2) y
    // with 1 - gamma_2:
    // | 1  0  0  1 |       |       psi.e0 + psi.e3  |
    // | 0  1 -1  0 | psi = |       psi.e1 - psi.e2  |
    // | 0 -1  1  0 |       |(-1)*( psi.e1 - psi.e2) |
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

    ///////////////////////////////////
    // mu = -2
    nn   = get_lower_neighbor(n, dir);
    plus = get_spinor_from_field(in, nn, t);
    U    = getSU3(field, get_link_pos(dir, nn, t));
    // in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
    bc_tmp.re = kappa_in * SPATIAL_RE;
    bc_tmp.im = kappa_in * MSPATIAL_IM;
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

    return out_tmp;
}

spinor dslash_local_3(__global const spinor* const restrict in,
                      __global const Matrixsu3StorageType* const restrict field, int n, int t, hmc_float kappa_in)
{
    spinor out_tmp, plus;
    int dir, nn;
    su3vec psi, phi;
    Matrixsu3 U;
    // this is used to save the BC-conditions...
    hmc_complex bc_tmp;
    out_tmp = set_spinor_zero();

    ///////////////////////////////////
    // mu = 3
    ///////////////////////////////////
    dir = 3;

    ///////////////////////////////////
    // mu = +3
    nn        = get_neighbor(n, dir);
    plus      = get_spinor_from_field(in, nn, t);
    U         = getSU3(field, get_link_pos(dir, n, t));
    bc_tmp.re = kappa_in * SPATIAL_RE;
    bc_tmp.im = kappa_in * SPATIAL_IM;
    ///////////////////////////////////
    // Calculate (1 - gamma_3) y
    // with 1 - gamma_3:
    // | 1  0  i  0 |        |       psi.e0 + i*psi.e2  |
    // | 0  1  0 -i |  psi = |       psi.e1 - i*psi.e3  |
    // |-i  0  1  0 |        | (-i)*(psi.e0 + i*psi.e2) |
    // | 0  i  0  1 |        | ( i)*(psi.e1 - i*psi.e3) |
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

    ///////////////////////////////////
    // mu = -3
    nn   = get_lower_neighbor(n, dir);
    plus = get_spinor_from_field(in, nn, t);
    U    = getSU3(field, get_link_pos(dir, nn, t));
    // in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
    bc_tmp.re = kappa_in * SPATIAL_RE;
    bc_tmp.im = kappa_in * MSPATIAL_IM;
    ///////////////////////////////////
    // Calculate (1 + gamma_3) y
    // with 1 + gamma_3:
    // | 1  0 -i  0 |       |       psi.e0 - i*psi.e2  |
    // | 0  1  0  i | psi = |       psi.e1 + i*psi.e3  |
    // | i  0  1  0 |       | ( i)*(psi.e0 - i*psi.e2) |
    // | 0 -i  0  1 |       | (-i)*(psi.e1 + i*psi.e3) |
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

    return out_tmp;
}
