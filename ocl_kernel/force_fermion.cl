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
 * @file kernel for the non-eo fermion force
 */

__kernel void fermion_force(__global const Matrixsu3StorageType* const restrict field,
                            __global const spinor* const restrict Y, __global const spinor* const restrict X,
                            __global aeStorageType* const restrict out, const hmc_float kappa_in)
{
    for (dir_idx dir = 0; dir < NDIM; ++dir) {
        PARALLEL_FOR (id_local, VOL4D_LOCAL) {
            st_index pos             = get_st_idx_from_site_idx(id_local);
            link_idx global_link_pos = get_link_idx(dir, pos);

            Matrixsu3 U;
            Matrix3x3 v1, v2, tmp;
            su3vec psia, psib, phia, phib;
            spinor y, plus;
            int nn;
            // this is used to save the BC-conditions...
            hmc_complex bc_tmp;
            int n = pos.space;
            int t = pos.time;

            ae ae_tmp = getAe(out, global_link_pos);

            if (dir == 0) {
                y = getSpinor(Y, get_pos(n, t));
                ///////////////////////////////////
                // Calculate gamma_5 y
                ///////////////////////////////////
                y = gamma5_local(y);

                ///////////////////////////////////
                // mu = 0
                ///////////////////////////////////
                nn = get_neighbor_temporal(t);

                // the 2 here comes from Tr(lambda_ij) = 2delta_ij
                bc_tmp.re = 2. * kappa_in * TEMPORAL_RE;
                bc_tmp.im = 2. * kappa_in * TEMPORAL_IM;

                ///////////////////////////////////
                // mu = +0
                plus = getSpinor(X, get_pos(n, nn));
                U    = get_matrixsu3(field, n, t, dir);
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
                // here, only the independent components are needed
                // the other ones are incorporated in the dirac trace
                // psia = 0. component of (1-gamma_0)plus
                psia = su3vec_acc(plus.e0, plus.e2);
                // psib = 1. component of (1-gamma_0)plus
                psib = su3vec_acc(plus.e1, plus.e3);
                phia = su3vec_acc(y.e0, y.e2);
                phib = su3vec_acc(y.e1, y.e3);
                // v1 = Tr(phi*psi_dagger)
                v1 = tr_v_times_u_dagger(phia, psia, phib, psib);
                // U*v1 = U*(phi_a)
                tmp = matrix_su3to3x3(U);
                v2  = multiply_matrix3x3_dagger(tmp, v1);
                v1  = multiply_matrix3x3_by_complex(v2, bc_tmp);

                ae_tmp = acc_factor_times_algebraelement(ae_tmp, 1., tr_lambda_u(v1));

                /////////////////////////////////////
                // mu = -0
                y    = getSpinor(Y, get_pos(n, nn));
                y    = gamma5_local(y);
                plus = getSpinor(X, get_pos(n, t));
                // if chemical potential is activated, U has to be multiplied by appropiate factor
                // this is the same as at mu=0 in the imag. case, since U is taken to be U^+ later:
                //  (exp(iq)U)^+ = exp(-iq)U^+
                // as it should be
                // in the real case, one has to take exp(q) -> exp(-q)
#if defined(_CP_REAL_) || defined(_CP_IMAG_)
                U = get_matrixsu3(field, n, t, dir);
#endif
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
                psia = su3vec_dim(plus.e0, plus.e2);
                psib = su3vec_dim(plus.e1, plus.e3);
                phia = su3vec_dim(y.e0, y.e2);
                phib = su3vec_dim(y.e1, y.e3);
                // CP: here is the difference with regard to +mu-direction: psi and phi interchanged!!
                // v1 = Tr(psi*phi_dagger)
                v1 = tr_v_times_u_dagger(psia, phia, psib, phib);
                // U*v1 = U*(phi_a)
                v2 = multiply_matrix3x3_dagger(tmp, v1);
                v1 = multiply_matrix3x3_by_complex(v2, bc_tmp);

                ae_tmp = acc_factor_times_algebraelement(ae_tmp, 1., tr_lambda_u(v1));

            } else {
                y = getSpinor(Y, get_pos(n, t));
                ///////////////////////////////////
                // Calculate gamma_5 y
                ///////////////////////////////////
                y = gamma5_local(y);

                // comments correspond to the mu=0 ones
                /////////////////////////////////
                // mu = 1
                /////////////////////////////////
                // this stays the same for all spatial directions at the moment
                bc_tmp.re = 2. * kappa_in * SPATIAL_RE;
                bc_tmp.im = 2. * kappa_in * SPATIAL_IM;

                /////////////////////////////////
                // mu = +1
                nn   = get_neighbor(n, dir);
                plus = getSpinor(X, get_pos(nn, t));
                U    = get_matrixsu3(field, n, t, dir);

                if (dir == 1) {
                    /////////////////////////////////
                    // Calculate (1 - gamma_1) y
                    // with 1 - gamma_1:
                    //| 1  0  0  i |       |       psi.e0 + i*psi.e3  |
                    //| 0  1  i  0 | psi = |       psi.e1 + i*psi.e2  |
                    //| 0  i  1  0 |       |(-i)*( psi.e1 + i*psi.e2) |
                    //| i  0  0  1 |       |(-i)*( psi.e0 + i*psi.e3) |
                    /////////////////////////////////
                    // psi = (1-gamma_mu)plus
                    psia = su3vec_acc_i(plus.e0, plus.e3);
                    psib = su3vec_acc_i(plus.e1, plus.e2);
                    phia = su3vec_acc_i(y.e0, y.e3);
                    phib = su3vec_acc_i(y.e1, y.e2);
                } else if (dir == 2) {
                    ///////////////////////////////////
                    // Calculate (1 - gamma_2) y
                    // with 1 - gamma_2:
                    // | 1  0  0  1 |       |       psi.e0 + psi.e3  |
                    // | 0  1 -1  0 | psi = |       psi.e1 - psi.e2  |
                    // | 0 -1  1  0 |       |(-1)*( psi.e1 + psi.e2) |
                    // | 1  0  0  1 |       |     ( psi.e0 + psi.e3) |
                    ///////////////////////////////////
                    psia = su3vec_acc(plus.e0, plus.e3);
                    psib = su3vec_dim(plus.e1, plus.e2);
                    phia = su3vec_acc(y.e0, y.e3);
                    phib = su3vec_dim(y.e1, y.e2);
                } else {  // dir == 3
                    ///////////////////////////////////
                    // Calculate (1 - gamma_3) y
                    // with 1 - gamma_3:
                    // | 1  0  i  0 |        |       psi.e0 + i*psi.e2  |
                    // | 0  1  0 -i |  psi = |       psi.e1 - i*psi.e3  |
                    // |-i  0  1  0 |        |   i *(psi.e0 + i*psi.e2) |
                    // | 0  i  0  1 |        | (-i)*(psi.e1 - i*psi.e3) |
                    ///////////////////////////////////
                    psia = su3vec_acc_i(plus.e0, plus.e2);
                    psib = su3vec_dim_i(plus.e1, plus.e3);
                    phia = su3vec_acc_i(y.e0, y.e2);
                    phib = su3vec_dim_i(y.e1, y.e3);
                }

                v1  = tr_v_times_u_dagger(phia, psia, phib, psib);
                tmp = matrix_su3to3x3(U);
                v2  = multiply_matrix3x3_dagger(tmp, v1);
                v1  = multiply_matrix3x3_by_complex(v2, bc_tmp);

                ae_tmp = acc_factor_times_algebraelement(ae_tmp, 1., tr_lambda_u(v1));

                ///////////////////////////////////
                // mu = -1
                y    = getSpinor(Y, get_pos(nn, t));
                y    = gamma5_local(y);
                plus = getSpinor(X, get_pos(n, t));

                if (dir == 1) {
                    ///////////////////////////////////
                    // Calculate (1 + gamma_1) y
                    // with 1 + gamma_1:
                    // | 1  0  0 -i |       |       psi.e0 - i*psi.e3  |
                    // | 0  1 -i  0 | psi = |       psi.e1 - i*psi.e2  |
                    // | 0  i  1  0 |       |(-i)*( psi.e1 - i*psi.e2) |
                    // | i  0  0  1 |       |(-i)*( psi.e0 - i*psi.e3) |
                    ///////////////////////////////////
                    psia = su3vec_dim_i(plus.e0, plus.e3);
                    psib = su3vec_dim_i(plus.e1, plus.e2);
                    phia = su3vec_dim_i(y.e0, y.e3);
                    phib = su3vec_dim_i(y.e1, y.e2);
                } else if (dir == 2) {
                    ///////////////////////////////////
                    // Calculate (1 + gamma_2) y
                    // with 1 + gamma_2:
                    // | 1  0  0 -1 |       |       psi.e0 - psi.e3  |
                    // | 0  1  1  0 | psi = |       psi.e1 + psi.e2  |
                    // | 0  1  1  0 |       |     ( psi.e1 + psi.e2) |
                    // |-1  0  0  1 |       |(-1)*( psi.e0 - psi.e3) |
                    ///////////////////////////////////
                    psia = su3vec_dim(plus.e0, plus.e3);
                    psib = su3vec_acc(plus.e1, plus.e2);
                    phia = su3vec_dim(y.e0, y.e3);
                    phib = su3vec_acc(y.e1, y.e2);
                } else {  // dir == 3
                    ///////////////////////////////////
                    // Calculate (1 + gamma_3) y
                    // with 1 + gamma_3:
                    // | 1  0 -i  0 |       |       psi.e0 - i*psi.e2  |
                    // | 0  1  0  i | psi = |       psi.e1 + i*psi.e3  |
                    // | i  0  1  0 |       | (-i)*(psi.e0 - i*psi.e2) |
                    // | 0 -i  0  1 |       |   i *(psi.e1 + i*psi.e3) |
                    ///////////////////////////////////
                    psia = su3vec_dim_i(plus.e0, plus.e2);
                    psib = su3vec_acc_i(plus.e1, plus.e3);
                    phia = su3vec_dim_i(y.e0, y.e2);
                    phib = su3vec_acc_i(y.e1, y.e3);
                }

                v1 = tr_v_times_u_dagger(psia, phia, psib, phib);
                v2 = multiply_matrix3x3_dagger(tmp, v1);
                v1 = multiply_matrix3x3_by_complex(v2, bc_tmp);

                ae_tmp = acc_factor_times_algebraelement(ae_tmp, 1., tr_lambda_u(v1));
            }
            putAe(out, global_link_pos, ae_tmp);
        }
    }
}
