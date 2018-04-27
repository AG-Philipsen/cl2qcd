/*
 * Copyright (c) 2013,2018 Alessandro Sciarra
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

/** @file
 * M normal staggered fermionmatrix
 *
 * In order to clarify notation and not to confuse things, let us recall some expressions.
 * In the staggered formulation, after the "staggering process", we remain with one component
 * field per site (an su3vec per site if we take into account the colour degree of freedom).
 * Since in an HMC algorithm we do not deal directly with such one-component field (in fact,
 * first we integrate over the fermionic field getting the fermionic determinant and then
 * we introduce a pseudofermionic field to rewrite the fermionic determinat in a suitable way),
 * here we will focus solely on the Dirac operator that we will indicate with M. Referring for
 * instance to the Rothe (Lattice Gauge Theories, Eq.(6.26), pag 94), we have:
 *
 * \f[
   M_{n,m}=\frac{1}{2} \sum_{\mu=0}^3 \eta_\mu(n)\Bigl[U_\mu(n)\,\delta_{n+\hat\mu,m}
                                                     - U^\dag_\mu(n-\hat\mu)\,\delta_{n-\hat\mu,m}\Bigr]
                                                     + m_0\,\delta_{n,m}
   \f]
 *
 * \internal
 *
 *               ___            _                                                     _
 *            1  \             |                                                       |
 * M_{n,m} = --- /__  eta_mu(n)| U_mu(n) delta_{n+mu,m} - U^dag_mu(n-mu) delta_{n-mu,m}| + m delta_{n,m}
 *            2   mu           |_                                                     _|
 *
 *
 * We refer to M as M=m+D_KS
 *
 * \endinternal
 *
 * Here in this kernel we make the M operator act on the field "in" getting the result in the
 * field "out". Note that the field "field" is the gauge field (i.e. the link variables) that
 * are needed in the functions dslash_local_x.
 *
 * @note This kernel takes into account the staggered phases eta_mu(n), since they
 *       will be NOT included in links (as in the phi-algorithm from Gottlieb and Toussaint).
 *       The reason why we make such a choice is that we can preserve gauge configurations
 *       (and all related code) as it is, without changing the links in the staggered algorithm.
 * \par
 * @note In this kernel \em no \em even-odd \em preconditioning is assumed: vectors "in" and "out" have
 *       VOL4D components.
 */

__kernel void M_staggered(__global const su3vec* const restrict in,
                          __global const Matrixsu3StorageType* const restrict field,
                          __global su3vec* const restrict out, hmc_float mass_in)
{
    int global_size = get_global_size(0);
    int id          = get_global_id(0);
    su3vec out_tmp, out_tmp2;

    for (int id_local = id; id_local < SPINORFIELDSIZE_LOCAL; id_local += global_size) {
        /** @todo this must be done more efficient */
        st_index pos = (id_local % 2 == 0) ? get_even_st_idx_local(id_local / 2) : get_odd_st_idx_local(id_local / 2);

        // From now on we adopt the notation M = m + D_KS

        // Diagonal part: m * in(n)
        out_tmp = get_su3vec_from_field(in, pos.space, pos.time);
        out_tmp = su3vec_times_real(out_tmp, mass_in);

        // Non-diagonal part: calc D_KS
        for (dir_idx dir = 0; dir < 4; ++dir) {
            out_tmp2 = D_KS_local(in, field, pos, dir);
            out_tmp  = su3vec_acc(out_tmp, out_tmp2);
        }

        put_su3vec_to_field(out_tmp, out, pos.space, pos.time);
    }
}
