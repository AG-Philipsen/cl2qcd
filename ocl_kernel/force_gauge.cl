/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

/** @file
 * Kernel for the standard Wilson-action gauge force calculation.
 * 
 * According to the Gattringer book (page 198) the gauge part of the force should be
 * @code
 * F_G = - \beta/6 * \sum_{k=1}^8 T_k Tr[i * T_k * (U*V - V^\dag*U^\dag)] 
 * 
 *     = - \beta/12 * i * (U*V - V^\dag*U^\dag)
 * @endcode
 * where V are the well known staples. Nevertheless this expression is not in
 * agreement with the Gottlieb and Toussaint work (Hybrid-molecular dynamics
 * algorithm for numerical simulation of QCD), where we found
 * @code
 * F_G = - \beta/6 * i * (U*V - V^\dag*U^\dag)
 * 
 *     = - \beta/3 * \sum_{k=1}^8 T_k Tr[i * T_k * (U*V - V^\dag*U^\dag)]  .
 * @endcode
 * Then there is a factor 1/2 of discrepancy. Actually, this is only apparent:
 * Gattringer and Gottlieb use different factors in the gaugemomenta part of the action
 * (compare eq.(8.35) of the first with eq.(4) of the second).
 * 
 * @note It is worth recalling that an overall factor in the force calculation does
 *       not affect the correctness of the (R)HMC algorithm. In fact (thinking to
 *       leapfrog, see page 197 Gattringer), this factor would mean to have a different
 *       time step size in the integration of the equations of motion. Of course if one
 *       wants to make the fields evolve for a fixed time in, let's say, 100 steps,
 *       he will deduce a certain time step. A missing overall factor in the force
 *       calculation would affect the time step and therefore all the fields would
 *       evolve for less or more time!
 * 
 * @attention Now, here, as one can see reading the code, the matrix U*V is calculated
 *            and then the function tr_lambda_u is called on it. This function returns
 *            Tr[i * T_k * (U*V - V^\dag*U^\dag)] that is multiplied ONLY by @f$ -\beta/3 @f$.
 *            We are then apparently coeherent with the Gottlieb and Toussaint result.
 *            Actually, this is not true because the gaugemomenta part of the action calculated
 *            in physics/algorithms/metropolis.cpp (in the Wilson case) is exactly that
 *            in eq.(8.35) of Gattringer book: as explained in the note above, there can be
 *            an overall factor in the force calculation that has to be balanced by an
 *            appropriate choice of the time interval (i.e. the parameter "tau" of Inputparameters),
 *            if one wants to integrate the MD equations of motion for a given interval (in
 *            other words, if there is an overall factor k in the force calculation, and the eq.
 *            of motion should be integrated for a time interval Dt, then tau must be Dt/k).
 *            In the staggered formulation, instead, all the code is completely coherent with
 *            the Gottlieb-Toussaint paper and this factors compensation has been avoided (the
 *            gaugemomenta part of the action calculated in physics/algorithms/metropolis.cpp
 *            is different from the Wilson case by a factor 0.5). The force calculated in the
 *            code is exactly "the right one", namely that which can be analitically derived
 *            from the action. If the eq. of motion should be integrated for a time interval Dt,
 *            then tau must be Dt.
 *            
 */

inline void gauge_force_per_link(__global const Matrixsu3StorageType * const restrict field, __global aeStorageType * const restrict out, const st_index pos, const dir_idx dir) {
	Matrix3x3 V = calc_staple(field, pos.space, pos.time, dir);
	Matrixsu3 U = get_matrixsu3(field, pos.space, pos.time, dir);
	V = multiply_matrix3x3 (matrix_su3to3x3(U), V);
	ae out_tmp = tr_lambda_u(V);

	hmc_float factor = -BETA / 3.;
#ifdef _USE_RECT_
	factor = factor * C0;
#endif
	int global_link_pos = get_link_idx(dir, pos);
	update_gaugemomentum(out_tmp, factor , global_link_pos, out);
}

__kernel void gauge_force(__global const Matrixsu3StorageType * const restrict field, __global aeStorageType * const restrict out)
{
	//Gauge force is factor*Im(Tr(T_i U V))
	//   with T_i being the SU3-Generator in i-th direction and V the staplematrix
	//   and the factor being beta / NC (for standard Wilson-action) and c0 * beta / NC (for tlSym)
	PARALLEL_FOR(id_local, VOL4D_LOCAL * NDIM) {
		//calc link-pos and mu out of the index
		//NOTE: this is not necessarily equal to the geometric  conventions, one just needs a one-to-one correspondence between thread-id and (n,t,mu) here
		const int pos_local = id_local % VOL4D_LOCAL;
		const int dir     = id_local / VOL4D_LOCAL;
		const st_index pos = (pos_local >= VOL4D_LOCAL / 2) ? get_even_st_idx_local(pos_local - (VOL4D_LOCAL / 2)) : get_odd_st_idx_local(pos_local);

		gauge_force_per_link(field, out, pos, dir);
	}
}
