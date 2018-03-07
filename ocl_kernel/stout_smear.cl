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
 * @file stout-smearing of the gaugefield
 */

//this is done after Phys. Rev. D69, 0545501 (Morningstar, Peardon) and the tmlqcd-analogue (stout_smear.c)
__kernel void stout_smear( __global Matrixsu3StorageType * in, __global Matrixsu3StorageType * out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	PARALLEL_FOR(id_local, VOL4D_LOCAL * NDIM) {
		//calc link-pos and mu out of the index
		//NOTE: this is not necessarily equal to the geometric conventions, one just needs a one-to-one correspondence between thread-id and (n,t,mu) here
		int2 pos_tmp;
		pos_tmp.x = id_local % VOL4D_LOCAL;
		//this is mu
		pos_tmp.y = id_local / VOL4D_LOCAL;

		st_index pos = (pos_tmp.x % 2 == 0) ? get_even_st_idx_local(pos_tmp.x / 2) : get_odd_st_idx_local(pos_tmp.x / 2);

		//calc staple
		Matrix3x3 staple;
		//NOTE: this calculates the staple as it appears in the Wilson
		//  gauge action. However, here one is interested in replacing
		//  a link at site x in direction mu by a staple that
		//  starts at x and in the end also "points" at x+mu.
		//  This means that one has to take the adjoint of the
		//  staple! (see also fig. 1 in the paper)
		staple = calc_staple(in, pos.space, pos.time, pos_tmp.y) ;

		//staple*rho
		staple = multiply_matrix3x3_by_real(staple, RHO);

		//omega = staple * u^dagger
		//  in our case this means then omega = staple^dagger * u^dagger
		Matrixsu3 u;
		u = get_matrixsu3(in, pos.space, pos.time, pos_tmp.y);

		Matrix3x3 omega;
		omega = multiply_matrix3x3_dagger_dagger(staple, matrix_su3to3x3(u));

		//project out anti-hermitian traceless part
		Matrixsu3 tmp;
		tmp = project_anti_herm(omega);
		Matrix3x3 tmp2 = matrix_su3to3x3(tmp);

		//exponentiate
		// 1/-2.0 is inserted to get su3 to su3adjoint consistency
		ae p = tr_lambda_u(tmp2) ;
		Matrixsu3 expp = build_su3matrix_by_exponentiation(p, -.5);

		//new_u = expp * u
		Matrixsu3 new_u = multiply_matrixsu3(expp, u);
		put_matrixsu3(out, new_u, pos.space, pos.time, pos_tmp.y);
	}
}
