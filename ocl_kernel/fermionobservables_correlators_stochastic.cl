/*
 * Copyright (c) 2012,2013 Matthias Bach
 * Copyright (c) 2012 Christopher Pinke
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
 @file fermion-observables: correlators calculated from stochastic sources
*/

/**
 * For stochastic sources one has:
 * D_u^-1(x,y)_(a,b) = 1/R Phi^r(x)_(a) * ( Xi^r(y)_(b) )^dagger
 * with a,b : colour-Dirac indices
 * with x,y : space-time indices
 * with r : index of source
 */

/**
 * Nf2 Pseudoscalar Correlator (Pion):
 * C(x,y)_(a,b) = <Pi^+(x)_a Pi^-(y)_b> = <u(x)_a gamma5 d^bar(x)_a u^bar(y)_b gamma5 d(y)_b >
 *              = <u(x)_a u^bar(y)_b gamma5 d^bar(x)_a d(y)_b gamma5 >
 *              = <D_u^-1(x,y)_(a,b) gamma5 D_d^-1(y,x)_(b,a) gamma5 >
 *              = <D_u^-1(x,y)_(a,b) (D_u^-1(x,y)_(a,b) )^dagger >
 *              = <|D_u^-1(x,y)_(a,b)|^2>
 */

/**
 * Nf2 Pion Norm:
 * |Pi|^2 = Sum_{x,y,a,b} C(x,y)_(a,b) = Sum_{x,y,a,b} <|D_u^-1(x,y)_(a,b)|^2>
 *        =
 */

/**
 * This is the pseudoscalar pion correlator in z-direction:
 * C(z) = Sum_{x,y,t,a,b} <|D_u^-1( x(z), y(z) )_(a,b)|^2>
 *      = Sum_{x,y,t,a,b,r} 1/R <| Phi^r(x)_(a) * ( Xi^r(y)_(b) )^dagger  |^2>
 */

inline hmc_float complex_squarenorm(const hmc_complex a)
{
	hmc_float res;
	res = a.re * a.re + a.im * a.im;
	return res;
}

void load_spinor_to_complex_array(spinor in, hmc_complex * out)
{
	out[0] = in.e0.e0;
	out[1] = in.e0.e1;
	out[2] = in.e0.e2;
	out[3] = in.e1.e0;
	out[4] = in.e1.e1;
	out[5] = in.e1.e2;
	out[6] = in.e2.e0;
	out[7] = in.e2.e1;
	out[8] = in.e2.e2;
	out[9] = in.e3.e0;
	out[10] = in.e3.e1;
	out[11] = in.e3.e2;
}


__kernel void correlator_ps_z(__global hmc_float * const restrict out, __global const spinor * const restrict phi, __global const spinor * const restrict b)
{
	const int local_size = get_local_size(0);
	const int global_size = get_global_size(0);
	const int id = get_global_id(0);
	const int loc_idx = get_local_id(0);
	const int num_groups = get_num_groups(0);
	const int group_id = get_group_id (0);

	//now, this should finally be the correct normalisation for the physical fields
	//one factor of 2*kappa per field and we construct the correlator from a multiplication of two fields phi
	const hmc_float fac = 2. * KAPPA * 2.* KAPPA / (NSPACE * NSPACE * NTIME_GLOBAL);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_float correlator = out[id_tmp];
		uint3 coord;
		uint3 coord2;
		//loop over first coordinate
		coord.z = id_tmp;
		for(int t    = 0; t < NTIME_LOCAL;  t++) {
			for(coord.x  = 0; coord.x <  NSPACE; coord.x++) {
				for(coord.y  = 0; coord.y <  NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor phi_tmp = phi[get_pos(nspace, t)];
					hmc_complex phi_arr[12];
					load_spinor_to_complex_array(phi_tmp, phi_arr);
					//loop over second coordinate
					for(int t2   = 0; t2 < NTIME_LOCAL;  t2++) { // TODO needs to be worked around for Multi-GPU
						for(coord2.x = 0; coord2.x < NSPACE; coord2.x++) {
							for(coord2.y = 0; coord2.y < NSPACE; coord2.y++) {
								for(coord2.z = 0; coord2.z < NSPACE; coord2.z++) {
									//coord2.z =z;// coord.z;//(z+coord.z)%NSPACE;
									int nspace2 = get_nspace(coord2);
									spinor b_tmp = b[get_pos(nspace2, t2)];
									hmc_complex b_arr[12];
									load_spinor_to_complex_array(b_tmp, b_arr);
									hmc_complex Dminus;
									//loop over spin-colour indices
									for(int a = 0; a < 12; a++) {
										for(int b = 0; b < 12; b++) {
											Dminus =  complexmult(phi_arr[a], complexconj(b_arr[b]));
											correlator += complex_squarenorm(Dminus) * fac;
										}
									}
								}
							}
						}
					}
				}
			}
		}

		out[id_tmp] = correlator;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the pseudoscalar pion correlator in t-direction from pointsources
__kernel void correlator_ps_t(__global hmc_float * const restrict out, __global const spinor * const restrict phi, __global const spinor * const restrict b)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NTIME threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME_LOCAL; id_tmp += global_size) {
		hmc_float correlator = 0.;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp = phi[get_pos(nspace, t)];
					correlator += spinor_squarenorm(tmp);
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[NTIME_OFFSET + id_tmp] += 2. * KAPPA * 2. * KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}


// //this is the scalar correlator in z-direction from pointsources
__kernel void correlator_sc_z(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict b1, __global const spinor * const restrict phi2, __global const spinor * const restrict b2, __global const spinor * const restrict phi3, __global const spinor * const restrict b3, __global const spinor * const restrict phi4, __global const spinor * const restrict b4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_float correlator = 0.;
		uint3 coord;
		coord.z = id_tmp;
		for(int t = 0; t < NTIME_LOCAL; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp;

					tmp = phi1[get_pos(nspace, t)];
					correlator += - su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

					tmp = phi2[get_pos(nspace, t)];
					correlator += - su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

					tmp = phi3[get_pos(nspace, t)];
					correlator += su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

					tmp = phi4[get_pos(nspace, t)];
					correlator += su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME_GLOBAL;
		out[id_tmp] += 2. * KAPPA * 2. * KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the scalar correlator in t-direction from pointsources
__kernel void correlator_sc_t(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict b1, __global const spinor * const restrict phi2, __global const spinor * const restrict b2, __global const spinor * const restrict phi3, __global const spinor * const restrict b3, __global const spinor * const restrict phi4, __global const spinor * const restrict b4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NTIME threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME_LOCAL; id_tmp += global_size) {
		hmc_float correlator = 0.;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp;

					tmp = phi1[get_pos(nspace, t)];
					correlator += - su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

					tmp = phi2[get_pos(nspace, t)];
					correlator += - su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

					tmp = phi3[get_pos(nspace, t)];
					correlator += su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

					tmp = phi4[get_pos(nspace, t)];
					correlator += su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[NTIME_OFFSET + id_tmp] += 2. * KAPPA * 2.*KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//endif

}


// //this is the vector correlator in z-direction from pointsources (gamma4*gamma1)
__kernel void correlator_vx_z(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict b1, __global const spinor * const restrict phi2, __global const spinor * const restrict b2, __global const spinor * const restrict phi3, __global const spinor * const restrict b3, __global const spinor * const restrict phi4, __global const spinor * const restrict b4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		coord.z = id_tmp;

		for(int t = 0; t < NTIME_LOCAL; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					tmp_a = phi1[get_pos(nspace, t)];
					tmp_b = phi2[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re += restmp.re;
					correlator.im += restmp.im;


					tmp_a = phi3[get_pos(nspace, t)];
					tmp_b = phi4[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re += restmp.re;
					correlator.im += restmp.im;
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME_GLOBAL;
		out[id_tmp] += 2. * KAPPA * 2. * KAPPA * 2.*correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the vector correlator in t-direction from pointsources (gamma4*gamma1)
__kernel void correlator_vx_t(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict b1, __global const spinor * const restrict phi2, __global const spinor * const restrict b2, __global const spinor * const restrict phi3, __global const spinor * const restrict b3, __global const spinor * const restrict phi4, __global const spinor * const restrict b4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME_LOCAL; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					tmp_a = phi1[get_pos(nspace, t)];
					tmp_b = phi2[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re += restmp.re;
					correlator.im += restmp.im;


					tmp_a = phi3[get_pos(nspace, t)];
					tmp_b = phi4[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re += restmp.re;
					correlator.im += restmp.im;
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[NTIME_OFFSET + id_tmp] += 2. * KAPPA * 2. * KAPPA * 2. * correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}


// //this is the vector correlator in z-direction from pointsources (gamma4*gamma2)
__kernel void correlator_vy_z(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict b1, __global const spinor * const restrict phi2, __global const spinor * const restrict b2, __global const spinor * const restrict phi3, __global const spinor * const restrict b3, __global const spinor * const restrict phi4, __global const spinor * const restrict b4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		coord.z = id_tmp;
		for(int t = 0; t < NTIME_LOCAL; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					tmp_a = phi1[get_pos(nspace, t)];
					tmp_b = phi2[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;


					tmp_a = phi3[get_pos(nspace, t)];
					tmp_b = phi4[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME_GLOBAL;
		out[id_tmp] += 2. * KAPPA * 2. * KAPPA * 2.*correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the vector correlator in t-direction from pointsources (gamma4*gamma2)
__kernel void correlator_vy_t(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict b1, __global const spinor * const restrict phi2, __global const spinor * const restrict b2, __global const spinor * const restrict phi3, __global const spinor * const restrict b3, __global const spinor * const restrict phi4, __global const spinor * const restrict b4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME_LOCAL; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					tmp_a = phi1[get_pos(nspace, t)];
					tmp_b = phi2[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;


					tmp_a = phi3[get_pos(nspace, t)];
					tmp_b = phi4[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[NTIME_OFFSET + id_tmp] += 2. * KAPPA * 2. * KAPPA * 2. * correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}



// //this is the vector correlator in z-direction from pointsources (gamma4*gamma3)
__kernel void correlator_vz_z(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict b1, __global const spinor * const restrict phi2, __global const spinor * const restrict b2, __global const spinor * const restrict phi3, __global const spinor * const restrict b3, __global const spinor * const restrict phi4, __global const spinor * const restrict b4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_float correlator = 0.0f;
		uint3 coord;
		coord.z = id_tmp;
		for(int t = 0; t < NTIME_LOCAL; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp;

					tmp = phi1[get_pos(nspace, t)];
					correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

					tmp = phi2[get_pos(nspace, t)];
					correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

					tmp = phi3[get_pos(nspace, t)];
					correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

					tmp = phi4[get_pos(nspace, t)];
					correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME_GLOBAL;
		out[id_tmp] += 2. * KAPPA * 2.*KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the vector correlator in t-direction from pointsources (gamma4*gamma3)
__kernel void correlator_vz_t(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict b1, __global const spinor * const restrict phi2, __global const spinor * const restrict b2, __global const spinor * const restrict phi3, __global const spinor * const restrict b3, __global const spinor * const restrict phi4, __global const spinor * const restrict b4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME_LOCAL; id_tmp += global_size) {
		hmc_float correlator = 0.0f;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp;

					tmp = phi1[get_pos(nspace, t)];
					correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

					tmp = phi2[get_pos(nspace, t)];
					correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

					tmp = phi3[get_pos(nspace, t)];
					correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

					tmp = phi4[get_pos(nspace, t)];
					correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[NTIME_OFFSET + id_tmp] += 2. * KAPPA * 2. * KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// axial vector correlators

// //this is the axial vector correlator in z-direction from pointsources (gamma5*gamma4*gamma1)
__kernel void correlator_ax_z(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict b1, __global const spinor * const restrict phi2, __global const spinor * const restrict b2, __global const spinor * const restrict phi3, __global const spinor * const restrict b3, __global const spinor * const restrict phi4, __global const spinor * const restrict b4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		coord.z = id_tmp;
		for(int t = 0; t < NTIME_LOCAL; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					tmp_a = phi1[get_pos(nspace, t)];
					tmp_b = phi2[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;


					tmp_a = phi3[get_pos(nspace, t)];
					tmp_b = phi4[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re += restmp.re;
					correlator.im += restmp.im;
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME_GLOBAL;
		out[id_tmp] += - 2. * KAPPA * 2.*KAPPA * 2.*correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the axial vector correlator in t-direction from pointsources (gamma5*gamma4*gamma1)
__kernel void correlator_ax_t(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict b1, __global const spinor * const restrict phi2, __global const spinor * const restrict b2, __global const spinor * const restrict phi3, __global const spinor * const restrict b3, __global const spinor * const restrict phi4, __global const spinor * const restrict b4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME_LOCAL; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					tmp_a = phi1[get_pos(nspace, t)];
					tmp_b = phi2[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;


					tmp_a = phi3[get_pos(nspace, t)];
					tmp_b = phi4[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re += restmp.re;
					correlator.im += restmp.im;
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[NTIME_OFFSET + id_tmp] += - 2. * KAPPA * 2.*KAPPA * 2. * correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}


// //this is the axial vector correlator in z-direction from pointsources (gamma5*gamma4*gamma2)
__kernel void correlator_ay_z(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict b1, __global const spinor * const restrict phi2, __global const spinor * const restrict b2, __global const spinor * const restrict phi3, __global const spinor * const restrict b3, __global const spinor * const restrict phi4, __global const spinor * const restrict b4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		coord.z = id_tmp;
		for(int t = 0; t < NTIME_LOCAL; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					tmp_a = phi1[get_pos(nspace, t)];
					tmp_b = phi2[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re += restmp.re;
					correlator.im += restmp.im;


					tmp_a = phi3[get_pos(nspace, t)];
					tmp_b = phi4[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME_GLOBAL;
		out[id_tmp] += - 2. * KAPPA * 2.*KAPPA * 2.*correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the axial vector correlator in t-direction from pointsources (gamma5*gamma4*gamma2)
__kernel void correlator_ay_t(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict b1, __global const spinor * const restrict phi2, __global const spinor * const restrict b2, __global const spinor * const restrict phi3, __global const spinor * const restrict b3, __global const spinor * const restrict phi4, __global const spinor * const restrict b4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME_LOCAL; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					tmp_a = phi1[get_pos(nspace, t)];
					tmp_b = phi2[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re += restmp.re;
					correlator.im += restmp.im;


					tmp_a = phi3[get_pos(nspace, t)];
					tmp_b = phi4[get_pos(nspace, t)];

					restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
					correlator.re += restmp.re;
					correlator.im += restmp.im;

					restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
					correlator.re -= restmp.re;
					correlator.im -= restmp.im;
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[NTIME_OFFSET + id_tmp] += - 2. * KAPPA * 2.*KAPPA * 2. * correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}



// //this is the axial vector correlator in z-direction from pointsources (gamma5*gamma4*gamma3)
__kernel void correlator_az_z(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict b1, __global const spinor * const restrict phi2, __global const spinor * const restrict b2, __global const spinor * const restrict phi3, __global const spinor * const restrict b3, __global const spinor * const restrict phi4, __global const spinor * const restrict b4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_float correlator = 0.0f;
		uint3 coord;
		coord.z = id_tmp;
		for(int t = 0; t < NTIME_LOCAL; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp;

					tmp = phi1[get_pos(nspace, t)];
					correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

					tmp = phi2[get_pos(nspace, t)];
					correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

					tmp = phi3[get_pos(nspace, t)];
					correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

					tmp = phi4[get_pos(nspace, t)];
					correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME_GLOBAL;
		out[id_tmp] += - 2. * KAPPA * 2.*KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the axial vector correlator in t-direction from pointsources (gamma5*gamma4*gamma3)
__kernel void correlator_az_t(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict b1, __global const spinor * const restrict phi2, __global const spinor * const restrict b2, __global const spinor * const restrict phi3, __global const spinor * const restrict b3, __global const spinor * const restrict phi4, __global const spinor * const restrict b4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NTIME threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME_LOCAL; id_tmp += global_size) {
		hmc_float correlator = 0.0f;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					spinor tmp;

					tmp = phi1[get_pos(nspace, t)];
					correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

					tmp = phi2[get_pos(nspace, t)];
					correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

					tmp = phi3[get_pos(nspace, t)];
					correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

					tmp = phi4[get_pos(nspace, t)];
					correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[NTIME_OFFSET + id_tmp] += - 2. * KAPPA * 2.* KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}
