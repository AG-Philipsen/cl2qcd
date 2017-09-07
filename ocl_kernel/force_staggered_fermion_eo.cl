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
 * Kernel for the eoprec staggered _@e partial fermion force.
 * @note Observe that the even-odd preconditioning does not directly affect
 *       the gaugemomenta since they are fields existing on the whole lattice as links.
 *       Here eoprec means that we calculate either the force on even sites or on odd ones.
 *       This kernel must then be called twice with different evenodd value to obtain the force
 *       on all sites.
 * 
 * Here, two cases are possible: @n
 *  @arg evenodd = EVEN means that the force will be calculated on even sites  @n
 *  @arg evenodd =  ODD means that the force will be calculated on  odd sites  @n
 * 
 * To make the code understandable, let us summarize the notation here. The force term here calculated
 * is not exactly the total fermion force (this explains the partial adjective), because of the rational
 * approximation. This kernel will be then used to reconstruct the whole force term (see fermion_force.cpp
 * file documentation in physics/algorithms to further informations).
 * In the RHMC the time derivative of the momenta conjugated to links is
 * @code
 *  i*Hdot_\mu(n)=[U_\mu(n)*\sum_{i=1}^k c_i Q^i_\mu(n)]_TA 
 * @endcode
 * where 
 * @code
 *               | +eta_\mu(n) (D_oe X_e^i)_{n+\mu} (X_e^i^\dag)_n     if evenodd = EVEN 
 *  Q^i_\mu(n) = | 
 *               | -eta_\mu(n) (X_e^i)_{n+\mu} ((D_oe X_e^i)^\dag)_n   if evenodd = ODD 
 * 
 * => this kernel will return -i*[U_\mu(n)*Q^i_\mu(n)]_TA
 * @endcode
 *
 * In this kernel only U_\mu(n)*Q^i_\mu(n) is evaluated and not the \sum. In the expression
 * above of Hdot_\mu(n) k is the order of rational approximation and c_i are the numerators.
 * Then for each i from 1 to k we will call this kernel twice. Finally, two different
 * spinorfield are passed to this kernel: A and B. Depending on evenodd they will be different
 * objects, but at this level it does not matter:
 * @code
 *               | +eta_\mu(n) (A)_{n+\mu} (B^\dag)_n     if evenodd = EVEN
 *  Q^i_\mu(n) = |
 *               | -eta_\mu(n) (A)_{n+\mu} (B^\dag)_n     if evenodd = ODD
 * 
 * => return -i*[U_\mu(n)*Q^i_\mu(n)]_TA
 * @endcode
 * 
 * If an imaginary chemical potential is introduced, the time derivative of the momenta conjugated to links
 * modifies due to a prefactor of the sum:
 * @code
 *  i*Hdot_\mu(n)=[U_\mu(n)* (1+\delta_{\mu,4}*(exp(\imath \mu_I)-1)) *\sum_{i=1}^k c_i Q^i_\mu(n)]_TA 
 * @endcode
 * that basically means that we have to multiply Q^i_\mu(n) by exp(\imath \mu_I) if dir==TDIR before 
 * taking the traceless_antihermitian_part. This is important because the factor is in general
 * not real. Then this kernel, if the imaginary chemical potential is used, returns
 * @code
 *  -i*[U_\mu(n)* (1+\delta_{\mu,4}*(exp(\imath \mu_I)-1)) *Q^i_\mu(n)]_TA
 * @endcode
 * 
 * @note Without imaginary chemical potential, Q^i_\mu(n) is evaluated, then the product
 *       U_\mu(n)*Q^i_\mu(n) is calculated and in the end the traceless_antihermitian_part
 *       is taken. If chem_pot_im is different from zero, the multiplication by the factor
 *       (1+\delta_{\mu,4}*(exp(\imath \mu_I)-1)) is added just before the TA operation.
 */

ae fermion_staggered_partial_force_eo_local(__global const Matrixsu3StorageType * const restrict field, __global const staggeredStorageType * const restrict A, __global const staggeredStorageType * const restrict B, const st_index pos, const dir_idx dir)
{
	Matrix3x3 tmp;
	Matrixsu3 U, aux;
	su3vec a, b;
	int n = pos.space;
	int t = pos.time;
	int nn_eo;
	st_index nn;
	//this is used to take into account the staggered phase and the BC-conditions
	hmc_complex eta_mod;
	
	////////////////////////////////////////////////
	nn = get_neighbor_from_st_idx(pos, dir);
	nn_eo = get_eo_site_idx_from_st_idx(nn); //transform normal indices to eoprec index
	U = get_matrixsu3(field, n, t, dir);
	a = get_su3vec_from_field_eo(A, nn_eo);
	
	eta_mod = get_modified_stagg_phase(n, dir);
	a = su3vec_times_complex(a, eta_mod);
	
	b = get_su3vec_from_field_eo(B, get_n_eoprec(n, t));
	tmp = multiply_matrix3x3(matrix_su3to3x3(U),u_times_v_dagger(a, b));
#ifdef _CP_IMAG_
	//Simplest code, if low performance try something like (dir==TDIR)*cpi_tmp to avoid the if.
	//Actually one could also think to include the imaginary chemical potential
	//in the staggered phases as done for the boundary conditions. In this case one
	//should move cpi_tmp to the file operations_staggered.cl
	if(dir == TDIR){
	  hmc_complex cpi_tmp = {COSCPI, SINCPI};
	  tmp = multiply_matrix3x3_by_complex(tmp, cpi_tmp);
	}
#endif
	tmp = traceless_antihermitian_part(tmp);
	aux = matrix_3x3tosu3(multiply_matrix3x3_by_complex(tmp, hmc_complex_minusi));
	////////////////////////////////////////////////
	
	return build_ae_from_su3(aux);  
}

////////////////////////////////////////////////////////////////////
//HERE, BOUNDARY CONDITIONS ARE INCLUDED IN STAGGERED PHASES!!!!! //
////////////////////////////////////////////////////////////////////
__kernel void fermion_staggered_partial_force_eo(__global const Matrixsu3StorageType * const restrict field, __global const staggeredStorageType * const restrict A, __global const staggeredStorageType * const restrict B, __global aeStorageType * const restrict out, int evenodd)
{
	//The following 2 lines were about the Wilson kernel. I do not know if they are still valid.
	// must include HALO, as we are updating neighbouring sites
	// -> not all local sites will fully updated if we don't calculate on halo indices, too
	PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
		//caculate (pos,time) out of id_local depending on evenodd
		st_index pos = (evenodd == EVEN) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);
		
		ae tmp;
		for(dir_idx dir = 0; dir < 4; ++dir) {
			tmp = fermion_staggered_partial_force_eo_local(field, A, B, pos, dir);
			//Depending on evenodd the sign in front of Q^i_\mu(n) is here taken into account
			if(evenodd == EVEN)
			  update_gaugemomentum(tmp, 1., get_link_idx(dir, pos), out);
			else
			  update_gaugemomentum(tmp, -1., get_link_idx(dir, pos), out);
		}
	}
}



