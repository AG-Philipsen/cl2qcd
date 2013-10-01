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
 * is not exactely the total fermion force (this explains the partial adjective), because of the rational
 * approximation. This kernel will be then used to reconstruct the whole force term. In the RHMC the
 * fermion force is
 * @code
 *  i*F_\mu(n)=[U_\mu(n)*\sum_{i=1}^k c_i Q^i_\mu(n)]_TA 
 * @endcode
 * where 
 * @code
 *               | +eta_\mu(n) (D_oe X_e^i)_{n+\mu} (X_e^i^\dag)_n     if evenodd = EVEN 
 *  Q^i_\mu(n) = | 
 *               | -eta_\mu(n) (X_e^i)_{n+\mu} ((D_oe X_e^i)^\dag)_n   if evenodd = ODD 
 * 
 * => return -i*[U_\mu(n)*Q^i_\mu(n)]_TA
 * @endcode
 *
 * In this kernel only U_\mu(n)*Q^i_\mu(n) is evaluated and not the \sum. In the expression
 * above of F_\mu(n) k is the order of rational approximation and c_i are the numerators.
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
 * @todo If a chemical potential is introduced, probably this kernel has to be modified!
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
	
	////////////////////////////////////////////////
	nn = get_neighbor_from_st_idx(pos, dir);
	nn_eo = get_eo_site_idx_from_st_idx(nn); //transform normal indices to eoprec index
	U = get_matrixsu3(field, n, t, dir);
	a = get_su3vec_from_field_eo(A, nn_eo);
	
	coord_spatial coord = get_coord_spatial(n);
	if(coord.x == (NSPACE-1) || coord.y == (NSPACE-1) || coord.z == (NSPACE-1) ||
	   t == (NTIME_GLOBAL-1)){
		hmc_complex eta; //staggered phase + BC
		eta = get_modified_stagg_phase(n, t, dir);
		a = su3vec_times_complex(a, eta);
	}else{
		int eta; //staggered phase
		if(dir != XDIR)
		  eta = get_staggered_phase(n, dir);
		else
		  eta = 1;
		a = su3vec_times_real(a, eta);
	}
	
	b = get_su3vec_from_field_eo(B, get_n_eoprec(n, t));
	tmp = multiply_matrix3x3(matrix_su3to3x3(U),u_times_v_dagger(a, b));
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
			if(evenodd == EVEN)
			  update_gaugemomentum(tmp, 1., get_link_idx(dir, pos), out);
			else
			  update_gaugemomentum(tmp, -1., get_link_idx(dir, pos), out);
		}
	}
}




//The following kernel is the same as that above, except that there are NOT included BC using
//modified staggered phases, i.e. BC are assumed to be periodic in all directions.
//Maybe it is useful for tests to avoid spilling...

/*
__kernel void fermion_staggered_partial_force_eo(__global const Matrixsu3StorageType * const restrict field, __global const staggeredStorageType * const restrict A, __global const staggeredStorageType * const restrict B, __global aeStorageType * const restrict out, int evenodd)
{
	//The following 2 lines were about the Wilson kernel. I do not know if they are still valid.
	// must include HALO, as we are updating neighbouring sites
	// -> not all local sites will fully updated if we don't calculate on halo indices, too
	PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
		//caculate (pos,time) out of id_local depending on evenodd
		st_index pos = (evenodd == EVEN) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);
		
		Matrix3x3 tmp;
		Matrixsu3 U, aux;
		su3vec a, b;
		ae out_tmp;
		int eta; //staggered phase
		int dir;
		int n = pos.space;
		int t = pos.time;
		int nn, nn_eo;

		//go through the different directions. Here we have only positive directions
		//because in the Q^i_\mu(n) we have only {n+\mu}
		///////////////////////////////////
		// mu = +0
		///////////////////////////////////
		dir = TDIR; //here in other parts of the code we have dir=0, but it should be changed
		            //to be coherent with the definitions of geometry. If one wanted to change
		            //conventions, it should be enough to change the operations_geometry.cl file.
		///////////////////////////////////
		nn = get_neighbor_temporal(t);
		nn_eo = get_n_eoprec(n, nn); //transform normal indices to eoprec index
		U = get_matrixsu3(field, n, t, dir);
		a = get_su3vec_from_field_eo(A, nn_eo);
		if(evenodd == EVEN) //Get staggered phase and take into account global sign
		  eta = get_staggered_phase(n, dir);
		else
		  eta = -1 * get_staggered_phase(n, dir);
		a = su3vec_times_real(a, eta);
		b = get_su3vec_from_field_eo(B, get_n_eoprec(n, t));
		tmp = multiply_matrix3x3(matrix_su3to3x3(U),u_times_v_dagger(a, b));
		tmp = traceless_antihermitian_part(tmp);
		aux = matrix_3x3tosu3(multiply_matrix3x3_by_complex(tmp, hmc_complex_minusi));
		out_tmp = build_ae_from_su3(aux);
		//Let's add out_tmp to out in the right site that is get_link_pos(dir, n, t)
		update_gaugemomentum(out_tmp, 1., get_link_pos(dir, n, t), out);

		///////////////////////////////////
		// mu = +1
		///////////////////////////////////
		dir = XDIR; //See comment for dir=TDIR
		///////////////////////////////////
		nn = get_neighbor_spatial(n, dir);
		nn_eo = get_n_eoprec(nn, t); //transform normal indices to eoprec index
		U = get_matrixsu3(field, n, t, dir);
		a = get_su3vec_from_field_eo(A, nn_eo);
		if(evenodd == ODD) //In XDIR the stagg. phase is +1 always, take into account only global sign
			a = su3vec_times_real(a, -1.);
		b = get_su3vec_from_field_eo(B, get_n_eoprec(n, t));
		tmp = multiply_matrix3x3(matrix_su3to3x3(U),u_times_v_dagger(a, b));
		tmp = traceless_antihermitian_part(tmp);
		aux = matrix_3x3tosu3(multiply_matrix3x3_by_complex(tmp, hmc_complex_minusi));
		out_tmp = build_ae_from_su3(aux);
		//Let's add out_tmp to out in the right site that is get_link_pos(dir, n, t)
		update_gaugemomentum(out_tmp, 1., get_link_pos(dir, n, t), out);
		
		///////////////////////////////////
		// mu = +2
		///////////////////////////////////
		dir = YDIR; //See comment for dir=TDIR
		///////////////////////////////////
		nn = get_neighbor_spatial(n, dir);
		nn_eo = get_n_eoprec(nn, t); //transform normal indices to eoprec index
		U = get_matrixsu3(field, n, t, dir);
		a = get_su3vec_from_field_eo(A, nn_eo);
		if(evenodd == EVEN) //Get staggered phase and take into account global sign
		  eta = get_staggered_phase(n, dir);
		else
		  eta = -1 * get_staggered_phase(n, dir);
		a = su3vec_times_real(a, eta);
		b = get_su3vec_from_field_eo(B, get_n_eoprec(n, t));
		tmp = multiply_matrix3x3(matrix_su3to3x3(U),u_times_v_dagger(a, b));
		tmp = traceless_antihermitian_part(tmp);
		aux = matrix_3x3tosu3(multiply_matrix3x3_by_complex(tmp, hmc_complex_minusi));
		out_tmp = build_ae_from_su3(aux);
		//Let's add out_tmp to out in the right site that is get_link_pos(dir, n, t)
		update_gaugemomentum(out_tmp, 1., get_link_pos(dir, n, t), out);
		
		///////////////////////////////////
		// mu = +3
		///////////////////////////////////
		dir = ZDIR; //See comment for dir=TDIR
		///////////////////////////////////
		nn = get_neighbor_spatial(n, dir);
		nn_eo = get_n_eoprec(nn, t); //transform normal indices to eoprec index
		U = get_matrixsu3(field, n, t, dir);
		a = get_su3vec_from_field_eo(A, nn_eo);
		if(evenodd == EVEN) //Get staggered phase and take into account global sign
		  eta = get_staggered_phase(n, dir);
		else
		  eta = -1 * get_staggered_phase(n, dir);
		a = su3vec_times_real(a, eta);
		b = get_su3vec_from_field_eo(B, get_n_eoprec(n, t));
		tmp = multiply_matrix3x3(matrix_su3to3x3(U),u_times_v_dagger(a, b));
		tmp = traceless_antihermitian_part(tmp);
		aux = matrix_3x3tosu3(multiply_matrix3x3_by_complex(tmp, hmc_complex_minusi));
		out_tmp = build_ae_from_su3(aux);
		//Let's add out_tmp to out in the right site that is get_link_pos(dir, n, t)
		update_gaugemomentum(out_tmp, 1., get_link_pos(dir, n, t), out);
	}
}
*/