/**
 * @file operations to calculate the forces
 */

__kernel void gauge_force(__global ocl_s_gaugefield * field, __global ae * out){
	
	int id = get_global_id(0);
	if(id == 0){
	
	
	hmc_float factor =-BETA/3.;
	int global_link_pos;
	Matrix3x3 V, tmp, tmp2;
	Matrixsu3 U;
	ae out_tmp;
	
	//Gauge force is factor*Im(i Tr(T_i U V))
	//   with T_i being the SU3-Generator in i-th direction and V the staplematrix
	//   and the factor being -beta/3. (for standard Wilson-action)	
	for(int t = 0; t < NTIME; t++){
		for(int n = 0; n < VOLSPACE; n++){
			for(int mu = 0; mu < NDIM; mu++){
				global_link_pos = get_global_link_pos(mu, n, t);
				V = calc_staple(field, n, t, mu);
				U = get_matrixsu3(field, n, t, mu);

				tmp2 = matrix_su3to3x3(U);
				tmp = multiply_matrix3x3 (tmp2, V);

				out_tmp = tr_lambda_u(tmp);

				update_gaugemomentum(out_tmp, factor, global_link_pos, out);
	}}}
	
	}
}

#ifdef _FERMIONS_

//CP: fermion_force = (gamma_5 Y)^dagger iT_i
//	it is assumed that the results can be added to out!!
__kernel void fermion_force(__global ocl_s_gaugefield * field,__global  spinorfield * Y,__global  spinorfield * X,__global  ae * out){
  Matrixsu3 U;
	Matrix3x3 v1, v2, tmp;
  su3vec psia,psib,phia,phib;
	spinor y, plus;
	int nn, nup, ndown;
	ae out_tmp;
	int global_link_pos;
	int global_link_pos_down;
	//the 2 here comes from Tr(lambda_ij) = 2delta_ij
	hmc_float factor = 2.*KAPPA;
	int dir;
	
	//main loop
	for(int t = 0; t<NTIME; t++){
		for(int n = 0; n<VOLSPACE; n++){
			y = get_spinor_from_field(Y, n, t);
			///////////////////////////////////
			// Calculate gamma_5 y
			///////////////////////////////////
			y = gamma5_local(y);

			//go through the different directions
			///////////////////////////////////
			// mu = 0
			///////////////////////////////////
			dir = 0;
			
			///////////////////////////////////
			//mu = +0
			global_link_pos = get_global_link_pos(dir, n, t);
			nn = (t+1)%NTIME;
			plus = get_spinor_from_field(X, n, nn);
			U = get_matrixsu3(field, n, t, dir);
			///////////////////////////////////
			// Calculate psi/phi = (1 - gamma_0) plus/y
			// with 1 - gamma_0:
			// | 1  0 -1  0 |        |       psi.e0 - psi.e2 |
			// | 0  1  0 -1 |  psi = |       psi.e1 - psi.e3 |
			// |-1  0  1  0 |        |(-1)*( psi.e1 - psi.e3 |
			// | 0 -1  0  1 |        |(-1)*( psi.e0 - psi.e2 |
			///////////////////////////////////
			//here, only the independent components are needed
			//the other ones are incorporated in the dirac trace
			// psia = 0. component of (1-gamma_0)plus
			psia= su3vec_dim(plus.e0, plus.e2);
			// psib = 1. component of (1-gamma_0)plus
			psib = su3vec_dim(plus.e1, plus.e3);
			phia= su3vec_dim(y.e0, y.e2);
			phib = su3vec_dim(y.e1, y.e3);
			// v1 = Tr(phi*psi_dagger)
			v1 = tr_v_times_u_dagger(phia, psia, phib, psib);
			//U*v1 = U*(phi_a)
			tmp = matrix_su3to3x3(U);
			v2 = multiply_matrix3x3 (tmp, v1);
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			out_tmp = tr_lambda_u(v2);
			update_gaugemomentum(out_tmp, factor, global_link_pos, out);
			
			/////////////////////////////////////
			//mu = -0
			nn = (t+NTIME-1)%NTIME;
			global_link_pos_down = get_global_link_pos(dir, n, nn);
			plus = get_spinor_from_field(X, n, nn);
			U = get_matrixsu3(field, n, nn, dir);
			///////////////////////////////////
			// Calculate psi/phi = (1 + gamma_0) y
			// with 1 + gamma_0:
			// | 1  0  1  0 |       | psi.e0 + i*psi.e2 |
			// | 0  1  0  1 | psi = | psi.e1 + i*psi.e3 |
			// | 1  0  1  0 |       | psi.e1 + i*psi.e2 |
			// | 0  1  0  1 |       | psi.e0 + i*psi.e3 |
			///////////////////////////////////
			psia = su3vec_acc(plus.e0, plus.e2);
			psib = su3vec_acc(plus.e1, plus.e3);
			phia = su3vec_acc(y.e0, y.e2);
			phib = su3vec_acc(y.e1, y.e3);
			//CP: here is the difference with regard to +mu-direction: psi and phi interchanged!!
			// v1 = Tr(psi*phi_dagger)
			v1 = tr_v_times_u_dagger(psia, phia, psib, phib);
			//U*v1 = U*(phi_a)
			tmp = matrix_su3to3x3(U);
			v2 = multiply_matrix3x3 (tmp, v1);
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			out_tmp = tr_lambda_u(v2);
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos_down, out);
			
			//comments correspond to the mu=0 ones
			///////////////////////////////////
			// mu = 1
			///////////////////////////////////
			dir = 1;
			
			///////////////////////////////////
			// mu = +1
			global_link_pos = get_global_link_pos(dir, n, t);
			nn = get_neighbor(n,t);
			plus = get_spinor_from_field(X, nn, t);
			U = get_matrixsu3(field, n, t, dir);
			///////////////////////////////////
			// Calculate (1 - gamma_1) y
			// with 1 - gamma_1:
			// | 1  0  0 -i |       |      psi.e0 - i*psi.e3  |
			// | 0  1 -i  0 | psi = |      psi.e1 - i*psi.e2  |
			// | 0  i  1  0 |       |(i)*( psi.e1 - i*psi.e2) |
			// | i  0  0  1 |       |(i)*( psi.e0 - i*psi.e3) |
			///////////////////////////////////
			//psi = (1-gamma_mu)plus
			psia = su3vec_dim_i(plus.e0, plus.e3);
			psib = su3vec_dim_i(plus.e1, plus.e2);
			phia = su3vec_dim_i(y.e0, y.e3);
			phib = su3vec_dim_i(y.e1, y.e2);


			v1 = tr_v_times_u_dagger(psia, phia, psib, phib);
			tmp = matrix_su3to3x3(U);
			v2 = multiply_matrix3x3 (tmp, v1);
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			out_tmp = tr_lambda_u(v2);
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos, out);
			///////////////////////////////////
			//mu = -1
			nn = get_lower_neighbor(n,dir);
			global_link_pos_down = get_global_link_pos(dir, nn,t);
			plus = get_spinor_from_field(X, nn, t);
			U = get_matrixsu3(field, nn, t, dir);
			///////////////////////////////////
			// Calculate (1 + gamma_1) y
			// with 1 + gamma_1:
			// | 1  0  0  i |       |       psi.e0 + i*psi.e3  |
			// | 0  1  i  0 | psi = |       psi.e1 + i*psi.e2  |
			// | 0 -i  1  0 |       |(-i)*( psi.e1 + i*psi.e2) |
			// |-i  0  0  1 |       |(-i)*( psi.e0 + i*psi.e3) |
			///////////////////////////////////
			psia = su3vec_acc_i(plus.e0, plus.e3);
			psib = su3vec_acc_i(plus.e1, plus.e2);
			phia = su3vec_acc_i(y.e0, y.e3);
			phib = su3vec_acc_i(y.e1, y.e2);
			
			v1 = tr_v_times_u_dagger(psia, phia, psib, phib);
			tmp = matrix_su3to3x3(U);
			v2 = multiply_matrix3x3 (tmp, v1);
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			out_tmp = tr_lambda_u(v2);
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos_down, out);
			
			///////////////////////////////////
			// mu = 2
			///////////////////////////////////
			dir = 2;
			
			///////////////////////////////////
			// mu = +2
			global_link_pos = get_global_link_pos(dir, n, t);
			nn = get_neighbor(n,t);
			plus = get_spinor_from_field(X, nn, t);
			U = get_matrixsu3(field, n, t, dir);
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
			
			v1 = tr_v_times_u_dagger(psia, phia, psib, phib);
			tmp = matrix_su3to3x3(U);
			v2 = multiply_matrix3x3 (tmp, v1);
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			out_tmp = tr_lambda_u(v2);
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos, out);
			
			///////////////////////////////////
			//mu = -2
			nn = get_lower_neighbor(n,dir);
			global_link_pos_down = get_global_link_pos(dir, nn,t);
			plus = get_spinor_from_field(X, nn, t);
			U = get_matrixsu3(field, nn, t, dir);
		
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
			
			v1 = tr_v_times_u_dagger(psia, phia, psib, phib);
			tmp = matrix_su3to3x3(U);
			v2 = multiply_matrix3x3 (tmp, v1);
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			out_tmp = tr_lambda_u(v2);
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos_down, out);
			
			///////////////////////////////////
			// mu = 3
			///////////////////////////////////
			dir = 3;
			
			///////////////////////////////////
			// mu = +3
			global_link_pos = get_global_link_pos(dir, n, t);
			nn = get_neighbor(n,t);
			plus = get_spinor_from_field(X, nn, t);
			U = get_matrixsu3(field, n, t, dir);

			///////////////////////////////////
			// Calculate (1 - gamma_3) y
			// with 1 - gamma_3:
			// | 1  0 -i  0 |        |       psi.e0 - i*psi.e2  |
			// | 0  1  0  i |  psi = |       psi.e1 + i*psi.e3  |
			// | i  0  1  0 |        |   i *(psi.e0 - i*psi.e2) |
			// | 0 -i  0  1 |        | (-i)*(psi.e1 + i*psi.e3) |
			///////////////////////////////////
			psia = su3vec_dim_i(plus.e0, plus.e2);
			psib = su3vec_acc_i(plus.e1, plus.e3);
			phia = su3vec_dim_i(y.e0, y.e2);
			phib = su3vec_acc_i(y.e1, y.e3);
			
			v1 = tr_v_times_u_dagger(psia, phia, psib, phib);
			tmp = matrix_su3to3x3(U);
			v2 = multiply_matrix3x3 (tmp, v1);
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			out_tmp = tr_lambda_u(v2);
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos, out);
			
			///////////////////////////////////
			//mu = -3
			nn = get_lower_neighbor(n,dir);
			global_link_pos_down = get_global_link_pos(dir, nn,t);
			plus = get_spinor_from_field(X, nn, t);
			U = get_matrixsu3(field, nn, t, dir);

			///////////////////////////////////
			// Calculate (1 + gamma_3) y
			// with 1 + gamma_3:
			// | 1  0  i  0 |       |       psi.e0 + i*psi.e2  |
			// | 0  1  0 -i | psi = |       psi.e1 - i*psi.e3  |
			// |-i  0  1  0 |       | (-i)*(psi.e0 + i*psi.e2) |
			// | 0  i  0  1 |       |   i *(psi.e1 - i*psi.e3) |
			///////////////////////////////////
			psia = su3vec_acc_i(plus.e0, plus.e2);
			psib = su3vec_dim_i(plus.e1, plus.e3);
			phia = su3vec_acc_i(y.e0, y.e2);
			phib = su3vec_dim_i(y.e1, y.e3);
			
			v1 = tr_v_times_u_dagger(psia, phia, psib, phib);
			tmp = matrix_su3to3x3(U);
			v2 = multiply_matrix3x3 (tmp, v1);
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			out_tmp = tr_lambda_u(v2);
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos_down, out);
			
		}}
		
}
#endif
			
			
//this has to go into a wrapper function!!
/*
//CP: this essentially calculates a hmc_gauge_momentum vector
//CP: if fermions are used, here is the point where the inversion has to be performed
hmc_error force(inputparameters * parameters, hmc_gaugefield * field
	#ifdef _FERMIONS_
	, hmc_spinor_field * phi, hmc_spinor_field * phi_inv
	#endif
	, hmc_algebraelement2 * out){
	cout << "\t\tstart calculating the force..." << endl;
	//CP: make sure that the output field is set to zero
	set_zero_gaugemomenta(out);
	//add contributions
	cout << "\t\tcalc gauge_force..." << endl;
	gauge_force(parameters, field, out);
#ifdef _FERMIONS_
	cout << "\t\tinvert fermion field..." << endl;
	//the algorithm needs two spinor-fields
	hmc_spinor_field* X = new hmc_spinor_field[SPINORFIELDSIZE];
	//CP: to begin with, consider only the cg-solver
	//source is at 0
	int k = 0;
	int use_cg = TRUE;
	//CP: at the moment, use_eo = 0 so that even-odd is not used!!!!!
	
	//debugging
	int err = 0;
	
	if(use_cg){
		if(!use_eo){
			//the inversion calculates Y = (QplusQminus)^-1 phi = phi_inv
			hmc_spinor_field b[SPINORFIELDSIZE];
			create_point_source(parameters,k,0,0,b);
			cout << "\t\t\tstart solver" << endl;
			err = solver(parameters, phi, b, field, use_cg, phi_inv);
			if (err != HMC_SUCCESS) cout << "\t\tsolver did not solve!!" << endl;
			else cout << "\t\tsolver solved!" << endl;
		}
		else{
			hmc_eoprec_spinor_field be[EOPREC_SPINORFIELDSIZE];
			hmc_eoprec_spinor_field bo[EOPREC_SPINORFIELDSIZE];
			
			create_point_source_eoprec(parameters, k,0,0, field, be,bo);
			solver_eoprec(parameters, phi, be, bo, field, use_cg, phi_inv);
		}
		cout << "\t\t\tcalc X" << endl;
		//X = Qminus Y = Qminus phi_inv 
		Qminus(parameters, phi_inv, field, X);
	}
	else{
		//here, one has first to invert (with BiCGStab) Qplus phi = X and then invert Qminus X => Qminus^-1 Qplus^-1 phi = (QplusQminus)^-1 phi = Y = phi_inv
	}
/** @todo control the fields here again!!! 
	fermion_force(parameters, field, phi_inv, X, out);
	
	delete [] X;
#endif
	return HMC_SUCCESS;
}
*/
