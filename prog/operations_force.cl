/**
 * @file operations to calculate the forces
 */

hmc_error gauge_force(inputparameters * parameters, hmc_gaugefield * field, hmc_algebraelement2 * out){
	hmc_float beta = (*parameters).get_beta();
	int global_link_pos;
	hmc_3x3matrix V;
	hmc_3x3matrix tmp;
	hmc_su3matrix U;
	
	//Gauge force is factor*Im(i Tr(T_i U V))
	//   with T_i being the SU3-Generator in i-th direction and V the staplematrix
	//   and the factor being -beta/3. (for standard Wilson-action)	
	for(int t = 0; t < NTIME; t++){
		for(int n = 0; n < VOLSPACE; n++){
			for(int mu = 0; mu < NDIM; mu++){
				global_link_pos = get_global_link_pos(mu, n, t);

				calc_staple(field, &V, n, t, mu);
				get_su3matrix(&U, field, n, t, mu);

				/** @TODO CP: this is not valid for REC12 */
				//when not using Reconstruct12 staplematrix and 3x3matrix are just the same!
				multiply_su3matrices (&tmp, &U, &V);

				hmc_algebraelement2 out_tmp;
				tr_lambda_u(tmp, &out_tmp);
				hmc_float factor =-beta/3.;
				update_gaugemomentum(out_tmp, factor, global_link_pos, out);
	}}}
	return HMC_SUCCESS;
}

#ifdef _FERMIONS_

//CP: fermion_force = (gamma_5 Y)^dagger iT_i
//	it is assumed that the results can be added to out!!
hmc_error fermion_force(inputparameters * parameters, hmc_gaugefield * field, hmc_spinor_field * Y, hmc_spinor_field * X, hmc_algebraelement2 * out){
  hmc_su3matrix U, U_up, U_down;
	hmc_3x3matrix v1, v2;
  hmc_su3vector psia,psib,phia,phib;
	hmc_full_spinor y, plus;
	int nn, nup, ndown;
	hmc_algebraelement2 out_tmp;
	int global_link_pos;
	int global_link_pos_down;
	//the additional 2 comes from Tr(lambda_ij) = 2delta_ij
	hmc_float factor = 2.*(*parameters).get_kappa();
	int dir;
	
	//main loop
// 	for(int t = 0; t<NTIME; t++){
// 		for(int n = 0; n<VOLSPACE; n++){
	for(int t = 0; t<NTIME; t++){
		for(int n = 0; n<VOLSPACE; n++){
			get_spinor_from_field(Y, y, n, t);
			///////////////////////////////////
			// Calculate gamma_5 y
			///////////////////////////////////
			gamma_5_spinor(y);

			//go through the different directions
			///////////////////////////////////
			// mu = 0
			///////////////////////////////////
			dir = 0;
			
			// mu = +0
			global_link_pos = get_global_link_pos(dir, n, t);
			nn = (t+1)%NTIME;
			get_spinor_from_field(X, plus, n, nn);
			dir = 1;
			get_su3matrix(&U, field, n, t, dir);
			dir = 0;
		
			//psi = (1-gamma_mu)plus
			spinproj_gamma0_a(plus, psia, -hmc_one_f);
			spinproj_gamma0_b(plus, psib, -hmc_one_f);
			//phi = (1-gamma_mu)y
			spinproj_gamma0_a(y, phia, -hmc_one_f);
			spinproj_gamma0_b(y, phib, -hmc_one_f);

			// v1 = Tr(phi*psi_dagger)
			tr_v_times_u_dagger(phia, psia, phib, psib, v1);

			//U*v1 = U*(phi_a)
			/** @todo what about REC12 and this call??*/
			multiply_3x3matrix (&v2, &U, &v1);
		
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			tr_lambda_u(v2, &out_tmp);
			
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos, out);

			//mu = -0
			nn = (t+NTIME-1)%NTIME;
			global_link_pos_down = get_global_link_pos(dir, n, nn);
			
			get_spinor_from_field(X, plus, n, nn);
			
			get_su3matrix(&U, field, n, nn, dir);
			
			//psi = (1+gamma_mu)plus
			spinproj_gamma0_a(plus, psia, hmc_one_f);
			spinproj_gamma0_b(plus, psib, hmc_one_f);
			//phi = (1+gamma_mu)y
			spinproj_gamma0_a(y, phia, hmc_one_f);
			spinproj_gamma0_b(y, phib, hmc_one_f);

			//CP: here is the difference with regard to +mu-direction: psi and phi interchanged!!
			// v1 = Tr(psi*phi_dagger)
			tr_v_times_u_dagger(psia, phia, psib, phib, v1);

			//U*v1 = U*(phi_a)
			/** @todo what about REC12 and this call??*/
			multiply_3x3matrix (&v2, &U, &v1);
		
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			tr_lambda_u(v2, &out_tmp);
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos_down, out);			
			///////////////////////////////////
			// mu = 1
			///////////////////////////////////
			dir = 1;
			
			// mu = +1
			global_link_pos = get_global_link_pos(dir, n, t);
			nn = get_neighbor(n,t);
			get_spinor_from_field(X, plus, nn, t);
			get_su3matrix(&U, field, n, t, dir);

			//psi = (1-gamma_mu)plus
			spinproj_gamma1_a(plus, psia, -hmc_one_f);
			spinproj_gamma1_b(plus, psib, -hmc_one_f);
			//phi = (1-gamma_mu)y
			spinproj_gamma1_a(y, phia, -hmc_one_f);
			spinproj_gamma1_b(y, phib, -hmc_one_f);

			// v1 = Tr(phi*psi_dagger)
			tr_v_times_u_dagger(phia, psia, phib, psib, v1);

			//U*v1 = U*(phi_a)
			/** @todo what about REC12 and this call??*/
			multiply_3x3matrix (&v2, &U, &v1);
		
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			tr_lambda_u(v2, &out_tmp);
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos, out);

			//mu = -1
			nn = get_lower_neighbor(n,dir);
			global_link_pos_down = get_global_link_pos(dir, nn,t);
			
			get_spinor_from_field(X, plus, nn, t);
			get_su3matrix(&U, field, nn, t, dir);

			//psi = (1+gamma_mu)plus
			spinproj_gamma1_a(plus, psia, hmc_one_f);
			spinproj_gamma1_b(plus, psib, hmc_one_f);
			//phi = (1+gamma_mu)y
			spinproj_gamma1_a(y, phia, hmc_one_f);
			spinproj_gamma1_b(y, phib, hmc_one_f);

			//CP: here is the difference with regard to +mu-direction: psi and phi interchanged!!
			// v1 = Tr(psi*phi_dagger)
			tr_v_times_u_dagger(psia, phia, psib, phib, v1);

			//U*v1 = U*(phi_a)
			/** @todo what about REC12 and this call??*/
			multiply_3x3matrix (&v2, &U, &v1);
		
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			tr_lambda_u(v2, &out_tmp);
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos_down, out);
			///////////////////////////////////
			// mu = 2
			///////////////////////////////////
			dir = 2;
			
			// mu = +2
			global_link_pos = get_global_link_pos(dir, n, t);
			nn = get_neighbor(n,t);
			get_spinor_from_field(X, plus, nn, t);
			get_su3matrix(&U, field, n, t, dir);
		
			//psi = (1-gamma_mu)plus
			spinproj_gamma2_a(plus, psia, -hmc_one_f);
			spinproj_gamma2_b(plus, psib, -hmc_one_f);
			//phi = (1-gamma_mu)y
			spinproj_gamma2_a(y, phia, -hmc_one_f);
			spinproj_gamma2_b(y, phib, -hmc_one_f);

			// v1 = Tr(phi*psi_dagger)
			tr_v_times_u_dagger(phia, psia, phib, psib, v1);

			//U*v1 = U*(phi_a)
			/** @todo what about REC12 and this call??*/
			multiply_3x3matrix (&v2, &U, &v1);
		
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			tr_lambda_u(v2, &out_tmp);
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos, out);

			//mu = -2
			nn = get_lower_neighbor(n,dir);
			global_link_pos_down = get_global_link_pos(dir, nn,t);
			
			get_spinor_from_field(X, plus, nn, t);
			get_su3matrix(&U, field, nn, t, dir);
		
			//psi = (1+gamma_mu)plus
			spinproj_gamma2_a(plus, psia, hmc_one_f);
			spinproj_gamma2_b(plus, psib, hmc_one_f);
			//phi = (1+gamma_mu)y
			spinproj_gamma2_a(y, phia, hmc_one_f);
			spinproj_gamma2_b(y, phib, hmc_one_f);

			//CP: here is the difference with regard to +mu-direction: psi and phi interchanged!!
			// v1 = Tr(psi*phi_dagger)
			tr_v_times_u_dagger(psia, phia, psib, phib, v1);

			//U*v1 = U*(phi_a)
			/** @todo what about REC12 and this call??*/
			multiply_3x3matrix (&v2, &U, &v1);
		
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			tr_lambda_u(v2, &out_tmp);
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos_down, out);
			///////////////////////////////////
			// mu = 3
			///////////////////////////////////
			dir = 3;
			
			// mu = +3
			global_link_pos = get_global_link_pos(dir, n, t);
			nn = get_neighbor(n,t);
			get_spinor_from_field(X, plus, nn, t);
			get_su3matrix(&U, field, n, t, dir);
		
			//psi = (1-gamma_mu)plus
			spinproj_gamma3_a(plus, psia, -hmc_one_f);
			spinproj_gamma3_b(plus, psib, -hmc_one_f);
			//phi = (1-gamma_mu)y
			spinproj_gamma3_a(y, phia, -hmc_one_f);
			spinproj_gamma3_b(y, phib, -hmc_one_f);

			// v1 = Tr(phi*psi_dagger)
			tr_v_times_u_dagger(phia, psia, phib, psib, v1);

			//U*v1 = U*(phi_a)
			/** @todo what about REC12 and this call??*/
			multiply_3x3matrix (&v2, &U, &v1);
		
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			tr_lambda_u(v2, &out_tmp);
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos, out);

			//mu = -3
			nn = get_lower_neighbor(n,dir);
			global_link_pos_down = get_global_link_pos(dir, nn,t);
			
			get_spinor_from_field(X, plus, nn, t);
			get_su3matrix(&U, field, nn, t, dir);
		
			//psi = (1+gamma_mu)plus
			spinproj_gamma3_a(plus, psia, hmc_one_f);
			spinproj_gamma3_b(plus, psib, hmc_one_f);
			//phi = (1+gamma_mu)y
			spinproj_gamma3_a(y, phia, hmc_one_f);
			spinproj_gamma3_b(y, phib, hmc_one_f);

			//CP: here is the difference with regard to +mu-direction: psi and phi interchanged!!
			// v1 = Tr(psi*phi_dagger)
			tr_v_times_u_dagger(psia, phia, psib, phib, v1);

			//U*v1 = U*(phi_a)
			/** @todo what about REC12 and this call??*/
			multiply_3x3matrix (&v2, &U, &v1);
		
			//TODO
	 		//ka0 is kappa*BC-factor
			//     _complex_times_su3(v1,ka0,v2);
			//this must become v1 if the above is included again
			tr_lambda_u(v2, &out_tmp);
			//what is the factor here??
			update_gaugemomentum(out_tmp, factor, global_link_pos_down, out);
		}}
	return HMC_SUCCESS;
}

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
