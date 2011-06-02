#include "host_hmc.h"

using namespace std;


void convert_ae_to_ae2(hmc_algebraelement in, hmc_algebraelement2 * out){
	(*out).e0 = in[0];
	(*out).e1 = in[1];
	(*out).e2 = in[2];
	(*out).e3 = in[3];
	(*out).e4 = in[4];
	(*out).e5 = in[5];
	(*out).e6 = in[6];
	(*out).e7 = in[7];
}

void convert_ae2_to_ae(hmc_algebraelement2 in, hmc_algebraelement out){
	out[0] = (in).e0;
	out[1] = (in).e1;
	out[2] = (in).e2;
	out[3] = (in).e3;
	out[4] = (in).e4;
	out[5] = (in).e5;
	out[6] = (in).e6;
	out[7] = (in).e7;
}

void convert_ae2_to_ae_global(hmc_algebraelement2 * in, hmc_gauge_momentum * out){
	for(int i = 0; i<GAUGEMOMENTASIZE2; i++){
		convert_ae2_to_ae(in[i], &(out[i*8]));
	}
}

void convert_ae_to_ae2_global(hmc_gauge_momentum * in, hmc_algebraelement2 * out){
	for(int i = 0; i<GAUGEMOMENTASIZE2; i++){
		convert_ae_to_ae2(&(in[i*8]), &out[i]);
	}
}

void update_gaugemomentum(hmc_algebraelement2 in, hmc_float factor, int global_link_pos, hmc_algebraelement2 * out){
			out[global_link_pos].e0 += factor*in.e0;
			out[global_link_pos].e1 += factor*in.e1;
			out[global_link_pos].e2 += factor*in.e2;
			out[global_link_pos].e3 += factor*in.e3;
			out[global_link_pos].e4 += factor*in.e4;
			out[global_link_pos].e5 += factor*in.e5;
			out[global_link_pos].e6 += factor*in.e6;
			out[global_link_pos].e7 += factor*in.e7;
}

void acc_factor_times_algebraelement(hmc_algebraelement2 * inout, hmc_float factor, hmc_algebraelement2 force_in){
	(*inout).e1+=factor*(force_in).e1;
	(*inout).e2+=factor*(force_in).e2;
	(*inout).e3+=factor*(force_in).e3;
	(*inout).e4+=factor*(force_in).e4;
	(*inout).e5+=factor*(force_in).e5;
	(*inout).e6+=factor*(force_in).e6;
	(*inout).e7+=factor*(force_in).e7;
	(*inout).e8+=factor*(force_in).e8; 
}

//CP: molecular dynamics update for the gauge momenta:
//p_out = p_in - eps/2 force(u_in, phi)
//it is assumed that the force term has already been computed. then one only has real-vectors and this is essentially adding one vector to another...
hmc_error md_update_gauge_momenta(hmc_float eps, hmc_algebraelement2 * p_inout, hmc_algebraelement2 * force_in){
	for(int i = 0; i<GAUGEMOMENTASIZE2; i++){
		acc_factor_times_algebraelement(&p_inout[i], -1.*eps, force_in[i]);
	}
	return HMC_SUCCESS;
}

//molecular dynamics update for the gaugefield:
//u_out = exp(i eps p_in) u_in
hmc_error md_update_gaugefield(hmc_float eps, hmc_gauge_momentum * p_in, hmc_gaugefield * u_inout){
	int index;
	for(int t = 0; t<NTIME; t++){
		for(int pos = 0; pos < VOLSPACE; pos++){
			for(int mu = 0; mu<NDIM; mu++){
			hmc_su3matrix tmp;
			hmc_su3matrix tmp2;
			index= get_global_link_pos(mu, pos, t);
			// an su3 algebra element has NC*NC-1 = 8 hmc_float entries
			// &(p_in[index*8]) should point to the right position for the pos-th element of the long gaugemomentum vector p_in
			build_su3matrix_by_exponentiation(&(p_in[index*8]), &tmp2, eps);
			get_su3matrix(&tmp, u_inout, pos, t, mu);
			accumulate_su3matrix_prod( &tmp2, &tmp);
			put_su3matrix(u_inout, &tmp2, pos, t, mu);
	}}}
	return HMC_SUCCESS;
}

#ifdef _FERMIONS_
//phi = Q+ chi
hmc_error md_update_spinorfield(hmc_spinor_field * in, hmc_spinor_field * out, hmc_gaugefield * field, inputparameters * parameters){
	Qplus(parameters, in, field, out);
	return HMC_SUCCESS;
}
#endif

// beta * sum_links sum_nu>mu ( 3 - Tr Re Plaquette )
//CP: since one is only interested in differences of s_gauge, the constant part can be left out!!
hmc_float s_gauge(hmc_gaugefield * field, hmc_float beta){
	/** @TODO CP: implement saving of plaquette measurement (and possibly t_plaq and s_plaq and also polyakov-loop??) */
	hmc_float plaq=0;
	//CP: alternative method: use already existing plaquette-functions
	hmc_float t_plaq;
	hmc_float s_plaq;
	plaq = plaquette(field, &t_plaq, &s_plaq);
	//plaq is normalized by factor of 2.0/(VOL4D*NDIM*(NDIM-1)*NC), so one has to divide by it again
	hmc_float factor = 2.0/static_cast<hmc_float>(VOL4D*NDIM*(NDIM-1)*NC);
	return beta/6./factor*( - plaq);
}

#ifdef _FERMIONS_
// sum_links phi*_i (M^+M)_ij^-1 phi_j
// here it is assumed that the rhs has already been computed... (because usually it has been during the leapfrog..)
hmc_complex s_fermion(hmc_spinor_field * phi, hmc_spinor_field * QplusQminusphi){
	return scalar_product(phi, QplusQminusphi);
}
#endif /* _FERMIONS_ */

//S_gauge + S_fermion + S_gaugemomenta
hmc_complex hamiltonian(hmc_gaugefield * field, hmc_float beta, hmc_gauge_momentum * p
		#ifdef _FERMIONS_
		, hmc_spinor_field * phi, hmc_spinor_field * phi_inv
		#endif
		){
	hmc_complex result;
	hmc_complex tmp;
	(result) = {0.,0.};
	(result).re += s_gauge(field, beta);
	#ifdef _FERMIONS_
	tmp  = s_fermion(phi, phi_inv);
	cout << "s_fermion is:\t" << tmp.re << "  " << tmp.im << endl;
	complexaccumulate(&result, &tmp);
	#endif
	//s_gm = 1/2*squarenorm(Pl)
	hmc_float s_gm;
	gaugemomenta_squarenorm(p, &s_gm);
	result.re += 0.5*s_gm;
	
	return result;
}

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
	hmc_float factor = (*parameters).get_kappa();
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
#endif

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
// 	gauge_force(parameters, field, out);
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
	hmc_complex tmp;
tmp  = s_fermion(phi, phi_inv);
cout << "s_fermion of inverted phi is: " << tmp.re << " " << tmp.im << endl;
/** @todo control the fields here again!!! */
	fermion_force(parameters, field, phi_inv, X, out);
	
	delete [] X;
#endif
	return HMC_SUCCESS;
}

hmc_error metropolis(hmc_float rndnumber, hmc_float beta, 
										 #ifdef _FERMIONS_
										 hmc_spinor_field * phi, hmc_spinor_field * phi_inv, hmc_spinor_field * phi_inv_orig, 
										 #endif
										 hmc_gaugefield * field,	hmc_gauge_momentum * p, hmc_gaugefield * new_field, hmc_gauge_momentum * new_p){
	// takes:
	//		phi and beta as constant
	//		new/old versions of gaugefield and of momenta
	//		and a random number
	//		if it has to be, performs the change old->new, and returns true if there are no failures.
	hmc_complex h_old = hamiltonian(field, beta, p
																	#ifdef _FERMIONS_
																	, phi, phi_inv_orig
																	#endif
																	);
	hmc_complex h_new = hamiltonian(new_field, beta, new_p
																	#ifdef _FERMIONS_
																	, phi, phi_inv
																	#endif
																	);
																	
																	
	if(h_old.im > projectioneps){
		printf("\n\tError: imaginary part in H_OLD [in function: metropolis(...)].\n");
		return HMC_COMPLEX_HAMILTONIANERROR;
	}
	if(h_new.im > projectioneps){
		printf("\n\tError: imaginary part in H_NEW [in function: metropolis(...)].\n");
		return HMC_COMPLEX_HAMILTONIANERROR;
	}
	/** @todo CP:  export h_diff */
	hmc_float h_diff = h_old.re - h_new.re;
	hmc_float compare_prob;
	if(h_diff<0){
		compare_prob = exp(h_diff);
	}else{
		compare_prob = 1.0;
	}
	// SL: the following can be tuned, whether it is more costly to draw always the rnd number even when compare_prob=1
	//     and whether the "if compare_prob==1" costs more or less than always evaluating the exp ...
	if(rndnumber <= compare_prob){
		// perform the change nonprimed->primed !
		copy_gaugefield(new_field, field);
		copy_gaugemomenta(new_p, p);
		// SL: this works as long as p and field are pointers to the *original* memory locations!
		cout << "new configuration accepted" << endl;
	}
	else{
		cout << "new configuration rejected" << endl;
	}
	return HMC_SUCCESS;
}

//it is assumed that gaugefield and gaugemomentum have been set to the old ones already
hmc_error leapfrog(inputparameters * parameters, 
									 #ifdef _FERMIONS_
									 hmc_spinor_field * phi, hmc_spinor_field * phi_inv, hmc_spinor_field * phi_inv_orig, 
									 #endif
									 hmc_gaugefield * u_out, hmc_gauge_momentum * p_out	){
	// CP: it operates directly on the fields p_out and u_out
	int steps = (*parameters).get_integrationsteps1() ;	
	hmc_float stepsize = ((*parameters).get_tau()) /((hmc_float) steps);
	int k;
	hmc_float stepsize_half = 0.5*stepsize;
	
	hmc_gauge_momentum* force_tmp = new hmc_gauge_momentum[GAUGEMOMENTASIZE];
	hmc_algebraelement2* force_tmp2 = new hmc_algebraelement2[GAUGEMOMENTASIZE2];
	hmc_algebraelement2* p2 = new hmc_algebraelement2[GAUGEMOMENTASIZE2];

	//initial step
	cout << "\tinitial step:" << endl;
	cout << p2[0].e2 << endl;
	convert_ae_to_ae2_global( p_out, p2);
	cout << p2[0].e2 << endl;
	//here, phi is inverted using the orig. gaugefield
	force(parameters, u_out ,
		#ifdef _FERMIONS_
		phi, phi_inv_orig, 
		#endif
		force_tmp2);
	md_update_gauge_momenta(stepsize_half, p2, force_tmp2);
	
	//intermediate steps
	if(steps > 1) cout << "\tperform " << steps << " intermediate steps " << endl;
	for(k = 1; k<steps; k++){
		convert_ae2_to_ae_global(p2, p_out);
		md_update_gaugefield(stepsize, p_out, u_out);
		force(parameters, u_out ,
			#ifdef _FERMIONS_
			phi, phi_inv, 
			#endif
			force_tmp2);
		md_update_gauge_momenta(stepsize, p2, force_tmp2);
	}
	
	//final step
	cout << "\tfinal step" << endl;
	convert_ae2_to_ae_global(p2, p_out);
	md_update_gaugefield(stepsize, p_out, u_out);
	force(parameters, u_out ,
		#ifdef _FERMIONS_
		phi, phi_inv, 
		#endif
		force_tmp2);
	md_update_gauge_momenta(stepsize_half, p2,force_tmp2); 
	convert_ae2_to_ae_global(p2, p_out);
	
	delete [] force_tmp;
	
	cout << "\tfinished leapfrog" << endl;
	return HMC_SUCCESS;
}

