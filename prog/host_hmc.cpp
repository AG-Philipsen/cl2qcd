#include "host_hmc.h"

//TODO CP: all return values have to be revisited, since they should be real numbers and not complex one in the end

//TODO in the end hmc_gaugemomenta has to be a vector of NC*NC matrices!! This means that also the declaration of such variables has to be revisited



//TODO the following three functions should later go into a new, seperate file

//molecular dynamics update for the gauge momenta:
//p_out = p_in - eps/2 force(u_in, phi)
//it is assumed that the force term has already been computed. then one only has real-vectors
hmc_error md_update_gauge_momenta(hmc_float eps, hmc_gauge_momentum * p_in, hmc_gauge_momentum * force_in, hmc_gauge_momentum * p_out){
	for(int i = 0; i<GAUGEMOMENTASIZE; i++){
		p_out[i] -= eps*force_in[i];
	}
	return HMC_SUCCESS;
}

//molecular dynamics update for the gaugefield:
//u_out = exp(i eps p_in) u_in
hmc_error md_update_gaugefield(hmc_float eps, hmc_gauge_momentum * p_in, hmc_gaugefield * u_in){
	for(int t = 0; t<NTIME; t++){
		for(int pos = 0; pos < VOLSPACE; pos++){
			for(int mu = 0; mu<NDIM; mu++){
			hmc_su3matrix tmp;
			hmc_su3matrix tmp2;
			//TODO CP: here one has to calculate exp (i eps p), which is in principle again a traceless, hermitian matrix. Let it be stored in tmp2, so one has to ensure that it fits the su3-format
			get_su3matrix(&tmp, u_in,pos, t, mu);
		
			//TODO check the order of multiplication again!! tmp U or Utmp?? does it matter??
			accumulate_su3matrix_prod( &tmp, &tmp2);

			put_su3matrix(u_in, &tmp, pos, t, mu);
	}}}
		
	return HMC_SUCCESS;
}

#ifdef _FERMIONS_
//phi = D chi
hmc_error md_update_spinorfield(hmc_spinor_field * in, hmc_spinor_field * out, hmc_gaugefield * field, inputparameters * parameters){
	//TODO extract needed parameters from paramters, or transform this into M itself
	hmc_float kappa, mu, theta, chem_pot_re, chem_pot_im;
	//TODO check again if it is M or Mdagger here
	Mdagger(in, out, field, kappa, mu, theta, chem_pot_re, chem_pot_im);
	
	return HMC_SUCCESS;
}
#endif

//TODO chech definition of beta again, is it 2/g^2??
// beta * sum_links sum_nu>mu ( 3 - Tr Re Plaquette )
hmc_float s_gauge(hmc_gaugefield * field, hmc_float beta){
	//TODO implement saving of plaquette measurement (and possibly t_plaq and s_plaq and also polyakov-loop??)
	hmc_float plaq=0;

// 	for(int t=0;t<NTIME;t++) {
// 		for(int n=0;n<VOLSPACE;n++) {
// 			for(int mu=0; mu<NDIM; mu++) {
// 				for(int nu=0;nu<mu; nu++) {
// 					hmc_su3matrix prod;
// 					local_plaquette(field, &prod, n, t, mu, nu );
// 					hmc_float tmpfloat = 3. - trace_su3matrix(&prod).re;
// 					plaq += tmpfloat;
// 				}}}}
// 	//normalize
// 	return beta*plaq*2.0/static_cast<hmc_float>(VOL4D*NDIM*(NDIM-1)*NC);
	
	//CP: alternative method: use already existing plaquette-functions
	hmc_float t_plaq;
	hmc_float s_plaq;
	
	plaq = plaquette(field, &t_plaq, &s_plaq);
	//plaq is normalized by factor of 2.0/(VOL4D*NDIM*(NDIM-1)*NC), so one has to multiply by it again
	return (beta*2.0)/(VOL4D*NDIM*(NDIM-1)*NC)*(1.- plaq);
}

#ifdef _FERMIONS_
// sum_links phi*_i (M^+M)_ij^-1 phi_j
// here it is assumed that the rhs has already been computed... (because usually it has been during the leapfrog..)
hmc_complex s_fermion(hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi){
	return scalar_product(phi, MdaggerMphi);
}
#endif /* _FERMIONS_ */

//S_gauge + S_fermion + S_gaugemomenta
hmc_complex hamiltonian(hmc_gaugefield * field, hmc_float beta, hmc_gauge_momentum * p, hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi){
	hmc_complex result;
	hmc_complex tmp;
	(result) = {0.,0.};
	(result).re += s_gauge(field, beta);
#ifdef _FERMIONS_
	tmp  = s_fermion(phi, MdaggerMphi);
	complexaccumulate(&result, &tmp);
#endif
	//s_gm = 1/2*squarenorm(Pl)
	hmc_float s_gm;
	gaugemomenta_squarenorm(p, &s_gm);
	result.re += 0.5*s_gm;
	
	return result;
}

#ifdef _FERMIONS_
hmc_error generate_gaussian_spinorfield(hmc_spinor_field * out){
	// SL: this is a layer that calls the all-purpose hmc_complex gaussianly-distributed vector
	// with appropriate length and variance, i.e. SPINORFIELDSIZE and 1/2
	return gaussianComplexVector((hmc_complex *)out, SPINORFIELDSIZE, 0.5);
	// SL: not yet tested
}
#endif

hmc_error generate_gaussian_gauge_momenta(hmc_gauge_momentum * out){
	// SL: this is a layer that calls the all-purpose hmc_complex gaussianly-distributed vector
	// with appropriate length and variance, i.e. GAUGEMOMENTASIZE and 1
	// CP: hmc_gauge_momentum should be a real vector, so one should use GAUGEMOMENTASIZE/2 ?!?
	return gaussianComplexVector((hmc_complex *)out, GAUGEMOMENTASIZE, 1.0);
	// SL: not yet tested
}

hmc_error gauge_force(inputparameters * parameters, hmc_gaugefield * field, hmc_gauge_momentum * out){
	hmc_float beta = (*parameters).get_beta();
	int globalpos;
	//out loop over sites
	for(int t = 0; t < NTIME; t++){
		for(int n = 0; n < VOLSPACE; n++){
// 			globalpos = get_global_pos(n,t);
			for(int mu = 0; mu < NDIM; mu ++){
				for(int nu=0;nu<mu; nu++) {
					for(int i = 0; i< (NC*NC-1); i++){
						hmc_su3matrix prod;
// 						local_plaquette(field, &prod, n, t, mu, nu );
// 						hmc_float tmpfloat = 3. - trace_su3matrix(&prod).re;
// 						plaq += tmpfloat;
					}
				}
			}
		}
	}
	
	
	//TODO gauge force
	
	return HMC_SUCCESS;
}

#ifdef _FERMIONS_
//CP: it is assumed that phi_inv has been computed already
hmc_error fermion_force(inputparameters * parameters, hmc_gaugefield * field, hmc_spinor_field * phi, hmc_spinor_field * phi_inv, hmc_gauge_momentum * out){
	//TODO fermion force
	
	return HMC_SUCCESS;
}
#endif

//CP: this essentially calculates a hmc_gauge_momentum vector
//CP: it is assumed that phi_inv has been computed already
hmc_error force(inputparameters * parameters, hmc_gaugefield * field
#ifdef _FERMIONS_
	, hmc_spinor_field * phi, hmc_spinor_field * phi_inv
#endif
	, hmc_gauge_momentum * out){
	//CP: make sure that the output field is set to zero
	set_zero_gaugemomenta(out);
	//add contributions
#ifdef _FERMIONS_
	fermion_force(parameters, field, phi, phi_inv, out);
#endif
	gauge_force(parameters, field, out);
	return HMC_SUCCESS;
}

hmc_error metropolis(hmc_float rndnumber, hmc_float beta, hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi, hmc_gaugefield * field,	hmc_gauge_momentum * p, hmc_gaugefield * new_field, hmc_gauge_momentum * new_p){
	// takes:
	//		phi and beta as constant
	//		new/old versions of gaugefield and of momenta
	//		and a random number
	//		if it has to be, performs the change old->new, and returns true if there are no failures.
	hmc_complex h_old = hamiltonian(field, beta, p, phi, MdaggerMphi);
	hmc_complex h_new = hamiltonian(new_field, beta, new_p, phi, MdaggerMphi);
	if(h_old.im > projectioneps){
		printf("\n\tError: imaginary part in H_OLD [in function: metropolis(...)].\n");
		return HMC_COMPLEX_HAMILTONIANERROR;
	}
	if(h_new.im > projectioneps){
		printf("\n\tError: imaginary part in H_NEW [in function: metropolis(...)].\n");
		return HMC_COMPLEX_HAMILTONIANERROR;
	}
	//TODO export h_diff
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
	}
	return HMC_SUCCESS;
}


//CP: as in Gattringer/Lang, QCD on the Lattice, 8.2, p. 197
//TODO lateron, a multi-step alg. with different stepsizes for gauge and fermion force should be implemented
hmc_error leapfrog(inputparameters * parameters, hmc_gaugefield * u_in, hmc_gauge_momentum * p_in
#ifdef _FERMIONS_
	, hmc_spinor_field * phi
#endif
	, hmc_gaugefield * u_out, hmc_gauge_momentum * p_out
#ifdef _FERMIONS_
	, hmc_spinor_field * phi_inv
#endif
	){
	//TODO get steps and stepsize from parameters (and perhaps give them more fancy names (tau...))
	hmc_float stepsize;
	int steps;	
	int k;
	hmc_float stepsize_half = 0.5*stepsize;
	hmc_float kappa, mu, theta, chem_pot_re, chem_pot_im;
	int cgmax;
	//intermediate states u_k, p_-1/2 and p_1/2 (CP: i guess one only needs one hmc_gaugefield for the update here, the second in the alg. is unused..)
	hmc_gaugefield u_next;
	
	//CP: one needs the cg solver here since one wants to invert MdaggerM
	int use_cg = TRUE;
	
	hmc_gauge_momentum* p_prev = new hmc_gauge_momentum[GAUGEMOMENTASIZE];
	hmc_gauge_momentum* p_next = new hmc_gauge_momentum[GAUGEMOMENTASIZE];
	hmc_gauge_momentum* force_tmp = new hmc_gauge_momentum[GAUGEMOMENTASIZE];
	//TODO check usage of source here
#ifdef _FERMIONS_
	hmc_spinor_field* source = new hmc_spinor_field[SPINORFIELDSIZE];
	set_zero_spinorfield(source);
#endif

	//initial step
	copy_gaugefield(u_in, &u_next);
	copy_gaugemomenta(p_in, p_prev);

#ifdef _FERMIONS_
	//calc phi_inv;
	solver(phi, phi_inv, source, &u_next, kappa, mu, theta, chem_pot_re, chem_pot_im, cgmax, use_cg);
#endif
	force(parameters, &u_next 
		#ifdef _FERMIONS_
		, phi, phi_inv
		#endif
		, force_tmp );
	md_update_gauge_momenta(stepsize_half, p_prev, force_tmp, p_prev);
	
	//intermediate steps
	for(k = 1; k<steps-1; k++){
		//calc u_next
		md_update_gaugefield(stepsize, p_prev, &u_next);
#ifdef _FERMIONS_
		//calc phi_inv;
		solver(phi, phi_inv, source, &u_next, kappa, mu, theta, chem_pot_re, chem_pot_im, cgmax, use_cg);
#endif
		//calc p_next
		force(parameters, &u_next 
			#ifdef _FERMIONS_
			, phi, phi_inv
			#endif
			, force_tmp );
		md_update_gauge_momenta(stepsize, p_prev,force_tmp, p_next);
		copy_gaugemomenta(p_next, p_prev);
	}
	
	//final step
	md_update_gaugefield(stepsize, p_prev, &u_next);
#ifdef _FERMIONS_
	//calc phi_inv;
	solver(phi, phi_inv, source, &u_next, kappa, mu, theta, chem_pot_re, chem_pot_im, cgmax, use_cg);
#endif
	force(parameters, &u_next 
		#ifdef _FERMIONS_
		, phi, phi_inv
		#endif
		, force_tmp);
	md_update_gauge_momenta(stepsize_half, p_prev,force_tmp, p_next);

	//copy final results
	copy_gaugefield(&u_next, u_out);
	copy_gaugemomenta(p_next, p_out);
	
	delete [] p_prev;
	delete [] p_next;
	delete [] force_tmp;
#ifdef _FERMIONS_
	delete [] source;
#endif
	
	return HMC_SUCCESS;
}

hmc_error construct_3x3_combination(hmc_float beta_0, hmc_float gamma_0, hmc_float beta[], hmc_float gamma[], hmc_3x3matrix out){
	// called by build_su3matrix_by_exponentiation in case of "smart" approach
	// takes the 2*(8+1) real parameters beta_0, gamma_0, beta[8], gamma[8] and compiles
	// all components of the generic 3x3 complex matrix that is the linear combination of identity+generators
	hmc_float redb8 = beta[7] * F_1_2S3;
	hmc_float redg8 = gamma[7]* F_1_2S3;
	out[0][0].re = beta_0 + 0.5*beta[2] + redb8;
	out[0][0].im = gamma_0 + 0.5*gamma[2] + redg8;
	out[0][1].re = 0.5*(beta[0]+gamma[1]);
	out[0][1].im = 0.5*(gamma[0]-beta[1]);
	out[0][2].re = 0.5*(beta[3]+gamma[4]);
	out[0][2].im = 0.5*(gamma[3]-beta[4]);
	out[1][0].re = 0.5*(beta[0]-gamma[1]);
	out[1][0].im = 0.5*(gamma[0]+beta[1]);
	out[1][1].re = beta_0 - 0.5*beta[2] + redb8;
	out[1][1].im = gamma_0 - 0.5*gamma[2] + redg8;
	out[1][2].re = 0.5*(beta[5]+gamma[6]);
	out[1][2].im = 0.5*(gamma[5]-beta[6]);
	out[2][0].re = 0.5*(beta[3]-gamma[4]);
	out[2][0].im = 0.5*(gamma[3]+beta[4]);
	out[2][1].re = 0.5*(beta[5]-gamma[6]);
	out[2][1].im = 0.5*(gamma[5]+beta[6]);
	out[2][2].re = beta_0 + 2*redb8;
	out[2][2].im = gamma_0 - 2*redg8;
	return HMC_SUCCESS;
}

hmc_error build_su3matrix_by_exponentiation(hmc_algebraelement in, hmc_su3matrix* out, hmc_float epsilon){
	// SL: this takes 8 real numbers and builds the su3 matrix exp(i*epsilon*p_i*T_i)
	//     either by truncated "smart" series expansion to order eps^2 or eps^3, or by direct evaluation of the series
	// NOTE: the code here is strictly SU(3)-specific! (and mostly generated by other code...)

	#ifdef _USE_MORNGINGSTAR_PEARDON_
	//TODO CP: this has to be implemented
	
	#else

	#ifndef _EXPONENTIATE_ALGEBRA_ALL_ORDERS_

	// Phase 1: calculate (number) coefficients for the reconstruction
	hmc_float beta_0, gamma_0;
	hmc_float beta[8], gamma[8];
	// the above are: coefficients for identity (beta_0+i*gamma_0), and coefficients for the T_l, namely T_L*(beta_l + i*gamma_l).
	#ifdef _EXPONENTIATE_ALGEBRA_ORDER_2_
		hmc_float pR, pQ[8], pPtilde[8];

		pR = in[0]*in[0]+in[1]*in[1]+in[2]*in[2]+in[3]*in[3]+in[4]*in[4]+in[5]*in[5]+in[6]*in[6]+in[7]*in[7];

		pPtilde[0] = (F_1_S3*in[7]*in[0])+(F_1_2*(in[5]*in[3]+in[6]*in[4]));
		pPtilde[1] = (F_1_S3*in[7]*in[1])+(F_1_2*(in[5]*in[4]-in[6]*in[3]));
		pPtilde[2] = (F_1_S3*in[7]*in[2]);
		pPtilde[3] = (F_1_2*(in[5]*in[0]-in[6]*in[1]+in[3]*in[2]))-(F_1_2S3*in[7]*in[3]);
		pPtilde[4] = (F_1_2*(in[6]*in[0]+in[5]*in[1]+in[4]*in[2]))-(F_1_2S3*in[7]*in[4]);
		pPtilde[5] = (F_1_2*(in[3]*in[0]+in[4]*in[1]-in[5]*in[2]))-(F_1_2S3*in[7]*in[5]);
		pPtilde[6] = (F_1_2*(in[4]*in[0]-in[3]*in[1]-in[6]*in[2]))-(F_1_2S3*in[7]*in[6]);
		pPtilde[7] = 0.0;

		pQ[2] = (F_1_2*(in[3]*in[3]+in[4]*in[4]-in[5]*in[5]-in[6]*in[6]));
		pQ[7] = (F_1_S3*(in[0]*in[0]+in[1]*in[1]+in[2]*in[2]-in[7]*in[7]))-(F_1_2S3*(in[3]*in[3]+in[4]*in[4]+in[5]*in[5]+in[6]*in[6]));
		pQ[0] = pQ[1] = pQ[3] = pQ[4] = pQ[5] = pQ[6] = 0.0;

		hmc_float eps_squared = epsilon*epsilon;
		beta_0  = 1.0 - (eps_squared*pR/12.);
		gamma_0 = 0.0;
		for(int il=0;il<8;il++){
			beta[il] = -0.25*eps_squared*(pQ[il]+2.0*pPtilde[il]);
			gamma[il]= epsilon*in[il];
		}
	#endif // EXPONENTIATE_ALGEBRA_ORDER_2
	#ifdef _EXPONENTIATE_ALGEBRA_ORDER_3_
		hmc_float pR, pQ[8], pPtilde[8];
		hmc_float pT, pS[8], pU[8];

		pR = in[0]*in[0]+in[1]*in[1]+in[2]*in[2]+in[3]*in[3]+in[4]*in[4]+in[5]*in[5]+in[6]*in[6]+in[7]*in[7];
		
		pPtilde[0] = (F_1_S3*in[7]*in[0])+(F_1_2*(in[5]*in[3]+in[6]*in[4]));
		pPtilde[1] = (F_1_S3*in[7]*in[1])+(F_1_2*(in[5]*in[4]-in[6]*in[3]));
		pPtilde[2] = (F_1_S3*in[7]*in[2]);
		pPtilde[3] = (F_1_2*(in[5]*in[0]-in[6]*in[1]+in[3]*in[2]))-(F_1_2S3*in[7]*in[3]);
		pPtilde[4] = (F_1_2*(in[6]*in[0]+in[5]*in[1]+in[4]*in[2]))-(F_1_2S3*in[7]*in[4]);
		pPtilde[5] = (F_1_2*(in[3]*in[0]+in[4]*in[1]-in[5]*in[2]))-(F_1_2S3*in[7]*in[5]);
		pPtilde[6] = (F_1_2*(in[4]*in[0]-in[3]*in[1]-in[6]*in[2]))-(F_1_2S3*in[7]*in[6]);
		pPtilde[7] = 0.0;

		pQ[2] = (F_1_2*(in[3]*in[3]+in[4]*in[4]-in[5]*in[5]-in[6]*in[6]));
		pQ[7] = (F_1_S3*(in[0]*in[0]+in[1]*in[1]+in[2]*in[2]-in[7]*in[7]))-(F_1_2S3*(in[3]*in[3]+in[4]*in[4]+in[5]*in[5]+in[6]*in[6]));
		pQ[0] = pQ[1] = pQ[3] = pQ[4] = pQ[5] = pQ[6] = 0.0;
		
		pT = (in[0]*(pQ[0]+2*pPtilde[0]))+(in[1]*(pQ[1]+2*pPtilde[1]))+(in[2]*(pQ[2]+2*pPtilde[2]))+(in[3]*(pQ[3]+2*pPtilde[3]))
			+(in[4]*(pQ[4]+2*pPtilde[4]))+(in[5]*(pQ[5]+2*pPtilde[5]))+(in[6]*(pQ[6]+2*pPtilde[6]))+(in[7]*(pQ[7]+2*pPtilde[7]));
			
		pS[0] = (F_1_2*(in[5]*(pQ[3]+2*pPtilde[3])+in[6]*(pQ[4]+2*pPtilde[4])))+(F_1_S3*in[7]*(pQ[0]+2*pPtilde[0]));
		pS[1] = (F_1_2*(in[5]*(pQ[4]+2*pPtilde[4])-in[6]*(pQ[3]+2*pPtilde[3])))+(F_1_S3*in[7]*(pQ[1]+2*pPtilde[1]));
		pS[2] = (F_1_S3*in[7]*(pQ[2]+2*pPtilde[2]));
		pS[3] = (F_1_2*(in[3]*(pQ[2]+2*pPtilde[2])+in[5]*(pQ[0]+2*pPtilde[0])-in[6]*(pQ[1]+2*pPtilde[1])))-(F_1_2S3*in[7]*(pQ[3]+2*pPtilde[3]));
		pS[4] = (F_1_2*(in[4]*(pQ[2]+2*pPtilde[2])+in[5]*(pQ[1]+2*pPtilde[1])+in[6]*(pQ[0]+2*pPtilde[0])))-(F_1_2S3*in[7]*(pQ[4]+2*pPtilde[4]));
		pS[5] = (F_1_2*(in[3]*(pQ[0]+2*pPtilde[0])+in[4]*(pQ[1]+2*pPtilde[1])-in[5]*(pQ[2]+2*pPtilde[2])))-(F_1_2S3*in[7]*(pQ[5]+2*pPtilde[5]));
		pS[6] = (F_1_2*(in[4]*(pQ[0]+2*pPtilde[0])-in[3]*(pQ[1]+2*pPtilde[1])-in[6]*(pQ[2]+2*pPtilde[2])))-(F_1_2S3*in[7]*(pQ[6]+2*pPtilde[6]));
		pS[7] = 0.0;

		pU[2] = (F_1_2*(in[3]*(pQ[3]+2*pPtilde[3])+in[4]*(pQ[4]+2*pPtilde[4])-in[5]*(pQ[5]+2*pPtilde[5])-in[6]*(pQ[6]+2*pPtilde[6])));
		pU[7] = (F_1_S3*(in[0]*(pQ[0]+2*pPtilde[0])+in[1]*(pQ[1]+2*pPtilde[1])+in[2]*(pQ[2]+2*pPtilde[2])-in[7]*(pQ[7]+2*pPtilde[7])))
			-(F_1_2S3*(in[3]*(pQ[3]+2*pPtilde[3])+in[4]*(pQ[4]+2*pPtilde[4])+in[5]*(pQ[5]+2*pPtilde[5])+in[6]*(pQ[6]+2*pPtilde[6])));
		pU[0] = pU[1] = pU[3] = pU[4] = pU[5] = pU[6] = 0.0;
		
		hmc_float eps_squared, eps_cubed;
		eps_squared = epsilon*epsilon;
		eps_cubed = eps_squared*epsilon;
		beta_0  = 1.0 - (eps_squared*pR/12.);
		gamma_0 = -eps_cubed*pT/72.;
		for(int il=0;il<8;il++){
			beta[il]  = -0.25*eps_squared*(pQ[il]+2.0*pPtilde[il]);
			gamma[il] = epsilon*in[il] - (eps_cubed/72.)*(2.*pR*in[il] + 3.*pU[il] + 6.*pS[il]);
		}
		#endif // EXPONENTIATE_ALGEBRA_ORDER_3
	
		// Phase 2: build a (generic 3x3) matrix from beta's and gamma's
		hmc_3x3matrix combination;
		construct_3x3_combination(beta_0, gamma_0, beta, gamma, combination);
		
		// Phase 3: project back onto SU(3). This requires knowledge whether one is using _RECONSTRUCT_TWELVE_ when filling "out"
		#ifdef _RECONSTRUCT_TWELVE_
			*out[0] = combination[0][0];
			*out[1] = combination[0][1];
			*out[2] = combination[0][2];
			*out[3] = combination[1][0];
			*out[4] = combination[1][1];
			*out[5] = combination[1][2];
		#else
			*out[0][0] = combination[0][0];
			*out[0][1] = combination[0][1];
			*out[0][2] = combination[0][2];
			*out[1][0] = combination[1][0];
			*out[1][1] = combination[1][1];
			*out[1][2] = combination[1][2];
			*out[2][0] = combination[2][0];
			*out[2][1] = combination[2][1];
			*out[2][2] = combination[2][2];
		#endif // _RECONSTRUCT_TWELVE_
		project_su3(out);
	
	#else  // (ifndef) _EXPONENTIATE_ALGEBRA_ALL_ORDERS_s
		// this is TODO ! the case where one actually evaluates -as matrices- many orders of exp(i*e*P)=1+i*e*P+(1/2)(i*e*P) + ...
	#endif // EXPONENTIATE_ALGEBRA_ALL_ORDERS
	
	#endif // _USE_MORNGINGSTAR_PEARDON_
	
	return HMC_SUCCESS;
}
