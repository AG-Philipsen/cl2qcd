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

//phi = D chi
hmc_error md_update_spinorfield(hmc_spinor_field * in, hmc_spinor_field * out, hmc_gaugefield * field, inputparameters * parameters){
	//TODO extract needed parameters from paramters, or transform this into M itself
	hmc_float kappa, mu, theta, chem_pot_re, chem_pot_im;
	//TODO check again if it is M or Mdagger here
	Mdagger(in, out, field, kappa, mu, theta, chem_pot_re, chem_pot_im);
	
	return HMC_SUCCESS;
}

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

// sum_links phi*_i (M^+M)_ij^-1 phi_j
// here it is assumed that the rhs has already been computed... (because usually it has been during the leapfrog..)
hmc_complex s_fermion(hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi){
	return scalar_product(phi, MdaggerMphi);
}

//S_gauge + S_fermion + S_gaugemomenta
hmc_complex hamiltonian(hmc_gaugefield * field, hmc_float beta, hmc_gauge_momentum * p, hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi){
	hmc_complex result;
	hmc_complex tmp;
	(result) = {0.,0.};
	(result).re += s_gauge(field, beta);
	tmp  = s_fermion(phi, MdaggerMphi);
	complexaccumulate(&result, &tmp);
	//s_gm = 1/2*squarenorm(Pl)
	hmc_float s_gm;
	gaugemomenta_squarenorm(p, &s_gm);
	result.re += 0.5*s_gm;
	
	return result;
}

hmc_error generate_gaussian_spinorfield(hmc_spinor_field * out){
	//TODO generate gaussian initial spinorfield
	
	return HMC_SUCCESS;
}

hmc_error generate_gaussian_gauge_momenta(hmc_gauge_momentum * out){
	//TODO generate gaussian initial gauge momenta
	
	return HMC_SUCCESS;
}

hmc_error gauge_force(inputparameters * parameters, hmc_gaugefield * field, hmc_gauge_momentum * out){
	//TODO gauge force
	
	return HMC_SUCCESS;
}

//CP: it is assumed that phi_inv has been computed already
hmc_error fermion_force(inputparameters * parameters, hmc_gaugefield * field, hmc_spinor_field * phi, hmc_spinor_field * phi_inv, hmc_gauge_momentum * out){
	//TODO fermion force
	
	return HMC_SUCCESS;
}

//CP: this essentially calculates a hmc_gauge_momentum vector
//CP: it is assumed that phi_inv has been computed already
hmc_error force(inputparameters * parameters, hmc_gaugefield * field, hmc_spinor_field * phi, hmc_gauge_momentum * out, hmc_spinor_field * phi_inv){
	//CP: make sure that the output field is set to zero
	set_zero_gaugemomenta(out);
	//add contributions
	gauge_force(parameters, field, out);
	fermion_force(parameters, field, phi, phi_inv, out);
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
hmc_error leapfrog(inputparameters * parameters, hmc_gaugefield * u_in, hmc_gauge_momentum * p_in, hmc_spinor_field * phi, hmc_gaugefield * u_out, hmc_gauge_momentum * p_out, hmc_spinor_field * phi_inv){
	//TODO get steps and stepsize from parameters (and perhaps give them more fancy names (tau...))
	hmc_float stepsize;
	int steps;	
	int k;
	hmc_float stepsize_half = 0.5*stepsize;
	//intermediate states u_k, p_-1/2 and p_1/2 (CP: i guess one only needs one hmc_gaugefield for the update here, the second in the alg. is unused..)
	hmc_gaugefield u_next;
	
	hmc_gauge_momentum* p_prev = new hmc_gauge_momentum[GAUGEMOMENTASIZE];
	hmc_gauge_momentum* p_next = new hmc_gauge_momentum[GAUGEMOMENTASIZE];
	hmc_gauge_momentum* force_tmp = new hmc_gauge_momentum[GAUGEMOMENTASIZE];
	
	//initial step
	copy_gaugefield(u_in, &u_next);
	copy_gaugemomenta(p_in, p_prev);

	//TODO calc phi_inv;
	force(parameters, &u_next, phi, force_tmp, phi_inv);
	md_update_gauge_momenta(stepsize_half, p_prev, force_tmp, p_prev);
	
	//intermediate steps
	for(k = 1; k<steps-1; k++){
		//calc u_next
		md_update_gaugefield(stepsize, p_prev, &u_next);
		//TODO calc phi_inv;
		//calc p_next
		force(parameters, &u_next, phi, force_tmp, phi_inv);
		md_update_gauge_momenta(stepsize, p_prev,force_tmp, p_next);
		copy_gaugemomenta(p_next, p_prev);
	}
	
	//final step
	md_update_gaugefield(stepsize, p_prev, &u_next);
	//TODO calc phi_inv;
	force(parameters, &u_next, phi, force_tmp, phi_inv);
	md_update_gauge_momenta(stepsize_half, p_prev,force_tmp, p_next);

	//copy final results
	copy_gaugefield(&u_next, u_out);
	copy_gaugemomenta(p_next, p_out);
	
	delete [] p_prev;
	delete [] p_next;
	delete [] force_tmp;
	
	return HMC_SUCCESS;
}





