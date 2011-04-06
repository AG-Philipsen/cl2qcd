#include "host_hmc.h"

//TODO CP: all return values have to be revisited, since they should be real numbers and not complex one in the end

//TODO in the end hmc_gaugemomenta has to be a vector of NC*NC matrices!! This means that also the declaration of such variables has to be revisited



//TODO the following three functions should later go into a new, seperate file

//molecular dynamics update for the gauge momenta:
//p_out = p_in - eps/2 force(u_in, phi)
//it is assumed that the force term has already been computed
hmc_error md_update_gauge_momenta(hmc_float eps, hmc_gauge_momentum * p_in, hmc_gauge_momentum * force_in, hmc_gauge_momentum * p_out){
	//TODO CP: this is a memory-intense, lazy implementation
	hmc_complex tmp;
	for(int i = 0; i<GAUGEMOMENTASIZE; i++){
		tmp = force_in[i];
		complexmult_real(&tmp, &eps);
		p_out[i] = complexsubtract(&(p_in[i]), &tmp);
	}
		
	return HMC_SUCCESS;
}

//molecular dynamics update for the gaugefield:
//u_out = exp(i eps p_in) u_in
hmc_error md_update_gaugefield(hmc_float eps, hmc_gauge_momentum * p_in, hmc_gaugefield * u_in){
	//TODO CP: this is a memory-intense, lazy implementation
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
	//TODO check again if it is M or Mdagger
	Mdagger(in, out, field, kappa, mu, theta, chem_pot_re, chem_pot_im);
	
	return HMC_SUCCESS;
}


// beta * sum_links sum_nu>mu ( 3 - Tr Re Plaquette )
hmc_float s_gauge(hmc_gaugefield * field, hmc_float beta){
	hmc_float plaq=0;

	for(int t=0;t<NTIME;t++) {
		for(int n=0;n<VOLSPACE;n++) {
			for(int mu=0; mu<NDIM; mu++) {
				for(int nu=0;nu<mu; nu++) {
					hmc_su3matrix prod;
					local_plaquette(field, &prod, n, t, mu, nu );
					hmc_float tmpfloat = 3. - trace_su3matrix(&prod).re;
					plaq += tmpfloat;
				}}}}
	//normalize
	return beta*plaq*2.0/static_cast<hmc_float>(VOL4D*NDIM*(NDIM-1)*NC);
}

// sum_links phi*_i (M^+M)_ij^-1 phi_j
// here it is assumed that the rhs has already been computed... (because usually it has been during the leapfrog..)
hmc_complex s_fermion(hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi){
	return scalar_product(phi, MdaggerMphi);
}

//S_gauge + S_fermion + 1/2*squarenorm(Pl)
hmc_complex hamiltonian(hmc_gaugefield * field, hmc_float beta, hmc_gauge_momentum * p, hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi){
	hmc_complex result;
	hmc_complex tmp;
	(result) = {0.,0.};
	(result).re += s_gauge(field, beta);
	tmp  = s_fermion(phi, MdaggerMphi);
	complexaccumulate(&result, &tmp);
	
	//TODO calc squarenorm of p, factor 1/2??
	
	return result;
}

//TODO generate gaussian initial spinorfield
hmc_error generate_gaussian_spinorfield(hmc_spinor_field * out){
	
	
	return HMC_SUCCESS;
}


//TODO generate gaussian initial gauge momenta
hmc_error generate_gaussian_gauge_momenta(hmc_gauge_momentum * out){
	
	
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
hmc_error leapfrog(inputparameters * parameters, hmc_gaugefield * u_in, hmc_gauge_momentum * p_in, hmc_spinor_field * phi, hmc_gaugefield * u_out, hmc_gauge_momentum * p_out, hmc_spinor_field * phi_inv){
	//TODO get steps and stepsize from parameters (and perhaps give them more fancy names (tau...))
	hmc_float stepsize;
	int steps;	
	int k;
	hmc_float stepsize_half = 0.5*stepsize;
	//intermediate states u_k, u_k-1, p_-1/2 and p_1/2
	hmc_gaugefield u_next;
	
	hmc_spinor_field* p_prev = new hmc_gauge_momentum[GAUGEMOMENTASIZE];
	hmc_spinor_field* p_next = new hmc_gauge_momentum[GAUGEMOMENTASIZE];
	hmc_spinor_field* force_tmp = new hmc_gauge_momentum[GAUGEMOMENTASIZE];
	
	//initial step
	copy_gaugefield(u_in, &u_next);
	copy_gaugemomenta(p_in, p_prev);
	
	//TODO calculating the force term has to be inserted at three points!!
	
	//CP: i guess one only needs one hmc_gaugefield for the update here...
	md_update_gauge_momenta(stepsize_half, p_prev, force_tmp, p_prev);
	
	//intermediate steps
	for(k = 1; k<steps-1; k++){
		//calc u_next
		md_update_gaugefield(stepsize, p_prev, &u_next);
		//calc p_next
		md_update_gauge_momenta(stepsize, p_prev,force_tmp, p_next);
		copy_gaugemomenta(p_next, p_prev);
	}
	
	//final step
	md_update_gaugefield(stepsize, p_prev, &u_next);
	md_update_gauge_momenta(stepsize_half, p_prev,force_tmp, p_next);

	//copy final results
	copy_gaugefield(&u_next, u_out);
	copy_gaugemomenta(p_next, p_out);
	
	delete [] p_prev;
	delete [] p_next;
	delete [] force_tmp;
	
	return HMC_SUCCESS;
}

//TODO force
//CP: this essentially calculates a hmc_gauge_momentum vector
//CP: one has to save the MdaggerM^-1 phi from here (in phi_inv) so that one does not have to do its calculation again later... (in the end, this should be optional!!)
hmc_error force(inputparameters * parameters, hmc_gaugefield * field, hmc_spinor_field * phi, hmc_gauge_momentum * out, hmc_spinor_field * phi_inv){
	//TODO get all variables needed from parameters
	
	
	
	return HMC_SUCCESS;
}



//TODO gauge force



//TODO fermion force


hmc_error perform_hybrid_monte_carlo(inputparameters * parameters){
	//init section
	//CP: here, one has to get all the parameters necessary from the inputparameters
	// suppose that in parameters the number of hmc_iterations is given, store them in a variable here...
	int hmc_iter; //= ...
	int iter;
	int err = 0;
	//beta has to be saved here to give it to the metropolis step, all other parameters can be given via parameters
	hmc_float beta;
	// TODO give seed meaningful value, perhaps read it in from parameters
	int seed = 89343894543895;
	Random rnd_gen (seed);
	hmc_float rnd;
	
	hmc_spinor_field* phi = new hmc_spinor_field[SPINORFIELDSIZE];
	hmc_spinor_field* phi_inv = new hmc_spinor_field[SPINORFIELDSIZE];
	hmc_spinor_field* chi = new hmc_spinor_field[SPINORFIELDSIZE];
	hmc_gauge_momentum* p = new hmc_gauge_momentum[GAUGEMOMENTASIZE];
	hmc_gauge_momentum* new_p = new hmc_gauge_momentum[GAUGEMOMENTASIZE];
	hmc_gaugefield field;
	hmc_gaugefield new_field;
	
	
	//CP: depending on the input parameters, the gaugefield should be initialized accordingly...
	set_gaugefield_cold(&field);
	
	
	//main hmc-loop
	for(iter = 0; iter < hmc_iter; iter ++){
		//init gauge_momenta
		//TODO perhaps write a wrapper that automatically evaluates the err's
		err = generate_gaussian_gauge_momenta(p);
		
		//init/update spinorfield phi
		err = generate_gaussian_spinorfield(chi);
		err = md_update_spinorfield(chi, phi, &field, parameters);
		
		//update gaugefield and gauge_momenta via leapfrog
		err = leapfrog(parameters, &field, p, phi, &new_field, new_p, phi_inv);
		
		//metropolis step: afterwards, the updated config is again in field and p
		//generate new random-number
		rnd = rnd_gen.doub();
		err = metropolis(rnd, beta, phi, phi_inv, &field, p, &new_field, new_p);
		
		//TODO if wanted, measurements can be added here...
	}
	
	delete [] phi;
	delete [] phi_inv;
	delete [] chi;
	delete [] p;
	delete [] new_p;
	
	return HMC_SUCCESS;
}

