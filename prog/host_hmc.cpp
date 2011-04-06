#include "host_hmc.h"

//TODO CP: all return values have to be revisited, since they should be real numbers and not complex one in the end

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
hmc_complex s_fermion(hmc_gaugefield * field, hmc_spinor_field * phi){
	hmc_complex result = {0., 0.};
	for(int t=0;t<NTIME;t++) {
		for(int n=0;n<VOLSPACE;n++) {
			//TODO do stuff
		}}
	return result;
}

//S_gauge + S_fermion + 1/2*squarenorm(Pl)
hmc_complex hamiltonian(hmc_gaugefield * field, hmc_float beta, hmc_gauge_momentum * p, hmc_spinor_field * phi){
	hmc_complex result;
	hmc_complex tmp;
	(result) = {0.,0.};
	(result).re += s_gauge(field, beta);
	tmp  = s_fermion(field, phi);
	complexaccumulate(&result, &tmp);
	
	//TODO calc squarenorm of p, factor 1/2??
	
	return result;
}

//TODO spinor update



//TODO generate gauge momenta



//TODO metropolis step
hmc_error metropolis(hmc_float rndnumber, hmc_float beta, hmc_spinor_field * phi, hmc_gaugefield * field,
	hmc_gauge_momentum * p, hmc_gaugefield * new_field, hmc_gauge_momentum * new_p){
	// takes:
	//		phi and beta as constant
	//		new/old versions of gaugefield and of momenta
	//		and a random number
	//		if it has to be, performs the change old->new, and returns true if there are no failures.
	hmc_complex h_old = hamiltonian(field, beta, p, phi);
	hmc_complex h_new = hamiltonian(new_field, beta, new_p, phi);
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


//TODO leapfrog integration
//CP: as in Gattringer/Lang, QCD on the Lattice, 8.2, p. 197
void leapfrog(const hmc_float stepsize, const int steps, hmc_gaugefield * u_in, hmc_gauge_momentum * p_in, hmc_spinor_field * phi, hmc_gaugefield * u_out, hmc_gauge_momentum * p_out){
	int k;
	hmc_float stepsize_half = 0.5*stepsize;
	//intermediate states u_k, u_k-1, p_-1/2 and p_1/2
	hmc_gaugefield u_next;
	hmc_gaugefield u_prev;
	hmc_gauge_momentum p_prev;
	hmc_gauge_momentum p_next;
	//initial step
	//TODO copy_gaugefield(u_prev = u_in);
	//TODO copy_gauge_momenta(p_prev = p_in);
	
	//TODO mc_update_gauge_momenta(stepsize_half, p_prev, u_prev, phi, p_prev);
	
	
	//intermediate steps
	for(k = 1; k<steps-1; k++){
		//calc u_next
		//TODO mc_update_gaugefield(stepsize, p_prev, u_prev, u_next);
		//calc p_next
		//TODO mc_update_gauge_momenta(stepsize, p_prev, u_next, phi, p_next);
		//TODO copy_gaugefield(u_prev = u_next);
		//TODO copy_gauge_momenta(p_prev = p_next);
	}
	
	//final step
	//TODO mc_update_gaugefield(stepsize, p_prev, u_prev, field_out);
	//TODO mc_update_gauge_momenta(stepsize_half, p_prev, field_out, phi, p_next);

	//copy results
	//TODO copy_gaugefield(u_next = u_out);
	//TODO copy_gauge_momenta(p_next = p_out);
	
	return;
}


//TODO gauge force



//TODO fermion force



//TODO force
hmc_error force(hmc_gaugefield field, hmc_spinor_field * phi, hmc_gauge_momentum * out){
	
	
	
	
	return HMC_SUCCESS;
}