#include "host_hmc.h"

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



//TODO gauge force



//TODO fermion force