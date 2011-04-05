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



//TODO leapfrog integration



//TODO gauge force



//TODO fermion force