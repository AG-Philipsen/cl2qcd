/**
 * @file operations used for the molecular dynamics update
 */

__kernel void md_update_gaugefield(){
	

//molecular dynamics update for the gaugefield:
//u_out = exp(i eps p_in) u_in
// hmc_error md_update_gaugefield(hmc_float eps, hmc_algebraelement2 * p_in, hmc_gaugefield * u_inout){
// 	int index;
// 	for(int t = 0; t<NTIME; t++){
// 		for(int pos = 0; pos < VOLSPACE; pos++){
// 			for(int mu = 0; mu<NDIM; mu++){
// 			hmc_su3matrix tmp;
// 			hmc_su3matrix tmp2;
// 			index= get_global_link_pos(mu, pos, t);
// 			// an su3 algebra element has NC*NC-1 = 8 hmc_float entries
// 			// &(p_in[index*8]) should point to the right position for the pos-th element of the long gaugemomentum vector p_in
// 			build_su3matrix_by_exponentiation((p_in[index]), &tmp2, eps);
// 			get_su3matrix(&tmp, u_inout, pos, t, mu);
// 			accumulate_su3matrix_prod( &tmp2, &tmp);
// 			put_su3matrix(u_inout, &tmp2, pos, t, mu);
// 	}}}
// 	return HMC_SUCCESS;
// }
	
	
}
	
__kernel void md_update_gaugemomenta(){
	
	//CP: molecular dynamics update for the gauge momenta:
//p_out = p_in - eps/2 force(u_in, phi)
//it is assumed that the force term has already been computed. then one only has real-vectors and this is essentially adding one vector to another...
// hmc_error md_update_gauge_momenta(hmc_float eps, hmc_algebraelement2 * p_inout, hmc_algebraelement2 * force_in){
// 	for(int i = 0; i<GAUGEMOMENTASIZE2; i++){
// 		acc_factor_times_algebraelement(&p_inout[i], -1.*eps, force_in[i]);
// 	}
// 	return HMC_SUCCESS;
// }
	
}

__kernel void s_gauge(){
	// beta * sum_links sum_nu>mu ( 3 - Tr Re Plaquette )
//CP: since one is only interested in differences of s_gauge, the constant part can be left out!!
// hmc_float s_gauge(hmc_gaugefield * field, hmc_float beta){
// 	/** @TODO CP: implement saving of plaquette measurement (and possibly t_plaq and s_plaq and also polyakov-loop??) */
// 	hmc_float plaq=0;
// 	//CP: alternative method: use already existing plaquette-functions
// 	hmc_float t_plaq;
// 	hmc_float s_plaq;
// 	plaq = plaquette(field, &t_plaq, &s_plaq);
// 	//plaq is normalized by factor of 2.0/(VOL4D*NDIM*(NDIM-1)*NC), so one has to divide by it again
// 	hmc_float factor = 2.0/static_cast<hmc_float>(VOL4D*NDIM*(NDIM-1)*NC);
// 	return beta/6./factor*( - plaq);
	
	
}

__kernel void s_fermion(){
	// sum_links phi*_i (M^+M)_ij^-1 phi_j
// here it is assumed that the rhs has already been computed... (because usually it has been during the leapfrog..)
// hmc_float s_fermion(inputparameters*parameters, hmc_gaugefield * field, hmc_spinor_field * phi, hmc_spinor_field * QplusQminusphi_inv){
// 	//CP: phi is not needed after this, so it can be used to store M (QplusQminus_inv)
// 	Qminus(parameters, QplusQminusphi_inv, field, phi);
// 	return global_squarenorm(phi);
	
}
