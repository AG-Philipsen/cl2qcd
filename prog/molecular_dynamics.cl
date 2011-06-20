/**
 * @file operations used for the molecular dynamics update
 */

//molecular dynamics update for the gaugefield:
//u_out = exp(i eps p_in) u_in
__kernel void md_update_gaugefield(hmc_float eps, __global ae * p_in, __global ocl_s_gaugefield * u_inout){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int n,t, index;
	Matrixsu3 tmp;
	Matrixsu3 tmp2;

	for(int id_tmp = id; id_tmp < GAUGEMOMENTASIZE; id_tmp += global_size) {	
		/** @todo this must be done more efficient */
		if(id_tmp%2 == 0) get_even_site(id_tmp/2, &n, &t);
		else get_odd_site(id_tmp/2, &n, &t);
		for(int mu = 0; mu<NDIM; mu++){
			index= get_global_link_pos(mu, n, t);
			// an su3 algebra element has NC*NC-1 = 8 hmc_float entries
 			// &(p_in[index*8]) should point to the right position for the pos-th element of the long gaugemomentum vector p_in
			//build_su3matrix_by_exponentiation((p_in[index]), &tmp2, eps);
			tmp = get_matrixsu3(u_inout, n, t, mu);
			tmp = multiply_matrixsu3( tmp2, tmp);
			put_matrixsu3(u_inout, tmp2, n, t, mu);
		}
	}	
}

//CP: molecular dynamics update for the gauge momenta:
//p_out = p_in - eps/2 force(u_in, phi)
//it is assumed that the force term has already been computed. then one only has real-vectors and this is essentially adding one vector to another...	
__kernel void md_update_gaugemomenta(hmc_float eps, __global ae * p_inout, __global ae* force_in){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	
	for(int id_tmp = id; id_tmp < GAUGEMOMENTASIZE; id_tmp += global_size) {	
		p_inout[id_tmp] = acc_factor_times_algebraelement(p_inout[id_tmp], -1.*eps, force_in[id_tmp]);
	}
	
}

//these are deprecated since they are not really needed and can be executed from outside...
/*
// beta * sum_links sum_nu>mu ( 3 - Tr Re Plaquette )
//CP: since one is only interested in differences of s_gauge, the constant part can be left out!!
__kernel void s_gauge(){
// hmc_float s_gauge(hmc_gaugefield * field, hmc_float beta){
// 	/** @TODO CP: implement saving of plaquette measurement (and possibly t_plaq and s_plaq and also polyakov-loop??) 
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

*/