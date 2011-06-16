/**
 * @file operations used for the molecular dynamics update
 */


#ifdef _FERMIONS_
// sum_links phi*_i (M^+M)_ij^-1 phi_j
// here it is assumed that the rhs has already been computed... (because usually it has been during the leapfrog..)
hmc_float s_fermion(inputparameters*parameters, hmc_gaugefield * field, hmc_spinor_field * phi, hmc_spinor_field * QplusQminusphi_inv){
	//CP: phi is not needed after this, so it can be used to store M (QplusQminus_inv)
	Qminus(parameters, QplusQminusphi_inv, field, phi);
	return global_squarenorm(phi);
}
#endif /* _FERMIONS_ */

//Without fermions here!!! H = S_gauge + S_gaugemomenta
hmc_float hamiltonian(hmc_gaugefield * field, hmc_float beta, hmc_algebraelement2 * p){
	hmc_float result;
	(result) = 0.;
	(result) += s_gauge(field, beta);
	//s_gm = 1/2*squarenorm(Pl)
	hmc_float s_gm;
	gaugemomenta_squarenorm(p, &s_gm);
	result += 0.5*s_gm;
	
	return result;
}


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
hmc_error md_update_gaugefield(hmc_float eps, hmc_algebraelement2 * p_in, hmc_gaugefield * u_inout){
	int index;
	for(int t = 0; t<NTIME; t++){
		for(int pos = 0; pos < VOLSPACE; pos++){
			for(int mu = 0; mu<NDIM; mu++){
			hmc_su3matrix tmp;
			hmc_su3matrix tmp2;
			index= get_global_link_pos(mu, pos, t);
			// an su3 algebra element has NC*NC-1 = 8 hmc_float entries
			// &(p_in[index*8]) should point to the right position for the pos-th element of the long gaugemomentum vector p_in
			build_su3matrix_by_exponentiation((p_in[index]), &tmp2, eps);
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







hmc_error metropolis(hmc_float rndnumber, inputparameters * parameter,  
										 #ifdef _FERMIONS_
										 hmc_spinor_field * phi, hmc_spinor_field * phi_inv, hmc_float energy_init, 
										 #endif
										 hmc_gaugefield * field,	hmc_algebraelement2 * p, hmc_gaugefield * new_field, hmc_algebraelement2 * new_p){
	// takes:
	//		phi and beta as constant
	//		new/old versions of gaugefield and of momenta
	//		and a random number
	//		if it has to be, performs the change old->new, and returns true if there are no failures.
	hmc_float h_old = hamiltonian(field, (*parameter).get_beta(), p);
	hmc_float h_new = hamiltonian(new_field, (*parameter).get_beta(), new_p);
	
	#ifdef _FERMIONS_
	hmc_float tmp;
	tmp  = s_fermion(parameter, field, phi, phi_inv);
	#endif
	h_new += tmp;
	h_old += energy_init;

	/** @todo CP:  export h_diff */
	hmc_float h_diff = h_old - h_new;
	hmc_float compare_prob;
	if(h_diff<0){
		compare_prob = exp(h_diff);
	}else{
		compare_prob = 1.0;
	}
	//CP: Debugging
	cout << "\t\th_old:\t"<<h_old<<"\tenergy_init:\t"<<energy_init<<endl;
	cout << "\t\th_new:\t"<<h_new<<"\tfermion con:\t"<<tmp<<endl;
	cout <<"\t\tacc_prob:\t"<<compare_prob<<endl;
	
	
	
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
