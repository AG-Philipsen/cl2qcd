/**
 * @file integrators needed by HMC
 */


//it is assumed that gaugefield and gaugemomentum have been set to the old ones already
hmc_error leapfrog(inputparameters * parameters, 
									 #ifdef _FERMIONS_
									 hmc_spinor_field * phi, hmc_spinor_field * phi_inv, 
									 #endif
									 hmc_gaugefield * u_out, hmc_algebraelement2 * p_out	){
	// CP: it operates directly on the fields p_out and u_out
	int steps = (*parameters).get_integrationsteps1() ;	
	hmc_float stepsize = ((*parameters).get_tau()) /((hmc_float) steps);
	int k;
	hmc_float stepsize_half = 0.5*stepsize;
	
	hmc_algebraelement2* force_vec = new hmc_algebraelement2[GAUGEMOMENTASIZE2];

	//initial step
	cout << "\tinitial step:" << endl;
	//here, phi is inverted using the orig. gaugefield
	force(parameters, u_out ,
		#ifdef _FERMIONS_
		phi, phi_inv, 
		#endif
		force_vec);
	hmc_float tmp;
	md_update_gauge_momenta(stepsize_half, p_out, force_vec);
	//intermediate steps
	if(steps > 1) cout << "\tperform " << steps << " intermediate steps " << endl;
	for(k = 1; k<steps; k++){
		md_update_gaugefield(stepsize, p_out, u_out);
		force(parameters, u_out ,
			#ifdef _FERMIONS_
			phi, phi_inv, 
			#endif
			force_vec);
		md_update_gauge_momenta(stepsize, p_out, force_vec);
	}
	
	//final step
	cout << "\tfinal step" << endl;
	md_update_gaugefield(stepsize, p_out, u_out);
	force(parameters, u_out ,
		#ifdef _FERMIONS_
		phi, phi_inv, 
		#endif
		force_vec);
	md_update_gauge_momenta(stepsize_half, p_out,force_vec); 
	
	delete [] force_vec;
	
	cout << "\tfinished leapfrog" << endl;
	return HMC_SUCCESS;
}


