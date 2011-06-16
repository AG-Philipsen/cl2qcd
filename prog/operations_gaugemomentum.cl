/**
 * @file operations on algebraelements
 */

ae set_zero_ae(){
	ae tmp;
	tmp.e0 = 0.;
	tmp.e1 = 0.;
	tmp.e2 = 0.;
	tmp.e3 = 0.;
	tmp.e4 = 0.;
	tmp.e5 = 0.;
	tmp.e6 = 0.;
	tmp.e7 = 0.;
	return tmp;
}

hmc_float ae_squarenorm(ae in){
	hmc_float result = 
		(in).e0*(in).e0 +
		(in).e1*(in).e1 +
		(in).e2*(in).e2 +
		(in).e3*(in).e3 +
		(in).e4*(in).e4 +
		(in).e5*(in).e5 +
		(in).e6*(in).e6 +
		(in).e7*(in).e7;
	return result;
}

ae acc_factor_times_algebraelement(ae in, hmc_float factor, ae force_in){
	ae tmp;
	tmp.e0 = (in).e0+factor*(force_in).e0; 
	tmp.e1 = (in).e1+factor*(force_in).e1;
	tmp.e2 = (in).e2+factor*(force_in).e2;
	tmp.e3 = (in).e3+factor*(force_in).e3;
	tmp.e4 = (in).e4+factor*(force_in).e4;
	tmp.e5 = (in).e5+factor*(force_in).e5;
	tmp.e6 = (in).e6+factor*(force_in).e6;
	tmp.e7 = (in).e7+factor*(force_in).e7;
	return tmp;
}



/**
 * @file operations used by gaugemomentum
 */


/** @todo memcpy ... */
/** @todo get rid of the workaround.. */
// hmc_error copy_gaugemomenta(hmc_gauge_momentum * source, hmc_gauge_momentum * dest){
// 	// copies source to destination within cpu memory, layer for momentum array
// 	return hmc_floatcopy((hmc_float *)source, (hmc_float *)dest, GAUGEMOMENTASIZE); // SL: not tested
// }

//deprecated, this can be done by in = out in the code...
/*
hmc_algebraelement copy_gaugemomenta(hmc_algebraelement2 * source, hmc_algebraelement2 * dest){
	hmc_algebraelement2 tmp;
	for(int i = 0; i<GAUGEMOMENTASIZE; i++){	
		(dest[i]).e0 = (source[i]).e0;
		(dest[i]).e1 = (source[i]).e1;
		(dest[i]).e2 = (source[i]).e2;
		(dest[i]).e3 = (source[i]).e3;
		(dest[i]).e4 = (source[i]).e4;
		(dest[i]).e5 = (source[i]).e5;
		(dest[i]).e6 = (source[i]).e6;
		(dest[i]).e7 = (source[i]).e7;
	}
	return HMC_SUCCESS;
}
*/
	
/** @todo add args for reduction... */
__kernel void gaugemomenta_squarenorm(__global ae * in){
	int id = get_global_id(0);
	if(id == 0){
		hmc_float result = 0.;
		for(int i = 0; i<GAUGEMOMENTASIZE; i++){
			result += ae_squarenorm(in[i]);
		}
		
		/** @todo add reduction.. */
		
	}
}

/** @todo memset... */
__kernel void set_zero_gaugemomentum(__global ae * in){
	int id = get_global_id(0);
	if(id == 0){
		for(int i = 0; i<GAUGEMOMENTASIZE; i++){
			in[i] = set_zero_ae();
		}
	}
}

/** @todo rewrite gaussianComplexVector to handle structs??*/
__kernel void generate_gaussian_gauge_momenta(__global ae * out){
	int id = get_global_id(0);
	if(id == 0){
		// SL: this is a layer that calls the all-purpose hmc_complex gaussianly-distributed vector
		// with appropriate length and variance, i.e. GAUGEMOMENTASIZE and 1
		// CP: hmc_gauge_momentum should be a real vector, so one should use GAUGEMOMENTASIZE/2 ?!?
		// CP: workaround for structs:
		//hmc_float * tmp = new hmc_float[GAUGEMOMENTASIZE*8];
		//gaussianComplexVector((hmc_complex *)tmp, GAUGEMOMENTASIZE/2, 1.0);
		// SL: not yet tested
// 		for(int i = 0; i<GAUGEMOMENTASIZE; i++){	
// 			(out[i]).e0 = tmp[i*8 + 0];
// 			(out[i]).e1 = tmp[i*8 + 1];
// 			(out[i]).e2 = tmp[i*8 + 2];
// 			(out[i]).e3 = tmp[i*8 + 3];
// 			(out[i]).e4 = tmp[i*8 + 4];
// 			(out[i]).e5 = tmp[i*8 + 5];
// 			(out[i]).e6 = tmp[i*8 + 6];
// 			(out[i]).e7 = tmp[i*8 + 7];
// 		}
// 		delete [] tmp;
	}
}



