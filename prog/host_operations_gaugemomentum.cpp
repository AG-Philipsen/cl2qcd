#include "host_operations_gaugemomentum.h"

//gauge-momenta operations

/** @todo memcpy ... */
/** @todo get rid of the workaround.. */
// hmc_error copy_gaugemomenta(hmc_gauge_momentum * source, hmc_gauge_momentum * dest){
// 	// copies source to destination within cpu memory, layer for momentum array
// 	return hmc_floatcopy((hmc_float *)source, (hmc_float *)dest, GAUGEMOMENTASIZE); // SL: not tested
// }

hmc_error copy_gaugemomenta(hmc_algebraelement2 * source, hmc_algebraelement2 * dest){
	for(int i = 0; i<GAUGEMOMENTASIZE2; i++){	
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
	
//gaugemomentum is just a hmc_float vector of length GAUGEMOMENTASIZE
hmc_error gaugemomenta_squarenorm(hmc_algebraelement2 * in, hmc_float * result){
	//make sure result is zero
	(*result) = 0.;
	hmc_float sum = 0.;
	for(int i = 0; i<GAUGEMOMENTASIZE2; i++){
		sum += (in[i]).e0*(in[i]).e0 +
		       (in[i]).e1*(in[i]).e1 +
		       (in[i]).e2*(in[i]).e2 +
		       (in[i]).e3*(in[i]).e3 +
		       (in[i]).e4*(in[i]).e4 +
		       (in[i]).e5*(in[i]).e5 +
		       (in[i]).e6*(in[i]).e6 +
		       (in[i]).e7*(in[i]).e7;
	}
	(*result) = sum;
	return HMC_SUCCESS;
}

/** @todo memset... */
hmc_error set_zero_gaugemomenta(hmc_algebraelement2 * in){
	for(int i = 0; i<GAUGEMOMENTASIZE2; i++){
		(in[i]).e0 = 0.;
		(in[i]).e1 = 0.;
		(in[i]).e2 = 0.;
		(in[i]).e3 = 0.;
		(in[i]).e4 = 0.;
		(in[i]).e5 = 0.;
		(in[i]).e6 = 0.;
		(in[i]).e7 = 0.;
	}
	return HMC_SUCCESS;
}

/** @todo rewrite gaussianComplexVector to handle structs??*/
hmc_error generate_gaussian_gauge_momenta(hmc_algebraelement2 * out){
	// SL: this is a layer that calls the all-purpose hmc_complex gaussianly-distributed vector
	// with appropriate length and variance, i.e. GAUGEMOMENTASIZE and 1
	// CP: hmc_gauge_momentum should be a real vector, so one should use GAUGEMOMENTASIZE/2 ?!?
	// CP: workaround for structs:
	hmc_error err = 0;
	hmc_float * tmp = new hmc_float[GAUGEMOMENTASIZE2*8];
	err =  gaussianComplexVector((hmc_complex *)tmp, GAUGEMOMENTASIZE/2, 1.0);
	// SL: not yet tested
	for(int i = 0; i<GAUGEMOMENTASIZE2; i++){	
		(out[i]).e0 = tmp[i*8 + 0];
		(out[i]).e1 = tmp[i*8 + 1];
		(out[i]).e2 = tmp[i*8 + 2];
		(out[i]).e3 = tmp[i*8 + 3];
		(out[i]).e4 = tmp[i*8 + 4];
		(out[i]).e5 = tmp[i*8 + 5];
		(out[i]).e6 = tmp[i*8 + 6];
		(out[i]).e7 = tmp[i*8 + 7];
	}
	
	
	delete [] tmp;
	return err;
}

void acc_factor_times_algebraelement(hmc_algebraelement2 * inout, hmc_float factor, hmc_algebraelement2 force_in){
	(*inout).e0+=factor*(force_in).e0; 
	(*inout).e1+=factor*(force_in).e1;
	(*inout).e2+=factor*(force_in).e2;
	(*inout).e3+=factor*(force_in).e3;
	(*inout).e4+=factor*(force_in).e4;
	(*inout).e5+=factor*(force_in).e5;
	(*inout).e6+=factor*(force_in).e6;
	(*inout).e7+=factor*(force_in).e7;
}