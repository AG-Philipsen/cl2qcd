#include "host_operations_gaugemomentum.h"

//gauge-momenta operations

/** @todo memcpy ... */
hmc_error copy_gaugemomenta(hmc_gauge_momentum * source, hmc_gauge_momentum * dest){
	// copies source to destination within cpu memory, layer for momentum array
	return hmc_floatcopy((hmc_complex *)source, (hmc_complex *)dest, GAUGEMOMENTASIZE); // SL: not tested
}

//gaugemomentum is just a hmc_float vector of length GAUGEMOMENTASIZE
hmc_error gaugemomenta_squarenorm(hmc_gauge_momentum * in, hmc_float * result){
	//make sure result is zero
	(*result) = 0.;
	hmc_float sum = 0.;
	for(int i = 0; i<GAUGEMOMENTASIZE; i++){
		sum += ((in)[i])*((in)[i]);
	}
	(*result) = sum;
	return HMC_SUCCESS;
}

/** @todo memset... */
hmc_error set_zero_gaugemomenta(hmc_gauge_momentum * in){
	for(int i = 0; i<GAUGEMOMENTASIZE; i++){
		(in[i]) = 0.;
	}
	return HMC_SUCCESS;
}

hmc_error generate_gaussian_gauge_momenta(hmc_gauge_momentum * out){
	// SL: this is a layer that calls the all-purpose hmc_complex gaussianly-distributed vector
	// with appropriate length and variance, i.e. GAUGEMOMENTASIZE and 1
	// CP: hmc_gauge_momentum should be a real vector, so one should use GAUGEMOMENTASIZE/2 ?!?
	return gaussianComplexVector((hmc_complex *)out, GAUGEMOMENTASIZE/2, 1.0);
	// SL: not yet tested
}
