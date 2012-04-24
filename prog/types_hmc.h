/**
 @file types needed in the HMC-algorithm
*/

#ifndef _TYPES_HMCH_
#define _TYPES_HMCH_

#ifndef _INKERNEL_
//CP: define struct for observables
struct hmc_observables {
	hmc_float plaq;
	hmc_float tplaq;
	hmc_float splaq;
	hmc_complex poly;
	hmc_float deltaH;
	hmc_float prob;
	int accept;
	hmc_float rectangles;
};
#endif

#endif // _TYPES_HMC

