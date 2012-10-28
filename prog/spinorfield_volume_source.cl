__kernel void create_volume_source(__global spinor * const restrict b, __global rngStateStorageType * const restrict rngStates)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

#ifdef _SAME_RND_NUMBERS_
	if(id > 0) return;
	global_size = 1;
#endif

	prng_state rnd;
	prng_loadState(&rnd, rngStates);

	/** @todo What is the correct norm here? */
	hmc_float sigma = 1. / ( VOL4D * 12. );
	spinor out_tmp;
	hmc_complex tmp;

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {
	  /** @todo this might be done more efficient */
	  st_index pos = (id_tmp < VOLSPACE * NTIME / 2) ? get_even_site(id_tmp) : get_odd_site(id_tmp - (VOLSPACE * NTIME / 2));
	  //CP: switch between source content
	  switch(SOURCE_CONTENT){
	  case 1:  //"one"
	    sigma = 1. / ( VOLSPACE * 12. );
	    out_tmp = set_spinor_cold();
	    break;
	   
	  case 2: //"z4"
	    sigma = 1. / ( VOLSPACE * 24. );
	    tmp = Z4_complex_number(&rnd);
	    out_tmp.e0.e0.re = tmp.re;
	    out_tmp.e0.e0.im = tmp.im;
	    tmp = Z4_complex_number(&rnd);
	    out_tmp.e0.e1.re = tmp.re;
	    out_tmp.e0.e1.im = tmp.im;
	    tmp = Z4_complex_number(&rnd);
	    out_tmp.e0.e2.re = tmp.re;
	    out_tmp.e0.e2.im = tmp.im;
	    tmp = Z4_complex_number(&rnd);
	    out_tmp.e1.e0.re = tmp.re;
	    out_tmp.e1.e0.im = tmp.im;
	    tmp = Z4_complex_number(&rnd);
	    out_tmp.e1.e1.re = tmp.re;
	    out_tmp.e1.e1.im = tmp.im;
	    tmp = Z4_complex_number(&rnd);
	    out_tmp.e1.e2.re = tmp.re;
	    out_tmp.e1.e2.im = tmp.im;
	    tmp = Z4_complex_number(&rnd);
	    out_tmp.e2.e0.re = tmp.re;
	    out_tmp.e2.e0.im = tmp.im;
	    tmp = Z4_complex_number(&rnd);
	    out_tmp.e2.e1.re = tmp.re;
	    out_tmp.e2.e1.im = tmp.im;
	    tmp = Z4_complex_number(&rnd);
	    out_tmp.e2.e2.re = tmp.re;
	    out_tmp.e2.e2.im = tmp.im;
	    tmp = Z4_complex_number(&rnd);
	    out_tmp.e3.e0.re = tmp.re;
	    out_tmp.e3.e0.im = tmp.im;
	    tmp = Z4_complex_number(&rnd);
	    out_tmp.e3.e1.re = tmp.re;
	    out_tmp.e3.e1.im = tmp.im;
	    tmp = Z4_complex_number(&rnd);
	    out_tmp.e3.e2.re = tmp.re;
	    out_tmp.e3.e2.im = tmp.im;
	    break;

	  case 3: //"gaussian"
	    /** @todo what is the norm here? */
	    sigma = 1. / ( VOL4D * 24. );
	    tmp = gaussianNormalPair(&rnd);
	    out_tmp.e0.e0.re = tmp.re;
	    out_tmp.e0.e0.im = tmp.im;
	    tmp = gaussianNormalPair(&rnd);
	    out_tmp.e0.e1.re = tmp.re;
	    out_tmp.e0.e1.im = tmp.im;
	    tmp = gaussianNormalPair(&rnd);
	    out_tmp.e0.e2.re = tmp.re;
	    out_tmp.e0.e2.im = tmp.im;
	    tmp = gaussianNormalPair(&rnd);
	    out_tmp.e1.e0.re = tmp.re;
	    out_tmp.e1.e0.im = tmp.im;
	    tmp = gaussianNormalPair(&rnd);
	    out_tmp.e1.e1.re = tmp.re;
	    out_tmp.e1.e1.im = tmp.im;
	    tmp = gaussianNormalPair(&rnd);
	    out_tmp.e1.e2.re = tmp.re;
	    out_tmp.e1.e2.im = tmp.im;
	    tmp = gaussianNormalPair(&rnd);
	    out_tmp.e2.e0.re = tmp.re;
	    out_tmp.e2.e0.im = tmp.im;
	    tmp = gaussianNormalPair(&rnd);
	    out_tmp.e2.e1.re = tmp.re;
	    out_tmp.e2.e1.im = tmp.im;
	    tmp = gaussianNormalPair(&rnd);
	    out_tmp.e2.e2.re = tmp.re;
	    out_tmp.e2.e2.im = tmp.im;
	    tmp = gaussianNormalPair(&rnd);
	    out_tmp.e3.e0.re = tmp.re;
	    out_tmp.e3.e0.im = tmp.im;
	    tmp = gaussianNormalPair(&rnd);
	    out_tmp.e3.e1.re = tmp.re;
	    out_tmp.e3.e1.im = tmp.im;
	    tmp = gaussianNormalPair(&rnd);
	    out_tmp.e3.e2.re = tmp.re;
	    out_tmp.e3.e2.im = tmp.im;
	    break;
	    
	  default:
	    if(id == 0) printf("Problem occured in source kernel: Selected sourcecontent not implemented! Fill with zero...\n");
	    out_tmp = set_spinor_zero();
	  }
	  //multiply by sigma
	  out_tmp = real_multiply_spinor(out_tmp, sqrt(sigma));	  
	  put_spinor_to_field(out_tmp, b, pos.space, pos.time);
	}

	prng_storeState(rngStates, &rnd);
}



