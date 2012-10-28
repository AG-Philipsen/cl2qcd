__kernel void create_timeslice_source(__global spinor * const restrict b, __global rngStateStorageType * const restrict rngStates, int const timeslice)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

#ifdef _SAME_RND_NUMBERS_
	if(id > 0) return;
	global_size = 1;
#endif

	prng_state rnd;
	prng_loadState(&rnd, rngStates);

	hmc_float sigma ;
	spinor out_tmp;

	for(int id_tmp = id; id_tmp < VOLSPACE; id_tmp += global_size) {
	  //CP: switch between source content
	  switch(SOURCE_CONTENT){
	  case 1:  //"one"
	    sigma = 1. / ( VOLSPACE * 12. );
	    out_tmp = set_spinor_cold();
	    break;
	   
	  case 2: //"z4"
	    sigma = 1. / ( VOLSPACE * 24. );
	    hmc_complex tmp;
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
	    
	  default:
	    if(id == 0) printf("Problem occured in source kernel: Selected sourcecontent not implemented! Fill with zero...\n");
	    out_tmp = set_spinor_zero();
	  }
	  //multiply by sigma
	  out_tmp = real_multiply_spinor(out_tmp, sqrt(sigma));	  
	  put_spinor_to_field(out_tmp, b, id_tmp, timeslice);
	}

	prng_storeState(rngStates, &rnd);

	return;
}
