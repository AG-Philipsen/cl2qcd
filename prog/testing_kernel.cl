__kernel void test(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int nsteps,__global hmc_float* check,__global hmc_ocl_gaugefield* gaugefield2){
  int id = get_global_id(0);
  hmc_complex testsum;
  testsum.re = 0;
  testsum.im = 0;
  hmc_complex ctmp;
  hmc_ocl_su3matrix prod[SU3SIZE];
  unit_su3matrix(prod);
  if(id==0) {
  for(int spacepos = 0; spacepos < VOLSPACE; spacepos++) {
    for(int mu = 0; mu<NDIM; mu++) {
      for(int t=0; t<NTIME; t++) {
	hmc_ocl_su3matrix tmp[SU3SIZE];
	get_su3matrix(tmp,gaugefield,spacepos,t,mu);
	hmc_ocl_su3matrix tmp2[SU3SIZE];
	copy_su3matrix(tmp2,tmp);

	//test matrix operations...
	hmc_ocl_su3matrix tmp3[SU3SIZE];
	multiply_su3matrices(tmp3,tmp2,tmp);
	adjoin_su3matrix(tmp3);
	accumulate_su3matrix_prod(prod, tmp3);
	ctmp = det_su3matrix(tmp3);
	ctmp.re -= hmc_one_f;
	hmc_complex ctmpconj = complexconj(&ctmp);
	hmc_complex square = complexmult(&ctmp,&ctmpconj);
	complexaccumulate(&testsum,&square);

	//	unit_su3matrix(tmp2);
	//	zero_su3matrix(tmp2);

	/*
	//test trace
	hmc_complex trace = trace_su3matrix(tmp2);
	trace.re -= 3.;
	complexaccumulate(&testsum,&trace);
	*/

	put_su3matrix(gaugefield,tmp2,spacepos,t,mu);
	/*	for(int n=0; n<2*SU3SIZE*VOLSPACE*NDIM*NTIME) {
	  gaugefield[n] = gf[n];
	  }*/
      }
    }
  }
  adjoin_su3(gaugefield,gaugefield2);
  adjoin_su3(gaugefield2,gaugefield);
  /*
  hmc_complex myctmp = global_trace_su3(gaugefield,2);
  testsum.re = myctmp.re - 3*VOLSPACE*NTIME;
  testsum.im = myctmp.im;
  */

  *check = testsum.re*testsum.re + testsum.im*testsum.im;
  *check += det_su3matrix(prod).re - hmc_one_f;
  }
  return;
}
