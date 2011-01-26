__kernel void test(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int nsteps,__global hmc_float* check,__global hmc_ocl_gaugefield* gaugefield2, __global hmc_ocl_ran * rnd, __global int * random_field_int, __global float *  random_field_float, __global hmc_float * su2mat, const int size_1, const int size_2, __global hmc_ocl_su3matrix * B , __global hmc_ocl_su3matrix * C){
  
  int id = get_global_id(0);
  hmc_complex testsum;
  testsum.re = 0;
  testsum.im = 0;
  hmc_complex ctmp;
  hmc_ocl_su3matrix prod[SU3SIZE];
hmc_ocl_staplematrix tester[STAPLEMATRIXSIZE];
hmc_ocl_staplematrix tester2[STAPLEMATRIXSIZE];
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

copy_staplematrix(tester, tester2);

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


  //random number test
  int order[3]; 
  for(int i=0;i<size_1/3;i++){
  	random_1_2_3(order, &rnd[id]);
  	random_field_int[3*i] = order[0];
  	random_field_int[3*i + 1] = order[1];
  	random_field_int[3*i + 2] = order[2];
  }
  for(int i=0;i<size_2;i++){
    random_field_float[i] = ocl_new_ran( &rnd[id] );
  }
  hmc_float sample[4];
  
 
  SU2Update(sample, 3.555, &rnd[id]);
  
  su2mat[0] = sample[0];
  su2mat[1] = sample[1];
  su2mat[2] = sample[2];
  su2mat[3] = sample[3];
  
  
  
  //overrelaxing test

/*
  
  hmc_ocl_su3matrix A[SU3SIZE];
    (A[0]).re = 0.042082;
    (A[0]).im = -0.203080;
    (A[1]).re =-0.911624;
    (A[1]).im =-1.019873;
    (A[2]).re =-0.541189;
    (A[2]).im =0.773009;
    (A[3]).re =-0.584805;
    (A[3]).im =-0.370283;
    (A[4]).re =0.066303;
    (A[4]).im =0.222891;
    (A[5]).re =0.239045;
    (A[5]).im =0.477723;
    (A[6]).re =0.435155;
    (A[6]).im =0.627326;
    (A[7]).re =0.042681;
    (A[7]).im =-0.258224;
    (A[8]).re =-0.162545;
    (A[8]).im =0.456923;
  
//    for(int i = 0; i<SU3SIZE; i++){
//      (C[i]).re = (A[i]).re;
//      (C[i]).im = (A[i]).im;
//    }
    
    
       hmc_ocl_su3matrix adj[SU3SIZE];
   copy_su3matrix(adj, A);
   adjoin_su3matrix(adj);
    
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
  hmc_ocl_staplematrix staplematrix[STAPLEMATRIXSIZE];
  
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k;
  //random_1_2_3(order, &rnd[id]);
  order[0] = 1; order[1] = 2; order[2] = 3;
    
  //get_su3matrix(U, gaugefield, pos, t, mu);

//   hmc_complex det = det_su3matrix(A);
//   hmc_complex detadj = complexconj(&det);
//   hmc_complex detsqnorm = complexmult(&det, &detadj);
  //if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(A);
  
  //calc_staple(gaugefield, staplematrix, pos, t, mu);
      for(int i = 0; i<NC; i++){
	for(int j = 0; j<NC; j++){
	  if(i==j) {(staplematrix[i + NC*j]).re = 6.;
	  	(staplematrix[i + NC*j]).im = 0.;}
	  else staplematrix[i + NC*j] = hmc_complex_zero;
      }}
      
          hmc_ocl_su3matrix extW[SU3SIZE]; 
  for(int i=0; i<1; i++)
  {
    multiply_staplematrix(W, A, staplematrix); 
    reduction(w, W, order[i]);

    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );

        hmc_float a[su2_entries];
//         a[0] = w_pauli[0]/k;
//         a[1] = w_pauli[1]/(-k);
//         a[2] = w_pauli[2]/(-k);
//         a[3] = w_pauli[3]/(-k);
	
	a[0] = w_pauli[0];
        a[1] = -w_pauli[1];
        a[2] = -w_pauli[2];
        a[3] = -w_pauli[3];
	
        //Square a and save in w
        w_pauli[0] = a[0]*a[0] - a[1]*a[1] - a[2]*a[2] - a[3]*a[3];
        w_pauli[1] = 2.*a[0]*a[1];
        w_pauli[2] = 2.*a[0]*a[2];
        w_pauli[3] = 2.*a[0]*a[3];
      
        //go back to a su2 matrix in standard basis
        w[0].re = w_pauli[0];
        w[0].im = w_pauli[3];
        w[1].re = w_pauli[2];
        w[1].im = w_pauli[1];
        w[2].re = -w_pauli[2];
        w[2].im = w_pauli[1];
        w[3].re = w_pauli[0];
        w[3].im = -w_pauli[3];
    
//     hmc_ocl_su3matrix extW[SU3SIZE]; 
    extend (extW, order[i], w); 
    //perhaps one can define another acc-func here and save one copying step
    accumulate_su3matrix_prod(extW, A);
    copy_su3matrix(A, extW);
  }
  //put_su3matrix(gaugefield, A, pos, t, mu);
  
   for(int i = 0; i<SU3SIZE; i++){
     (B[i]).re = (A[i]).re;
     (B[i]).im = (A[i]).im;
   }
  
      for(int i = 0; i<SU3SIZE; i++){
     (C[i]).re = (adj[i]).re;
     (C[i]).im = (adj[i]).im;
   }
*/
  return;
}
