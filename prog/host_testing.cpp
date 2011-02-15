#include "host_testing.h"


void testing_correlator(hmc_gaugefield* gf, inputparameters* parameters) {
  simple_correlator(gf, (*parameters).get_kappa(),(*parameters).get_mu(), (*parameters).get_theta_fermion(), (*parameters).get_cgmax());
 return;
}

void testing_fermionmatrix(){
  hmc_color_vector* vec = new hmc_color_vector[NC];
  hmc_color_vector* out = new hmc_color_vector[NC];

  hmc_su3matrix u;
  unit_su3matrix(&u);

  vec[0].re=1;
  vec[0].im=2;
  vec[1].re=3;
  vec[1].im=4;
  vec[2].re=5;
  vec[2].im=6;

  su3matrix_times_colorvector(&u,vec,out);

  for(int a=0; a<NC; a++) {
    printf("(%f,%f)\n",out[a].re,out[a].im);
  }

  delete [] vec;
  delete [] out;
  return;
}

void testing_eoprec_spinor() {
  hmc_spinor_field* test = new hmc_spinor_field[SPINORFIELDSIZE];

  hmc_eoprec_spinor_field* eventest = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* oddtest = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];

  printf("global squarenorm: %f\n",global_squarenorm(test));
  convert_to_eoprec(eventest,oddtest,test);
  printf("eo global squarenorm: %f -- %f\n",global_squarenorm_eoprec(eventest),global_squarenorm_eoprec(oddtest));
  convert_from_eoprec(eventest,oddtest,test);
  printf("global squarenorm: %f\n",global_squarenorm(test));

  set_zero_spinor(test);

  printf("global squarenorm: %f\n",global_squarenorm(test));
  convert_to_eoprec(eventest,oddtest,test);
  printf("eo global squarenorm: %f -- %f\n",global_squarenorm_eoprec(eventest),global_squarenorm_eoprec(oddtest));
  convert_from_eoprec(eventest,oddtest,test);
  printf("global squarenorm: %f\n",global_squarenorm(test));

  fill_with_one(test, 0, 0, 0, 0);

  printf("global squarenorm: %f\n",global_squarenorm(test));
  convert_to_eoprec(eventest,oddtest,test);
  printf("eo global squarenorm: %f -- %f\n",global_squarenorm_eoprec(eventest),global_squarenorm_eoprec(oddtest));
  convert_from_eoprec(eventest,oddtest,test);
  printf("global squarenorm: %f\n",global_squarenorm(test));

  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int a=0; a<NSPIN; a++) {
	for(int j=0; j<NC; j++) {
	  fill_with_one(test, n, t, a, j);
	}
      }
    }
  }


  printf("global squarenorm: %f\n",global_squarenorm(test));
  convert_to_eoprec(eventest,oddtest,test);
  printf("eo global squarenorm: %f -- %f\n",global_squarenorm_eoprec(eventest),global_squarenorm_eoprec(oddtest));
  convert_from_eoprec(eventest,oddtest,test);
  printf("global squarenorm: %f\n",global_squarenorm(test));

  set_zero_spinor(test);

  printf("global squarenorm: %f\n",global_squarenorm(test));
  convert_to_eoprec(eventest,oddtest,test);
  printf("eo global squarenorm: %f -- %f\n",global_squarenorm_eoprec(eventest),global_squarenorm_eoprec(oddtest));
  convert_from_eoprec(eventest,oddtest,test);
  printf("global squarenorm: %f\n",global_squarenorm(test));

  printf("fill only even sites:\n");
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int a=0; a<NSPIN; a++) {
	for(int j=0; j<NC; j++) {
	  int check = t + get_spacecoord(n,1) + get_spacecoord(n,2) + get_spacecoord(n,3);
	  if(check%2 == 0) fill_with_one(test, n, t, a, j);
	}
      }
    }
  }


  printf("global squarenorm: %f\n",global_squarenorm(test));
  convert_to_eoprec(eventest,oddtest,test);
  printf("eo global squarenorm: %f -- %f\n",global_squarenorm_eoprec(eventest),global_squarenorm_eoprec(oddtest));
  convert_from_eoprec(eventest,oddtest,test);
  printf("global squarenorm: %f\n",global_squarenorm(test));

  set_zero_spinor(test);

  printf("fill only odd sites:\n");
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int a=0; a<NSPIN; a++) {
	for(int j=0; j<NC; j++) {
	  int check = t + get_spacecoord(n,1) + get_spacecoord(n,2) + get_spacecoord(n,3);
	  if(check%2 == 1) fill_with_one(test, n, t, a, j);
	}
      }
    }
  }


  printf("global squarenorm: %f\n",global_squarenorm(test));
  convert_to_eoprec(eventest,oddtest,test);
  printf("eo global squarenorm: %f -- %f\n",global_squarenorm_eoprec(eventest),global_squarenorm_eoprec(oddtest));
  convert_from_eoprec(eventest,oddtest,test);
  printf("global squarenorm: %f\n",global_squarenorm(test));


  printf("global squarenorm: %f\n",global_squarenorm(test));
  convert_to_eoprec(eventest,oddtest,test);
  printf("eo global squarenorm: %f -- %f\n",global_squarenorm_eoprec(eventest),global_squarenorm_eoprec(oddtest));
  convert_to_kappa_format(oddtest,4.214523);
  printf("eo global squarenorm: %f -- %f\n",global_squarenorm_eoprec(eventest),global_squarenorm_eoprec(oddtest));
  convert_from_kappa_format(oddtest, oddtest,4.214523);
  printf("eo global squarenorm: %f -- %f\n",global_squarenorm_eoprec(eventest),global_squarenorm_eoprec(oddtest));
  convert_from_eoprec(eventest,oddtest,test);
  printf("global squarenorm: %f\n",global_squarenorm(test));

  delete [] test;
  delete [] eventest;
  delete [] oddtest;

  return;
}

void print_fullspinorfield(hmc_spinor* in){
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int a=0; a<NSPIN; a++) {
	printf("[");
	for(int j=0; j<NC; j++) {
	  printf("(%.3f,%.3f) ", (*in).re, (*in).im );
	}
      }printf("]\n");
    }printf("\n");
  }printf("\n");
}

void testing_spinor() {
  hmc_spinor_field* test = new hmc_spinor_field[SPINORFIELDSIZE];


  printf("global squarenorm: %f\n",global_squarenorm(test));
  set_zero_spinor(test);


  printf("global squarenorm: %f\n",global_squarenorm(test));

  fill_with_one(test, 0, 0, 0, 0);

  printf("global squarenorm: %f\n",global_squarenorm(test));


  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int a=0; a<NSPIN; a++) {
	for(int j=0; j<NC; j++) {
	  fill_with_one(test, n, t, a, j);
	}
      }
    }
  }


  printf("global squarenorm: %f\n",global_squarenorm(test));

  delete [] test;

  return;
}

void print_su3mat(hmc_su3matrix* A){
#ifdef _RECONSTRUCT_TWELVE_
  printf("\n| (%f,%f)\t(%f,%f)\t(%f,%f) |\n",(*A)[0].re,(*A)[0].im,(*A)[2].re,(*A)[2].im,(*A)[4].re,(*A)[4].im);
  printf("| (%f,%f)\t(%f,%f)\t(%f,%f) |\n",(*A)[1].re,(*A)[1].im,(*A)[3].re,(*A)[3].im,(*A)[5].re,(*A)[5].im);
  hmc_complex ca = reconstruct_su3(A,0);
  hmc_complex cb = reconstruct_su3(A,1);
  hmc_complex cc = reconstruct_su3(A,2);
  printf("| (%f,%f)\t(%f,%f)\t(%f,%f) |\n",ca.re,ca.im,cb.re,cb.im,cc.re,cc.im);
#else
  printf("\n");
  for(int a = 0; a<NC; a++) {
  printf("| (%f,%f)\t(%f,%f)\t(%f,%f) |\n",(*A)[a][0].re,(*A)[a][0].im,(*A)[a][1].re,(*A)[a][1].im,(*A)[a][2].re,(*A)[a][2].im);}
  printf("\n");
#endif
  return;
}

void print_staplemat(hmc_staplematrix* A){
#ifdef _RECONSTRUCT_TWELVE_
  printf("\n| (%f,%f)\t(%f,%f)\t(%f,%f) |\n",(*A)[0].re,(*A)[0].im,(*A)[2].re,(*A)[2].im,(*A)[4].re,(*A)[4].im);
  printf("| (%f,%f)\t(%f,%f)\t(%f,%f) |\n",(*A)[1].re,(*A)[1].im,(*A)[3].re,(*A)[3].im,(*A)[5].re,(*A)[5].im);
  printf("| (%f,%f)\t(%f,%f)\t(%f,%f) |\n",(*A)[6].re,(*A)[6].im,(*A)[7].re,(*A)[7].im,(*A)[8].re,(*A)[8].im);
#else
  printf("\n");
  for(int a = 0; a<NC; a++) 
  printf("| (%f,%f)\t(%f,%f)\t(%f,%f) |\n",(*A)[a][0].re,(*A)[a][0].im,(*A)[a][1].re,(*A)[a][1].im,(*A)[a][2].re,(*A)[a][2].im);
  printf("\n");
#endif
  return;
}

void testing_su3mat(){
  hmc_su3matrix A;
#ifdef _RECONSTRUCT_TWELVE_
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC-1; b++) {
      int n = b + a*(NC-1);
      A[n].re = b;
      A[n].im = a;
    }
  }
  print_su3mat(&A);

  hmc_complex trace = trace_su3matrix(&A);
  printf("trace: (%f,%f)\n",trace.re,trace.im);

  printf("adjoint:\n");
  adjoin_su3matrix(&A);
  print_su3mat(&A);
  trace = trace_su3matrix(&A);
  printf("trace: (%f,%f)\n",trace.re,trace.im);

  unit_su3matrix(&A);
  print_su3mat(&A);
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      A[a][b].re = a;
      A[a][b].im = b;
    }
  }
  print_su3mat(&A);

  hmc_complex trace = trace_su3matrix(&A);
  printf("trace: (%f,%f)\n",trace.re,trace.im);

  printf("adjoint:\n");
  adjoin_su3matrix(&A);
  print_su3mat(&A);
  trace = trace_su3matrix(&A);
  printf("trace: (%f,%f)\n",trace.re,trace.im);

  unit_su3matrix(&A);
  print_su3mat(&A);
#endif
return;
}

void testing_gaugefield(){
  hmc_gaugefield test;

  hmc_complex trace;

  for(int mu=0; mu<NDIM; mu++) {
    trace = global_trace_su3(&test,mu);
    printf("global trace mu = %d:  (%f,%f)\n",mu,trace.re,trace.im);
  }
  set_gaugefield_cold(&test);
  for(int mu=0; mu<NDIM; mu++) {
    trace = global_trace_su3(&test,mu);
    printf("global trace mu = %d:  (%f,%f)\n",mu,trace.re,trace.im);
  }
  
  hmc_gaugefield out;
  adjoin_su3(&test,&out);
  for(int mu=0; mu<NDIM; mu++) {
    trace = global_trace_su3(&out,mu);
    printf("global trace mu = %d:  (%f,%f)\n",mu,trace.re,trace.im);
  }
  return;
}

void testing_geometry(){
  int coord[NDIM];
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      printf("%d\t%d\n",t,n);
      coord[0]=t;
      coord[1]=get_spacecoord(n,1);
      coord[2]=get_spacecoord(n,2);
      coord[3]=get_spacecoord(n,3);
      printf("%d\t%d\t%d\t%d\n",coord[0],coord[1],coord[2],coord[3]);
      int nspace = get_nspace(coord);
      printf("%d\t%d\n",t,nspace);
      printf("\n");
    }
  }

  int nspace=3;
  int neighbor;
  printf("%d\t%d\t%d\t%d\n",0,get_spacecoord(nspace,1),get_spacecoord(nspace,2),get_spacecoord(nspace,3));
  neighbor = get_neighbor(nspace,1);
  printf("%d\t%d\t%d\t%d\n",0,get_spacecoord(neighbor,1),get_spacecoord(neighbor,2),get_spacecoord(neighbor,3));
  neighbor = get_neighbor(nspace,2);
  printf("%d\t%d\t%d\t%d\n",0,get_spacecoord(neighbor,1),get_spacecoord(neighbor,2),get_spacecoord(neighbor,3));
  neighbor = get_neighbor(nspace,3);
  printf("%d\t%d\t%d\t%d\n",0,get_spacecoord(neighbor,1),get_spacecoord(neighbor,2),get_spacecoord(neighbor,3));

  nspace=1;
  printf("%d\t%d\t%d\t%d\n",0,get_spacecoord(nspace,1),get_spacecoord(nspace,2),get_spacecoord(nspace,3));
  neighbor = get_neighbor(nspace,1);
  printf("%d\t%d\t%d\t%d\n",0,get_spacecoord(neighbor,1),get_spacecoord(neighbor,2),get_spacecoord(neighbor,3));
  neighbor = get_neighbor(nspace,2);
  printf("%d\t%d\t%d\t%d\n",0,get_spacecoord(neighbor,1),get_spacecoord(neighbor,2),get_spacecoord(neighbor,3));
  neighbor = get_neighbor(nspace,3);
  printf("%d\t%d\t%d\t%d\n",0,get_spacecoord(neighbor,1),get_spacecoord(neighbor,2),get_spacecoord(neighbor,3));

  nspace=27;
  printf("%d\t%d\t%d\t%d\n",0,get_spacecoord(nspace,1),get_spacecoord(nspace,2),get_spacecoord(nspace,3));
  neighbor = get_neighbor(nspace,1);
  printf("%d\t%d\t%d\t%d\n",0,get_spacecoord(neighbor,1),get_spacecoord(neighbor,2),get_spacecoord(neighbor,3));
  neighbor = get_neighbor(nspace,2);
  printf("%d\t%d\t%d\t%d\n",0,get_spacecoord(neighbor,1),get_spacecoord(neighbor,2),get_spacecoord(neighbor,3));
  neighbor = get_neighbor(nspace,3);
  printf("%d\t%d\t%d\t%d\n",0,get_spacecoord(neighbor,1),get_spacecoord(neighbor,2),get_spacecoord(neighbor,3));

  return;
}


void testing_su3matrix(hmc_gaugefield * in, int spacepos, int timepos){
  hmc_su3matrix test;
  hmc_complex trace, det;

  for(int mu=0; mu<NDIM; mu++) {
    get_su3matrix(&test, in, spacepos, timepos, mu);
    trace = trace_su3matrix(&test);
    //printf("trace mu = %d:  (%f,%f)\n",mu,trace.re,trace.im);
  }

  for(int mu=0; mu<NDIM; mu++) {
    get_su3matrix(&test, in, spacepos, timepos, mu);
    det = det_su3matrix(&test);
    if(det.re < 0.99) printf("det mu = %d:  (%f,%f)\n",mu,det.re,det.im);
    if(det.re > 1.01) printf("det mu = %d:  (%f,%f)\n",mu,det.re,det.im);
  }

  return;
}

void testing_matrix_ops(hmc_gaugefield * in){
  int spacepos = 23, timepos = 3;
  hmc_su3matrix test, test2, unity, naiv;
printf("\n\ntestmatrix:\n");
  get_su3matrix(&test, in, spacepos, timepos, 2); 
print_su3mat(&test);
printf("\nadjoined:\n");
  get_su3matrix(&test2, in, spacepos, timepos, 2);
  adjoin_su3matrix(&test2);
print_su3mat(&test2);
printf("\n multiplication: mat * adj(mat)\n");
    multiply_su3matrices(&naiv, &test, &test2);
    print_su3mat(&naiv);

printf("\naccumulate mat and adj(mat) \n");
  accumulate_su3matrix_prod(&test2, &test);
print_su3mat(&test2);

printf("\ncopy that into some other matrix\n");
  copy_su3matrix( &unity, &test2);
  print_su3mat(&unity);
  
printf("\n");
  return;
}

void testing_adjoin(hmc_gaugefield * in, int spacepos, int timepos){
  hmc_su3matrix test, test2, unity;
  hmc_complex trace, det;

  for(int mu=0; mu<NDIM; mu++) {
    get_su3matrix(&test, in, spacepos, timepos, mu);
    get_su3matrix(&test2, in, spacepos, timepos, mu);
    adjoin_su3matrix(&test2);
    trace = trace_su3matrix(&test2);
    det = det_su3matrix(&test2);
    printf("trace mu = %d:  (%f,%f)\n",mu,trace.re,trace.im);
    printf("det mu = %d:  (%f,%f)\n",mu,det.re,det.im);
    multiply_su3matrices(&unity, &test, &test2);
    print_su3mat(&unity);
  }

  return;
}


void testing_det_global(hmc_gaugefield * in){
  int cter = 0, cter2 = 0;
  hmc_complex det;
  hmc_su3matrix test;
  for (int t = 0; t < NTIME; t++){
    for (int pos = 0; pos< VOLSPACE; pos++){
     for (int mu = 0; mu < NDIM; mu++){
       get_su3matrix(&test, in, pos, t, mu);
       det = det_su3matrix(&test);
       if(det.re < 0.9999999) cter++;
       if(det.re > 1.0000001) cter++;
       if(det.im < -0.0000001) cter++;
       if(det.im > 0.0000001) cter++;
       adjoin_su3matrix(&test);
       det = det_su3matrix(&test);
       if(det.re < 0.9999999) {cter2++;}
       if(det.re > 1.0000001) {cter2++;}
       if(det.im < -0.0000001) {cter2++;}
       if(det.im > 0.0000001) {cter2++;}
  }}}
  printf("there were %i wrong matrices and %i wrong adjoint matrices\n" ,cter, cter2);
return;
}

void print_linkplusstaplematrix(hmc_gaugefield * in, int pos, int t, int mu){
  hmc_su3matrix dummy;
  get_su3matrix(&dummy, in, pos,t,mu);

  printf("link at %i %i %i:\n", pos, t, mu);
  print_su3mat(&dummy);

  hmc_staplematrix dummy2;
  calc_staple(in, &dummy2, pos, t,  mu);
  printf("staplematrix at %i %i %i:\n", pos, t, mu);
  print_staplemat(&dummy2);  
  
  return;
}


void testing_heatbath_norandommat_no123(hmc_su3matrix * in, hmc_staplematrix * staple_in, hmc_su3matrix * out, hmc_float beta){
  hmc_staplematrix W;
  hmc_su3matrix su3_in;
  copy_su3matrix(&su3_in, in);
      
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
//   random_1_2_3(order);
  order[0] = 1; order[1] = 2; order[2] = 3;
  
  hmc_complex det = det_su3matrix(&su3_in);
      hmc_complex detadj = complexconj(&det);
      hmc_complex detsqnorm = complexmult(&det, &detadj);
      if( (detsqnorm.re - hmc_one_f) <= projectioneps)
	project_su3(&su3_in);
      
  for(int i=0; i<3; i++)
  {
    //Produkt aus U und staple, saved in W
    multiply_staplematrix(&W, &su3_in, staple_in);
    //Reduktion auf SU(2) submatrix
    reduction (w, W, order[i]);
    //Darstellung von W in Paulibasis
    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);

    //Berechnung der Normierung k, entspricht Determinante in normaler Basis
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
    w_pauli[0] = hmc_float(1)/k * w_pauli[0];
    w_pauli[1] = hmc_float(-1)/k * w_pauli[1];
    w_pauli[2] = hmc_float(-1)/k * w_pauli[2];
    w_pauli[3] = hmc_float(-1)/k * w_pauli[3];
	
    //beta' = 2Beta/Nc*k
    hmc_float beta_neu =  2.*beta / hmc_float(NC)*k;
      
    //Neuer Link in Paulibasis mittels Kennedy-Pendleton-Algorithmus
    //SU2Update(r_pauli, beta_neu);
    //random su2-matrix
    r_pauli[0] = 0.748102;
    r_pauli[1] = 0.293055;
    r_pauli[2] = 0.591048;
    r_pauli[3] = -0.0715786;

    //Multipliziere neuen Link mit w, alles in Paulibasis
    hmc_float su2_tmp[su2_entries];
    su2_tmp[0] = r_pauli[0]*w_pauli[0] - r_pauli[1]*w_pauli[1] - r_pauli[2]*w_pauli[2] - r_pauli[3]*w_pauli[3] ;
    su2_tmp[1] = w_pauli[0]*r_pauli[1] + w_pauli[1]*r_pauli[0] - r_pauli[2]*w_pauli[3] + r_pauli[3]*w_pauli[2] ;
    su2_tmp[2] = w_pauli[0]*r_pauli[2] + w_pauli[2]*r_pauli[0] - r_pauli[3]*w_pauli[1] + r_pauli[1]*w_pauli[3] ;
    su2_tmp[3] = w_pauli[0]*r_pauli[3] + w_pauli[3]*r_pauli[0] - r_pauli[1]*w_pauli[2] + r_pauli[2]*w_pauli[1] ;
    r_pauli[0] = su2_tmp[0];      
    r_pauli[1] = su2_tmp[1];
    r_pauli[2] = su2_tmp[2];
    r_pauli[3] = su2_tmp[3];

      
    //go back to a su2 matrix in standard basis
    w[0].re = r_pauli[0];
    w[0].im = r_pauli[3];
    w[1].re = r_pauli[2];
    w[1].im = r_pauli[1];
    w[2].re = -r_pauli[2];
    w[2].im = r_pauli[1];
    w[3].re = r_pauli[0];
    w[3].im = -r_pauli[3];
 
    //extend to SU3
    hmc_su3matrix extW, tmp;
    extend (&extW, order[i], w);
	
    multiply_su3matrices(&tmp, &extW, &su3_in);
    copy_su3matrix(&su3_in, &tmp);
  }
  copy_su3matrix(out, &su3_in);
  
  return;
}

void testing_heatbath_no123(hmc_su3matrix * in, hmc_staplematrix * staple_in, hmc_su3matrix * out, hmc_float beta, hmc_float * rnd_array, int * cter_out){
  hmc_staplematrix W;
  hmc_su3matrix su3_in;
  copy_su3matrix(&su3_in, in);
      
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
//   random_1_2_3(order);
  order[0] = 1; order[1] = 2; order[2] = 3;
  
  hmc_complex det = det_su3matrix(&su3_in);
      hmc_complex detadj = complexconj(&det);
      hmc_complex detsqnorm = complexmult(&det, &detadj);
      if( (detsqnorm.re - hmc_one_f) <= projectioneps)
	project_su3(&su3_in);
      
  int cter = 0;
  for(int i=0; i<3; i++)
  {
    //Produkt aus U und staple, saved in W
    multiply_staplematrix(&W, &su3_in, staple_in);
    //Reduktion auf SU(2) submatrix
    reduction (w, W, order[i]);
    //Darstellung von W in Paulibasis
    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);

    //Berechnung der Normierung k, entspricht Determinante in normaler Basis
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
    w_pauli[0] = hmc_float(1)/k * w_pauli[0];
    w_pauli[1] = hmc_float(-1)/k * w_pauli[1];
    w_pauli[2] = hmc_float(-1)/k * w_pauli[2];
    w_pauli[3] = hmc_float(-1)/k * w_pauli[3];
	
    //beta' = 2Beta/Nc*k
    hmc_float beta_neu =  2.*beta / hmc_float(NC)*k;
    
    
    //Neuer Link in Paulibasis mittels Kennedy-Pendleton-Algorithmus
    //function  SU2Update(r_pauli, beta_neu):
    hmc_float delta;
    hmc_float a0 ;
    hmc_float eta ;
    hmc_float alpha = beta_neu;

    do
    {
      hmc_float rnd1 = rnd_array[cter];
      hmc_float rnd2 = rnd_array[cter+1];
      hmc_float rnd3 = rnd_array[cter+2];
      hmc_float rnd4 = rnd_array[cter+3];
      delta = -log(rnd1)/alpha*pow(cos(2. * PI * rnd2), 2.) -log(rnd3)/alpha;
      a0 = 1.-delta;
      eta = rnd4;
      cter += 4;
    }while ( (1.-0.5*delta) < eta*eta);
    hmc_float rnd5 = rnd_array[cter];
    hmc_float rnd6 = rnd_array[cter+1];
    cter += 2;
    hmc_float phi = 2.*PI*rnd5;
    hmc_float theta = asin(2.*rnd6 - 1.);
    r_pauli[0] = a0;
    r_pauli[1] = sqrt(1.-a0 * a0)*cos(theta) * cos(phi);
    r_pauli[2] = sqrt(1.-a0 * a0)*cos(theta) * sin(phi);
    r_pauli[3] = sqrt(1.-a0 * a0)*sin(theta);
    
    
    //Multipliziere neuen Link mit w, alles in Paulibasis
    hmc_float su2_tmp[su2_entries];
    su2_tmp[0] = r_pauli[0]*w_pauli[0] - r_pauli[1]*w_pauli[1] - r_pauli[2]*w_pauli[2] - r_pauli[3]*w_pauli[3] ;
    su2_tmp[1] = w_pauli[0]*r_pauli[1] + w_pauli[1]*r_pauli[0] - r_pauli[2]*w_pauli[3] + r_pauli[3]*w_pauli[2] ;
    su2_tmp[2] = w_pauli[0]*r_pauli[2] + w_pauli[2]*r_pauli[0] - r_pauli[3]*w_pauli[1] + r_pauli[1]*w_pauli[3] ;
    su2_tmp[3] = w_pauli[0]*r_pauli[3] + w_pauli[3]*r_pauli[0] - r_pauli[1]*w_pauli[2] + r_pauli[2]*w_pauli[1] ;
    r_pauli[0] = su2_tmp[0];      
    r_pauli[1] = su2_tmp[1];
    r_pauli[2] = su2_tmp[2];
    r_pauli[3] = su2_tmp[3];

      
    //go back to a su2 matrix in standard basis
    w[0].re = r_pauli[0];
    w[0].im = r_pauli[3];
    w[1].re = r_pauli[2];
    w[1].im = r_pauli[1];
    w[2].re = -r_pauli[2];
    w[2].im = r_pauli[1];
    w[3].re = r_pauli[0];
    w[3].im = -r_pauli[3];
 
    //extend to SU3
    hmc_su3matrix extW, tmp;
    extend (&extW, order[i], w);
	
    multiply_su3matrices(&tmp, &extW, &su3_in);
    copy_su3matrix(&su3_in, &tmp);
  }
  copy_su3matrix(out, &su3_in);
  *cter_out = cter;
  
  return;
}

void testing_heatbath(hmc_su3matrix * in, hmc_staplematrix * staple_in, hmc_su3matrix * out, hmc_float beta, hmc_float * rnd_array, int * cter_out){
  hmc_staplematrix W;
  hmc_su3matrix su3_in;
  copy_su3matrix(&su3_in, in);
            
  int cter = 0;
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
  order[0] = (int)( rnd_array[cter] * 3 ) + 1;
  cter ++;
  do
    {order[1] = (int)( rnd_array[cter] * 3 ) + 1;
     cter ++;}
  while (order[1] == order[0]);
  order[2] = 6 - order[1] - order[0];
  
  hmc_complex det = det_su3matrix(&su3_in);
      hmc_complex detadj = complexconj(&det);
      hmc_complex detsqnorm = complexmult(&det, &detadj);
      if( (detsqnorm.re - hmc_one_f) <= projectioneps)
	project_su3(&su3_in);

  for(int i=0; i<3; i++)
  {
    //Produkt aus U und staple, saved in W
    multiply_staplematrix(&W, &su3_in, staple_in);
    //Reduktion auf SU(2) submatrix
    reduction (w, W, order[i]);
    //Darstellung von W in Paulibasis
    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);

    //Berechnung der Normierung k, entspricht Determinante in normaler Basis
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
    w_pauli[0] = hmc_float(1)/k * w_pauli[0];
    w_pauli[1] = hmc_float(-1)/k * w_pauli[1];
    w_pauli[2] = hmc_float(-1)/k * w_pauli[2];
    w_pauli[3] = hmc_float(-1)/k * w_pauli[3];
	
    //beta' = 2Beta/Nc*k
    hmc_float beta_neu =  2.*beta / hmc_float(NC)*k;
    
    
    //Neuer Link in Paulibasis mittels Kennedy-Pendleton-Algorithmus
    //function  SU2Update(r_pauli, beta_neu):
    hmc_float delta;
    hmc_float a0 ;
    hmc_float eta ;
    hmc_float alpha = beta_neu;

    do
    {
      hmc_float rnd1 = rnd_array[cter];
      hmc_float rnd2 = rnd_array[cter+1];
      hmc_float rnd3 = rnd_array[cter+2];
      hmc_float rnd4 = rnd_array[cter+3];
      delta = -log(rnd1)/alpha*pow(cos(2. * PI * rnd2), 2.) -log(rnd3)/alpha;
      a0 = 1.-delta;
      eta = rnd4;
      cter += 4;
    }while ( (1.-0.5*delta) < eta*eta);
    hmc_float rnd5 = rnd_array[cter];
    hmc_float rnd6 = rnd_array[cter+1];
    cter += 2;
    hmc_float phi = 2.*PI*rnd5;
    hmc_float theta = asin(2.*rnd6 - 1.);
    r_pauli[0] = a0;
    r_pauli[1] = sqrt(1.-a0 * a0)*cos(theta) * cos(phi);
    r_pauli[2] = sqrt(1.-a0 * a0)*cos(theta) * sin(phi);
    r_pauli[3] = sqrt(1.-a0 * a0)*sin(theta);
    
    
    //Multipliziere neuen Link mit w, alles in Paulibasis
    hmc_float su2_tmp[su2_entries];
    su2_tmp[0] = r_pauli[0]*w_pauli[0] - r_pauli[1]*w_pauli[1] - r_pauli[2]*w_pauli[2] - r_pauli[3]*w_pauli[3] ;
    su2_tmp[1] = w_pauli[0]*r_pauli[1] + w_pauli[1]*r_pauli[0] - r_pauli[2]*w_pauli[3] + r_pauli[3]*w_pauli[2] ;
    su2_tmp[2] = w_pauli[0]*r_pauli[2] + w_pauli[2]*r_pauli[0] - r_pauli[3]*w_pauli[1] + r_pauli[1]*w_pauli[3] ;
    su2_tmp[3] = w_pauli[0]*r_pauli[3] + w_pauli[3]*r_pauli[0] - r_pauli[1]*w_pauli[2] + r_pauli[2]*w_pauli[1] ;
    r_pauli[0] = su2_tmp[0];      
    r_pauli[1] = su2_tmp[1];
    r_pauli[2] = su2_tmp[2];
    r_pauli[3] = su2_tmp[3];

      
    //go back to a su2 matrix in standard basis
    w[0].re = r_pauli[0];
    w[0].im = r_pauli[3];
    w[1].re = r_pauli[2];
    w[1].im = r_pauli[1];
    w[2].re = -r_pauli[2];
    w[2].im = r_pauli[1];
    w[3].re = r_pauli[0];
    w[3].im = -r_pauli[3];
 
    //extend to SU3
    hmc_su3matrix extW, tmp;
    extend (&extW, order[i], w);
	
    multiply_su3matrices(&tmp, &extW, &su3_in);
    copy_su3matrix(&su3_in, &tmp);
  }
  copy_su3matrix(out, &su3_in);
  *cter_out = cter;
  
  return;
}





