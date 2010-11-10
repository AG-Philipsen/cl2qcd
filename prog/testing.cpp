#include "testing.h"

void testing_spinor() {
  hmc_full_spinor_field test;
  printf("global squarenorm: %f\n",global_squarenorm(&test));
  set_zero_spinor(&test);

  printf("global squarenorm: %f\n",global_squarenorm(&test));

  fill_with_one(&test, 0, 0, 0);

  printf("global squarenorm: %f\n",global_squarenorm(&test));

  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int j=0; j<NSPIN*NC; j++) {
	fill_with_one(&test, n, t, j);
      }
    }
  }

  printf("global squarenorm: %f\n",global_squarenorm(&test));
  return;
}

void print_su3mat(hmc_su3matrix* A){
#ifdef _RECONSTRUCT_TWELVE_
  printf("| (%f,%f)\t(%f,%f)\t(%f,%f) |\n",(*A)[0].re,(*A)[0].im,(*A)[2].re,(*A)[2].im,(*A)[4].re,(*A)[4].im);
  printf("| (%f,%f)\t(%f,%f)\t(%f,%f) |\n",(*A)[1].re,(*A)[1].im,(*A)[3].re,(*A)[3].im,(*A)[5].re,(*A)[5].im);
  hmc_complex ca = reconstruct_su3(A,0);
  hmc_complex cb = reconstruct_su3(A,1);
  hmc_complex cc = reconstruct_su3(A,2);
  printf("| (%f,%f)\t(%f,%f)\t(%f,%f) |\n",ca.re,ca.im,cb.re,cb.im,cc.re,cc.im);
#else
  for(int a = 0; a<NC; a++) 
  printf("| (%f,%f)\t(%f,%f)\t(%f,%f) |\n",(*A)[a][0].re,(*A)[a][0].im,(*A)[a][1].re,(*A)[a][1].im,(*A)[a][2].re,(*A)[a][2].im);
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
