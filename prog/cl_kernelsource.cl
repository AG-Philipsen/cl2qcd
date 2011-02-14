#ifdef _USEDOUBLEPREC_
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
//#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

#include "globaldefs.h" //NDIM, NSPIN, NC
#include "types.h"

//for hmc_ocl_su3matrix
#ifdef _RECONSTRUCT_TWELVE_
#define SU3SIZE NC*(NC-1)
#define STAPLEMATRIXSIZE NC*NC
#else
#define SU3SIZE NC*NC
#define STAPLEMATRIXSIZE NC*NC
#endif


#ifdef _USEDOUBLEPREC_
hmc_float const projectioneps = 10.e-12;
#else
hmc_float const projectioneps = 10.e-6;
#endif


//site = pos + VOLSPACE*t =  x + y*NSPACE + z*NSPACE*NSPACE + VOLSPACE*t
//idx = mu + NDIM*site
int inline ocl_gaugefield_element(int c, int a, int b, int mu, int spacepos, int t){
#ifdef _RECONSTRUCT_TWELVE_
  return c + 2*a + 2*(NC-1)*b+2*NC*(NC-1)*mu+2*NC*(NC-1)*NDIM*spacepos+2*NC*(NC-1)*NDIM*VOLSPACE*t;
#else
  return c + 2*a + 2*NC*b+2*NC*NC*mu+2*NC*NC*NDIM*spacepos+2*NC*NC*NDIM*VOLSPACE*t;
#endif
}

int inline ocl_su3matrix_element(int a, int b){
#ifdef _RECONSTRUCT_TWELVE_
  return a + (NC-1)*b;
#else
  return a + NC*b;
#endif
}

//it is assumed that idx iterates only over half the number of sites
void inline get_even_site(int idx, int * out_space, int * out_t){
  int x,y,z,t;
  x = idx;
  t = convert_int(idx/(VOLSPACE/2));
  x -= t*VOLSPACE/2;
  z = convert_int(x/(NSPACE*NSPACE/2));
  x -= z*NSPACE*NSPACE/2;
  y = convert_int(x/NSPACE);
  x -= y*NSPACE;
  (*out_space) =  convert_int((z+t)%2)*(1 + 2*x - convert_int (2*x/NSPACE)) + convert_int((t+z+1)%2)*(2*x + convert_int (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
  (*out_t) = t;
}

//it is assumed that idx iterates only over half the number of sites
void inline get_odd_site(int idx, int * out_space, int * out_t){
  int x,y,z,t;
  x = idx;
  t = convert_int(idx/(VOLSPACE/2));
  x -= t*VOLSPACE/2;
  z = convert_int(x/(NSPACE*NSPACE/2));
  x -= z*NSPACE*NSPACE/2;
  y = convert_int(x/NSPACE);
  x -= y*NSPACE;
  (*out_space) =  convert_int((z+t+1)%2)*(1 + 2*x - convert_int (2*x/NSPACE)) + convert_int((t+z)%2)*(2*x + convert_int (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
  (*out_t) = t;
}

int inline get_global_pos(int spacepos, int t){
  return spacepos + VOLSPACE * t;
}

hmc_complex complexconj(__private hmc_complex *in){
  hmc_complex tmp;
  tmp.re = (*in).re;
  tmp.im = -(*in).im;
  return tmp;
}

hmc_complex complexmult(__private hmc_complex *a,__private hmc_complex *b){
  hmc_complex res;
  res.re = (*a).re*(*b).re - (*a).im*(*b).im;
  res.im = (*a).im*(*b).re + (*a).re*(*b).im;
  return res;
}

hmc_complex complexadd(__private hmc_complex *a,__private hmc_complex *b){
  hmc_complex res;
  res.re = (*a).re + (*b).re;
  res.im = (*a).im + (*b).im;
  return res;
}

hmc_complex complexsubtract(__private hmc_complex *a,__private hmc_complex *b){
  hmc_complex res;
  res.re = (*a).re - (*b).re;
  res.im = (*a).im - (*b).im;
  return res;
}

void complexaccumulate(__private hmc_complex *inout,__private hmc_complex *incr){
  (*inout).re += (*incr).re;
  (*inout).im += (*incr).im;
  return;
}

void get_su3matrix(__private hmc_ocl_su3matrix* out, __global hmc_ocl_gaugefield * in, int spacepos, int timepos, int mu){
#ifdef _RECONSTRUCT_TWELVE_
  for(int a=0; a<NC-1; a++) {
#else
    for(int a=0; a<NC; a++) {
#endif
      for(int b=0; b<NC; b++) { 
	out[ocl_su3matrix_element(a,b)].re = in[ocl_gaugefield_element(0,a,b,mu,spacepos,timepos)];
	out[ocl_su3matrix_element(a,b)].im = in[ocl_gaugefield_element(1,a,b,mu,spacepos,timepos)];
      }
#ifdef _RECONSTRUCT_TWELVE_
    }
#else
  }
#endif
  return;
}

void put_su3matrix(__global hmc_ocl_gaugefield * field, __private hmc_ocl_su3matrix * in, int spacepos, int timepos, int mu){
#ifdef _RECONSTRUCT_TWELVE_
  for(int a=0; a<NC-1; a++) {
#else
    for(int a=0; a<NC; a++) {
#endif
      for(int b=0; b<NC; b++) { 
	field[ocl_gaugefield_element(0,a,b,mu,spacepos,timepos)] = in[ocl_su3matrix_element(a,b)].re;
	field[ocl_gaugefield_element(1,a,b,mu,spacepos,timepos)] = in[ocl_su3matrix_element(a,b)].im;
      }
#ifdef _RECONSTRUCT_TWELVE_
    }
#else
  }
#endif
  return;
}

void copy_su3matrix(__private hmc_ocl_su3matrix * out, __private hmc_ocl_su3matrix * in){
  for(int n = 0; n<SU3SIZE; n++) {
    (out[n]).re = (in[n]).re;
    (out[n]).im = (in[n]).im;
  }
  return;
}

void unit_su3matrix(__private hmc_ocl_su3matrix *out){
#ifdef _RECONSTRUCT_TWELVE_
  for(int a = 0; a < NC-1; a++) {
#else 
    for(int a = 0; a < NC; a++) {
#endif
      for(int b = 0; b < NC; b++) {
	if(a==b) {
	  out[ocl_su3matrix_element(a,b)].re = hmc_one_f;
	} else {
	  out[ocl_su3matrix_element(a,b)].re = 0;
	}
	out[ocl_su3matrix_element(a,b)].im = 0;
      }
#ifdef _RECONSTRUCT_TWELVE_
    }
#else 
  }
#endif
  return;
}

void zero_su3matrix(__private hmc_ocl_su3matrix *out){
  for(int n=0; n<SU3SIZE; n++) {
    out[n].re = 0;
    out[n].im = 0;
  }
  return;
}

#ifdef _RECONSTRUCT_TWELVE_
hmc_complex reconstruct_su3(__private hmc_ocl_su3matrix *in, int ncomp){
  int jplusone = (ncomp+1)%NC;
  int jplustwo = (ncomp+2)%NC;
  hmc_complex first = complexmult(&(in[(NC-1)*jplusone]),&(in[1+(NC-1)*jplustwo]));
  hmc_complex second = complexmult(&(in[(NC-1)*jplustwo]),&(in[1+(NC-1)*jplusone]));
  hmc_complex result = complexsubtract(&first,&second);
  return complexconj(&result);
}
#endif

void multiply_su3matrices(__private hmc_ocl_su3matrix *out,__private hmc_ocl_su3matrix *p,__private hmc_ocl_su3matrix *q){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
      out[n].re=0;
      out[n].im=0;
      for(int j=0;j<NC;j++) {
	int k = (int)(n/(NC-1));
	int i = n - (NC-1)*k;
	int np = i + (NC-1)*j;
	hmc_complex qcomponent;
	if(j==2) {
	  qcomponent = reconstruct_su3(q,k);
	} else {
	  int nq = j + (NC-1)*k;
	  qcomponent = q[nq];
	}
	hmc_complex tmp = complexmult(&p[np],&qcomponent);
	complexaccumulate(&out[n],&tmp);
      }
    }
#else
    out[0].re = p[0].re*q[0].re + p[3].re*q[1].re + p[6].re*q[2].re -
	       (p[0].im*q[0].im + p[3].im*q[1].im + p[6].im*q[2].im );
    out[3].re = p[0].re*q[3].re + p[3].re*q[4].re + p[6].re*q[5].re-
	       (p[0].im*q[3].im + p[3].im*q[4].im + p[6].im*q[5].im);
    out[6].re = p[0].re*q[6].re + p[3].re*q[7].re + p[6].re*q[8].re -
	       (p[0].im*q[6].im + p[3].im*q[7].im + p[6].im*q[8].im);
    out[1].re = p[1].re*q[0].re + p[4].re*q[1].re + p[7].re*q[2].re -
	       (p[1].im*q[0].im + p[4].im*q[1].im + p[7].im*q[2].im);
    out[4].re = p[1].re*q[3].re + p[4].re*q[4].re + p[7].re*q[5].re -
	       (p[1].im*q[3].im + p[4].im*q[4].im + p[7].im*q[5].im);
    out[7].re = p[1].re*q[6].re + p[4].re*q[7].re + p[7].re*q[8].re -
	       (p[1].im*q[6].im + p[4].im*q[7].im + p[7].im*q[8].im);
    out[2].re = p[2].re*q[0].re + p[5].re*q[1].re + p[8].re*q[2].re -
	       (p[2].im*q[0].im + p[5].im*q[1].im + p[8].im*q[2].im);
    out[5].re = p[2].re*q[3].re + p[5].re*q[4].re + p[8].re*q[5].re -
	       (p[2].im*q[3].im + p[5].im*q[4].im + p[8].im*q[5].im);
    out[8].re = p[2].re*q[6].re + p[5].re*q[7].re + p[8].re*q[8].re -
	       (p[2].im*q[6].im + p[5].im*q[7].im + p[8].im*q[8].im); 
  
    out[0].im = p[0].re*q[0].im + p[3].re*q[1].im + p[6].re*q[2].im + 
		p[0].im*q[0].re + p[3].im*q[1].re + p[6].im*q[2].re;
    out[3].im = p[0].re*q[3].im + p[3].re*q[4].im + p[6].re*q[5].im + 
		p[0].im*q[3].re + p[3].im*q[4].re + p[6].im*q[5].re;
    out[6].im = p[0].re*q[6].im + p[3].re*q[7].im + p[6].re*q[8].im + 
		p[0].im*q[6].re + p[3].im*q[7].re + p[6].im*q[8].re;
    out[1].im = p[1].re*q[0].im + p[4].re*q[1].im + p[7].re*q[2].im + 
		p[1].im*q[0].re + p[4].im*q[1].re + p[7].im*q[2].re;
    out[4].im = p[1].re*q[3].im + p[4].re*q[4].im + p[7].re*q[5].im + 
		p[1].im*q[3].re + p[4].im*q[4].re + p[7].im*q[5].re ;
    out[7].im = p[1].re*q[6].im + p[4].re*q[7].im + p[7].re*q[8].im + 
		p[1].im*q[6].re + p[4].im*q[7].re + p[7].im*q[8].re;
    out[2].im = p[2].re*q[0].im + p[5].re*q[1].im + p[8].re*q[2].im + 
		p[2].im*q[0].re + p[5].im*q[1].re + p[8].im*q[2].re;
    out[5].im = p[2].re*q[3].im + p[5].re*q[4].im + p[8].re*q[5].im + 
		p[2].im*q[3].re + p[5].im*q[4].re + p[8].im*q[5].re;
    out[8].im = p[2].re*q[6].im + p[5].re*q[7].im + p[8].re*q[8].im + 
		p[2].im*q[6].re + p[5].im*q[7].re + p[8].im*q[8].re;
		

//     for(int i=0; i<NC; i++) {
//       for(int k=0; k<NC; k++) {
// 	out[ocl_su3matrix_element(i,k)].re=0;
// 	out[ocl_su3matrix_element(i,k)].im=0;
// 	for(int j=0;j<NC;j++) {
// 	  hmc_complex tmp = complexmult(&p[ocl_su3matrix_element(i,j)],&q[ocl_su3matrix_element(j,k)]);
// 	  if(i==k) complexaccumulate(&out[ocl_su3matrix_element(i,k)],&tmp);
// 	}
//       }
//     }
#endif
  return;
}

void accumulate_su3matrix_prod(__private hmc_ocl_su3matrix *acc,__private hmc_ocl_su3matrix *mult){
  hmc_ocl_su3matrix tmp[SU3SIZE];
  multiply_su3matrices(tmp,acc,mult);
  copy_su3matrix(acc,tmp);
  return;
}


void adjoin_su3matrix(__private hmc_ocl_su3matrix * mat){
#ifdef _RECONSTRUCT_TWELVE_
  hmc_ocl_su3matrix tmp[SU3SIZE];
  copy_su3matrix(tmp, mat);
  for(int n=0; n<NC*(NC-1); n++) {
    int j = (int)(n/(NC-1));
    int i = n - (NC-1)*j;
    hmc_complex element;
    if ( j==2 ) {
      element = reconstruct_su3(tmp,i);
    } else {
      int nnew = j + (NC-1)*i;
      element = tmp[nnew];
    }
    mat[n] = complexconj(&element);
  }
#else
  hmc_ocl_su3matrix tmp[SU3SIZE];
  copy_su3matrix(tmp, mat);
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      mat[ocl_su3matrix_element(a,b)] =
      complexconj(&(tmp[ocl_su3matrix_element(b,a)]));
    }
  }
#endif
  return;
}

hmc_complex trace_su3matrix(__private hmc_ocl_su3matrix * mat){
  hmc_complex trace;
#ifdef _RECONSTRUCT_TWELVE_
  trace = reconstruct_su3(mat,NC-1);
  for(int n=0; n<(NC-1); n++) complexaccumulate(&trace,&(mat[NC*n]));
#else
  trace.re=0;
  trace.im=0;
  for(int a=0; a<NC; a++) complexaccumulate(&trace,&(mat[ocl_su3matrix_element(a,a)]));;
#endif
  return trace;
}

hmc_complex det_su3matrix(__private hmc_ocl_su3matrix * U){
#ifdef _RECONSTRUCT_TWELVE_
  hmc_complex det;
  det.re=0;
  det.im=0;
  hmc_complex subdet;
  subdet.re=0;
  subdet.im=0;
  hmc_complex tmp1;
  hmc_complex tmp2;
  hmc_complex tmp3;
  tmp1 = complexmult( &U[0], &U[3] );
  tmp2 = reconstruct_su3(U,2);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&det,&tmp3);

  tmp1 = complexmult( &U[2], &U[5] );
  tmp2 = reconstruct_su3(U,0);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&det,&tmp3);

  tmp1 = complexmult( &U[4], &U[1] );
  tmp2 = reconstruct_su3(U,1);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&det,&tmp3);

  tmp1 = complexmult( &U[3], &U[4] );
  tmp2 = reconstruct_su3(U,0);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&subdet,&tmp3);

  tmp1 = complexmult( &U[5], &U[0] );
  tmp2 = reconstruct_su3(U,1);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&subdet,&tmp3);

  tmp1 = complexmult( &U[1], &U[2] );
  tmp2 = reconstruct_su3(U,2);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&subdet,&tmp3);

  det.re -= subdet.re;
  det.im -= subdet.im;

#else
  hmc_complex det, det1, det2, det3, det4, det5, det6, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
  det.re=0;
  det.im=0;
  tmp1 = complexmult( &U[ocl_su3matrix_element(1,1)], &U[ocl_su3matrix_element(2,2)] );
  det1 = complexmult( &U[ocl_su3matrix_element(0,0)] , &tmp1);
  tmp2 = complexmult( &U[ocl_su3matrix_element(1,2)], &U[ocl_su3matrix_element(2,0)] );
  det2 = complexmult( &U[ocl_su3matrix_element(0,1)] , &tmp2);
  tmp3 = complexmult( &U[ocl_su3matrix_element(1,0)], &U[ocl_su3matrix_element(2,1)] );
  det3 = complexmult( &U[ocl_su3matrix_element(0,2)] , &tmp3);
  tmp4 = complexmult( &U[ocl_su3matrix_element(1,1)], &U[ocl_su3matrix_element(2,0)] );
  det4 = complexmult( &U[ocl_su3matrix_element(0,2)] , &tmp4);
  tmp5 = complexmult( &U[ocl_su3matrix_element(1,0)], &U[ocl_su3matrix_element(2,2)] );
  det5 = complexmult( &U[ocl_su3matrix_element(0,1)] , &tmp5);
  tmp6 = complexmult( &U[ocl_su3matrix_element(1,2)], &U[ocl_su3matrix_element(2,1)] );
  det6 = complexmult( &U[ocl_su3matrix_element(0,0)] , &tmp6);

  det.re = det1.re + det2.re + det3.re - det4.re - det5.re - det6.re;
  det.im = det1.im + det2.im + det3.im - det4.im - det5.im - det6.im;

#endif
  return det;
}

void project_su3(__private hmc_ocl_su3matrix *U){
  hmc_complex det = det_su3matrix(U);
  hmc_float detsqunorm = det.re*det.re + det.im*det.im;

  hmc_float phi;
  if(det.re*det.re<projectioneps) { 
  phi = PI/2.;
  } else {
  phi = atan(det.im/det.re);
  if(det.re<0) phi += PI;
  }

  hmc_complex norm;
  norm.re = pow(detsqunorm,hmc_one_f/6.)*cos(phi/3.);
  norm.im = pow(detsqunorm,hmc_one_f/6.)*sin(phi/3.);
  
  hmc_float normsqunorm = norm.re*norm.re+norm.im*norm.im;

  #ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) { 
    hmc_complex tmp = (U)[n];
    (U)[n].re = (tmp.re*norm.re+tmp.im*norm.im)/normsqunorm;
    (U)[n].im = (tmp.im*norm.re-tmp.re*norm.im)/normsqunorm;
  }
  #else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      hmc_complex tmp = (U)[ocl_su3matrix_element(a,b)]; 
      (U)[ocl_su3matrix_element(a,b)].re = (tmp.re*norm.re+tmp.im*norm.im)/normsqunorm;
      (U)[ocl_su3matrix_element(a,b)].im = (tmp.im*norm.re-tmp.re*norm.im)/normsqunorm;
    }
  }
  #endif
  return;
}

void adjoin_su3(__global hmc_ocl_gaugefield * in,__global hmc_ocl_gaugefield * out){
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int mu=0; mu<NDIM; mu++) {
	hmc_ocl_su3matrix tmp[SU3SIZE];
	get_su3matrix(tmp, in, n, t, mu);
	adjoin_su3matrix(tmp);
	put_su3matrix(out, tmp, n, t, mu);
      }
    }
  }
  return;
}

hmc_complex global_trace_su3(__global hmc_ocl_gaugefield * field, int mu) {
  hmc_complex sum;
  sum.re = 0;
  sum.im = 0;
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      hmc_ocl_su3matrix tmp[SU3SIZE];
      get_su3matrix(tmp, field, n, t, mu);
      sum.re += trace_su3matrix(tmp).re;
      sum.im += trace_su3matrix(tmp).im;
    }
  }
  return sum;
}

//LZ: everything above has been tested...

void copy_staplematrix(__private hmc_ocl_staplematrix *out,__private hmc_ocl_staplematrix *in){
  for(int n=0; n<NC*NC; n++) {
    out[n] = in[n];
  }
  return;
}

void zero_staplematrix(__private hmc_ocl_staplematrix * u){
  for(int n=0; n<STAPLEMATRIXSIZE; n++) {
    u[n].re = 0;
    u[n].im = 0;
  }
  return;
}

void unit_staplematrix(__private hmc_ocl_staplematrix * u){
  u[0].re = 1.;
  u[0].im = 0;
  u[4].re = 1.;
  u[4].im = 0;
  u[8].re = 1.;
  u[8].im = 0;
  
  return;
}

void multiply_staplematrix(__private hmc_ocl_staplematrix *out, __private hmc_ocl_su3matrix *p,__private  hmc_ocl_staplematrix *q){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
      out[n].re=0;
      out[n].im=0;
      for(int j=0;j<NC;j++) {
	int k = (int)(n/(NC-1));
	int i = n - (NC-1)*k;
	int np = i + (NC-1)*j;
	hmc_complex qcomponent;
	if(j==2) {
// 	  qcomponent = reconstruct_su3(q,k);
          qcomponent = q[NC*(NC-1)+k];
	} else {
	  int nq = j + (NC-1)*k;
	  qcomponent = q[nq];
	}
	hmc_complex tmp = complexmult(&p[np],&qcomponent);
	complexaccumulate(&out[n],&tmp);
      }
    }
    //the left components:
    hmc_complex X = reconstruct_su3(p,0);
    hmc_complex Y = reconstruct_su3(p,1);
    hmc_complex Z = reconstruct_su3(p,2);
    hmc_complex tmp;
    out[6].re=0;
    out[6].im=0;
    out[7].re=0;
    out[7].im=0;
    out[8].re=0;
    out[8].im=0;
    
    tmp = complexmult(&X,&q[0]);
    complexaccumulate(&out[6],&tmp);
    tmp = complexmult(&Y,&q[1]);
    complexaccumulate(&out[6],&tmp);
    tmp = complexmult(&Z,&q[6]);
    complexaccumulate(&out[6],&tmp);

    tmp = complexmult(&X,&q[2]);
    complexaccumulate(&out[7],&tmp);
    tmp = complexmult(&Y,&q[3]);
    complexaccumulate(&out[7],&tmp);
    tmp = complexmult(&Z,&q[7]);
    complexaccumulate(&out[7],&tmp);

    tmp = complexmult(&X,&q[4]);
    complexaccumulate(&out[8],&tmp);
    tmp = complexmult(&Y,&q[5]);
    complexaccumulate(&out[8],&tmp);
    tmp = complexmult(&Z,&q[8]);
    complexaccumulate(&out[8],&tmp);
    
#else
    multiply_su3matrices(out, p, q);
    /*
    for(int i=0; i<NC; i++) {
      for(int k=0; k<NC; k++) {
	out[ocl_su3matrix_element(i,k)].re=0;
	out[ocl_su3matrix_element(i,k)].im=0;
	for(int j=0;j<NC;j++) {
	  hmc_complex tmp = complexmult(&p[ocl_su3matrix_element(i,j)],&q[ocl_su3matrix_element(j,k)]);
	  complexaccumulate(&out[ocl_su3matrix_element(i,k)],&tmp);
	}
      }
    }
    */
#endif
  return;
}

void accumulate_su3matrices_add(__private hmc_ocl_staplematrix *p,__private hmc_ocl_su3matrix *q){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
    complexaccumulate(&(p[n]), &(q[n]));
  }
  for(int n=NC*(NC-1);  n<NC*NC; n++) {
    hmc_complex tmp = reconstruct_su3(q, n-NC*(NC-1)); 
    complexaccumulate(&(p[n]), &(tmp));
  }  
#else
  for(int k=0; k<NC*NC; k++) {
    complexaccumulate(&(p[k]),&(q[k]));
  }
#endif
  return;
}

void reduction (hmc_complex dest[su2_entries],__private hmc_ocl_staplematrix* src, const int rand){
#ifdef _RECONSTRUCT_TWELVE_
  if(rand == 1)
  { 
    dest[0] = src[0]; 
    dest[1] = src[2]; 
    dest[2] = src[1]; 
    dest[3] = src[3]; 
  } 
  else if (rand==2) 
  { 
    dest[0] = src[3]; 
    dest[1] = src[5]; 
    dest[2] = src[7];
    dest[3] = src[8];
  } 
  else if (rand==3)
  { 
    dest[0] = src[0]; 
    dest[1] = src[4]; 
    dest[2] = src[6];
    dest[3] = src[8];
  } 
#else
  if(rand == 1)
  {
    dest[0] = src[ocl_su3matrix_element(0,0)];
    dest[1] = src[ocl_su3matrix_element(0,1)];
    dest[2] = src[ocl_su3matrix_element(1,0)];
    dest[3] = src[ocl_su3matrix_element(1,1)];
  }
  else if (rand==2)
  {
    dest[0] = src[ocl_su3matrix_element(1,1)];
    dest[1] = src[ocl_su3matrix_element(1,2)];
    dest[2] = src[ocl_su3matrix_element(2,1)];
    dest[3] = src[ocl_su3matrix_element(2,2)];
  }
  else if (rand==3)
  {
    dest[0] = src[ocl_su3matrix_element(0,0)];
    dest[1] = src[ocl_su3matrix_element(0,2)];
    dest[2] = src[ocl_su3matrix_element(2,0)];
    dest[3] = src[ocl_su3matrix_element(2,2)];
  }
#endif
}

// return an SU2 matrix (std basis) extended to SU3 (std basis)
void extend (__private hmc_ocl_su3matrix * dest, const int random, hmc_complex src[su2_entries]){
#ifdef _RECONSTRUCT_TWELVE_
  if (random == 1){
    dest[0] = src[0];
    dest[2] = src[1];
    dest[4] = hmc_complex_zero;
    dest[1] = src[2];
    dest[3] = src[3];
    dest[5] = hmc_complex_zero;
  }
  else if (random == 2){
    dest[0] = hmc_complex_one;
    dest[2] = hmc_complex_zero;
    dest[4] = hmc_complex_zero;
    dest[1] = hmc_complex_zero;
    dest[3] = src[0];
    dest[5] = src[1];
  }
  else if (random == 3){
    dest[0] = src[0];
    dest[2] = hmc_complex_zero;
    dest[4] = src[1];
    dest[1] = hmc_complex_zero;
    dest[3] = hmc_complex_one;
    dest[5] = hmc_complex_zero;
  }

#else
  if (random == 1){
    dest[ocl_su3matrix_element(0,0)] = src[0];
    dest[ocl_su3matrix_element(0,1)] = src[1];
    dest[ocl_su3matrix_element(0,2)] = hmc_complex_zero;
    dest[ocl_su3matrix_element(1,0)] = src[2];
    dest[ocl_su3matrix_element(1,1)] = src[3];
    dest[ocl_su3matrix_element(1,2)] = hmc_complex_zero;
    dest[ocl_su3matrix_element(2,0)] = hmc_complex_zero;
    dest[ocl_su3matrix_element(2,1)] = hmc_complex_zero;
    dest[ocl_su3matrix_element(2,2)] = hmc_complex_one;
  }
  else if (random == 2){
    dest[ocl_su3matrix_element(0,0)] = hmc_complex_one;
    dest[ocl_su3matrix_element(0,1)] = hmc_complex_zero;
    dest[ocl_su3matrix_element(0,2)] = hmc_complex_zero;
    dest[ocl_su3matrix_element(1,0)] = hmc_complex_zero;
    dest[ocl_su3matrix_element(1,1)] = src[0];
    dest[ocl_su3matrix_element(1,2)] = src[1];
    dest[ocl_su3matrix_element(2,0)] = hmc_complex_zero;
    dest[ocl_su3matrix_element(2,1)] = src[2];
    dest[ocl_su3matrix_element(2,2)] = src[3];
  }
  else if (random == 3){
    dest[ocl_su3matrix_element(0,0)] = src[0];
    dest[ocl_su3matrix_element(0,1)] = hmc_complex_zero;
    dest[ocl_su3matrix_element(0,2)] = src[1];
    dest[ocl_su3matrix_element(1,0)] = hmc_complex_zero;
    dest[ocl_su3matrix_element(1,1)] = hmc_complex_one;
    dest[ocl_su3matrix_element(1,2)] = hmc_complex_zero;
    dest[ocl_su3matrix_element(2,0)] = src[2];
    dest[ocl_su3matrix_element(2,1)] = hmc_complex_zero;
    dest[ocl_su3matrix_element(2,2)] = src[3];
  }
#endif
}

int get_nspace(const int* coord){
  int n = coord[1] +  NSPACE*coord[2] + NSPACE*NSPACE*coord[3];
  return n;
}

int get_spacecoord(const int nspace, const int dir){
  int res = convert_int(nspace/(NSPACE*NSPACE));
  if(dir==3) return res;
  int acc = res;
  res = convert_int(nspace/(NSPACE)) - NSPACE*acc;
  if(dir==2) return res;
  acc = NSPACE*acc + res;
  res = nspace - NSPACE*acc;
  return res;
}

void get_allspacecoord(const int nspace, int coord[NDIM]){
  int res = convert_int(nspace/(NSPACE*NSPACE));
  coord[3] = res;
  int acc = res;
  res = convert_int(nspace/(NSPACE)) - NSPACE*acc;
  coord[2] = res;
  acc = NSPACE*acc + res;
  res = nspace - NSPACE*acc;
  coord[1] = res;
}

int get_neighbor(const int nspace,const int dir) {
  int coord[NDIM];
//   coord[0]=0;
//   for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(nspace,j);
  get_allspacecoord(nspace, coord);
  coord[dir] = (coord[dir] + 1)%NSPACE;
  return get_nspace(coord);
}

int get_lower_neighbor(const int nspace, int const dir) {
  int coord[NDIM];
//   coord[0]=0;
//   for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(nspace,j);
  get_allspacecoord(nspace, coord);
  coord[dir] = (coord[dir] - 1 + NSPACE)%NSPACE;
  return get_nspace(coord);
}
//CP: old generator
/*
typedef uint4 hmc_ocl_ran;

inline unsigned _ocl_taus_step( unsigned *z, int S1, int S2, int S3, unsigned M)
{
	unsigned b = ( ( ( (*z) << S1 )^(*z) ) >> S2 );
	return *z = ( ( ( (*z) & M ) << S3 )^b );
}

inline unsigned _ocl_LCG_step( unsigned *z, unsigned A, unsigned C)
{
	return *z= ( A * (*z) + C );
}

inline float ocl_new_ran( __global hmc_ocl_ran* state )
{
	uint x,y,z,w;
	x = (*state).x;
	y = (*state).y;
	z = (*state).z;
	w = (*state).w;
	float res = 2.328306436538696e-10f*( _ocl_taus_step( &x, 13, 19, 12, 4294967294ul)^
	                                     _ocl_taus_step( &y, 2, 25, 4,   4294967288ul)^
	                                     _ocl_taus_step( &z, 3, 11, 17,  4294967280ul)^
	                                     _ocl_LCG_step(  &w, 1664525,    1013904223ul) );
	*state = (hmc_ocl_ran)(x,y,z,w);
	return res;
}
*/

//CP: new one

/**
 * RNG state for the NR3
 */
// typedef ulong4 nr3_state;
typedef ulong4 hmc_ocl_ran;


/**
 * Calculate the next random number as described in NR3
 */
// inline ulong nr3_int64( nr3_state * state ) {
inline ulong nr3_int64(__global hmc_ocl_ran * state ) {
	(*state).x = (*state).x * 2862933555777941757L + 7046029254386353087L;
	(*state).y ^= (*state).y >> 17; (*state).y ^= (*state).y << 31; (*state).y ^= (*state).y >> 8;
	(*state).z = 4294957665U*((*state).z & 0xffffffff) + ((*state).z >> 32);
	ulong tmp = (*state).x ^ ((*state).x << 21); tmp ^= tmp >> 35; tmp ^= tmp << 4;
	return (tmp + (*state).y) ^ (*state).z;
}

/**
 * Calculate the next random number and return as float
 * FIXME this conversion is probably broken
 */
// inline float nr3_float( nr3_state * state )
inline float ocl_new_ran(__global hmc_ocl_ran * state )
{
	return 5.42101086242752217E-20f * nr3_int64( state );
}

/**
 * Calculate the next random number and return as int32
 */
// inline uint nr3_int32( nr3_state * state )
inline uint nr3_int32(__global hmc_ocl_ran * state )
{
	return (uint) nr3_int64( state );
}




int random_int( int range, __global hmc_ocl_ran* taus_state )
{
	return convert_int( ocl_new_ran( taus_state ) * range );
}

//returns 1,2,3 in a random way
void random_1_2_3 (int rand[3], __global hmc_ocl_ran * rnd) { 

	// first value can be any value 1..3
	rand[ 0 ] = random_int( 3, rnd ) + 1;

	// second value must be 1..3 but not equal to the first value
	do
	{
	  rand[ 1 ] = random_int( 3, rnd ) + 1;
	} while( rand[ 0 ] == rand[ 1 ] );

	// third value must take the remaining value from 1..3
	rand[ 2 ] = 6 - rand[ 1 ] - rand[ 0 ];
}

hmc_float calc_delta (__global hmc_ocl_ran* rnd , hmc_float alpha){
  return  (-log(ocl_new_ran(rnd))/alpha*pow(cos(2. * PI * ocl_new_ran(rnd)), 2.) -log(ocl_new_ran(rnd))/alpha );
}

hmc_float calc_eta (__global hmc_ocl_ran* rnd ){
  return  (ocl_new_ran(rnd) );
}

// Construct new SU2 matrix using improved alg by Kennedy Pendleton
void SU2Update(__private hmc_float dst [su2_entries], const hmc_float alpha, __global hmc_ocl_ran * rnd)
{
  hmc_float delta;
  hmc_float a0 ;
  hmc_float eta ; 
  do
  {
//     delta = -log(ocl_new_ran(rnd))/alpha*pow(cos(2. * PI * ocl_new_ran(rnd)), 2.) -log(ocl_new_ran(rnd))/alpha;
    delta = calc_delta(rnd, alpha);
    a0 = 1.-delta;
//     eta = ocl_new_ran(rnd);
    eta = calc_eta(rnd);
  }while ( (1.-0.5*delta) < eta*eta); 
  hmc_float phi = 2.*PI*ocl_new_ran(rnd);
  hmc_float theta = asin(2.*ocl_new_ran(rnd) - 1.);
  dst[0] = a0;
  dst[1] = sqrt(1.-a0 * a0)*cos(theta) * cos(phi);
  dst[2] = sqrt(1.-a0 * a0)*cos(theta) * sin(phi);
  dst[3] = sqrt(1.-a0 * a0)*sin(theta);
}

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


void calc_staple(__global hmc_ocl_gaugefield* field,__private hmc_ocl_staplematrix* dest, const int pos, const int t, const int mu_in){
  hmc_ocl_su3matrix prod[SU3SIZE];
  hmc_ocl_su3matrix prod2[SU3SIZE];
  hmc_ocl_su3matrix tmp[SU3SIZE];
  hmc_ocl_staplematrix dummy[STAPLEMATRIXSIZE];
  int nu, newpos, newt;
  
  zero_staplematrix(dummy);
  
  //iterate through the three directions other than mu
  for(int i = 1; i<NDIM; i++){
    
    zero_su3matrix(prod);
    nu = (mu_in + i)%NDIM;
    //first staple
    //u_nu(x+mu)
    if(mu_in==0) {
      newt = (t+1)%NTIME;
      get_su3matrix(tmp,field,pos,newt,nu);
    } else {
      get_su3matrix(tmp,field,get_neighbor(pos,mu_in),t,nu);
    }
//     get_su3matrix(tmp,field,1,0,0);
    
    copy_su3matrix(prod, tmp);
    //adjoint(u_mu(x+nu))
    if(nu==0) {
      newt = (t+1)%NTIME;
      get_su3matrix(tmp,field,pos,newt,mu_in);
    } else {
      get_su3matrix(tmp,field,get_neighbor(pos,nu),t,mu_in);
    }
    adjoin_su3matrix(tmp);
    accumulate_su3matrix_prod(prod,tmp);
    //adjoint(u_nu(x))
    get_su3matrix(tmp,field,pos,t,nu);
    adjoin_su3matrix(tmp);
    accumulate_su3matrix_prod(prod,tmp);  
    //second staple
    //adjoint (u_nu(x+mu-nu))
    //newpos is "pos-nu" (spatial)
    newpos = get_lower_neighbor(pos, nu);
    if(mu_in==0) {
      newt = (t+1)%NTIME;
      get_su3matrix(tmp,field,newpos,newt,nu);
    } else if (nu == 0){
        newt = (t-1+NTIME)%NTIME;
        get_su3matrix(tmp,field,get_neighbor(pos,mu_in),newt,nu);
    }
    else{
	get_su3matrix(tmp,field,get_neighbor(newpos,mu_in),t,nu);
    } 
    adjoin_su3matrix(tmp);
    copy_su3matrix(prod2, tmp);
    //adjoint(u_mu(x-nu))
    if(mu_in==0) {
      get_su3matrix(tmp,field,newpos,t,mu_in);
    } 
    else if (nu == 0){
        newt = (t-1+NTIME)%NTIME;
        get_su3matrix(tmp,field,pos,newt,mu_in);
    }
    else{
	get_su3matrix(tmp,field,newpos,t,mu_in);
    }
    adjoin_su3matrix(tmp);
    accumulate_su3matrix_prod(prod2,tmp);
    //adjoint(u_nu(x-nu))
    if(mu_in==0) {
      get_su3matrix(tmp,field,newpos,t,nu);
    } 
    else if (nu == 0){
        newt = (t-1+NTIME)%NTIME;
        get_su3matrix(tmp,field,pos,newt,nu);
    }
    else{
	get_su3matrix(tmp,field,newpos,t,nu);
    }
    accumulate_su3matrix_prod(prod2,tmp); 
   
    accumulate_su3matrices_add(dummy, prod);
    accumulate_su3matrices_add(dummy, prod2);
  }
  copy_staplematrix(dest, dummy);
}

__kernel void heatbath_even(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd){
  //it is assumed that there are VOL4D/2 threads
  int t, pos;
  int id = get_global_id(0);
  get_even_site(id, &pos, &t);
  
  hmc_ocl_su3matrix U[SU3SIZE];
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
  hmc_ocl_staplematrix staplematrix[STAPLEMATRIXSIZE];
  
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
  random_1_2_3(order, &rnd[id]);
  get_su3matrix(U, gaugefield, pos, t, mu);
  
  hmc_complex det = det_su3matrix(U);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(U);

  calc_staple(gaugefield, staplematrix, pos, t, mu);

  for(int i=0; i<3; i++)
  {
    multiply_staplematrix(W, U, staplematrix); 
    reduction(w, W, order[i]);

    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );

    hmc_float beta_neu =  2.*beta / NC*k;
    SU2Update(r_pauli, beta_neu, &rnd[id]);
    
    w[0].re = (r_pauli[0]*w_pauli[0] + r_pauli[1]*w_pauli[1] + r_pauli[2]*w_pauli[2] + r_pauli[3]*w_pauli[3] )/k;
    w[0].im = (w_pauli[0]*r_pauli[3] - w_pauli[3]*r_pauli[0] + r_pauli[1]*w_pauli[2] - r_pauli[2]*w_pauli[1] )/k;
    w[1].re = (w_pauli[0]*r_pauli[2] - w_pauli[2]*r_pauli[0] + r_pauli[3]*w_pauli[1] - r_pauli[1]*w_pauli[3] )/k;
    w[1].im = (w_pauli[0]*r_pauli[1] - w_pauli[1]*r_pauli[0] + r_pauli[2]*w_pauli[3] - r_pauli[3]*w_pauli[2] )/k;
    w[2].re = -(w_pauli[0]*r_pauli[2] - w_pauli[2]*r_pauli[0] + r_pauli[3]*w_pauli[1] - r_pauli[1]*w_pauli[3] )/k;
    w[2].im = (w_pauli[0]*r_pauli[1] - w_pauli[1]*r_pauli[0] + r_pauli[2]*w_pauli[3] - r_pauli[3]*w_pauli[2] )/k;
    w[3].re = (r_pauli[0]*w_pauli[0] + r_pauli[1]*w_pauli[1] + r_pauli[2]*w_pauli[2] + r_pauli[3]*w_pauli[3] )/k;
    w[3].im = -(w_pauli[0]*r_pauli[3] - w_pauli[3]*r_pauli[0] + r_pauli[1]*w_pauli[2] - r_pauli[2]*w_pauli[1] )/k;
    
    hmc_ocl_su3matrix extW[SU3SIZE]; 
    extend (extW, order[i], w); 
    //perhaps one can define another acc-func here and save one copying step
    accumulate_su3matrix_prod(extW, U);
    copy_su3matrix(U, extW);
  }
  put_su3matrix(gaugefield, U, pos, t, mu);

  return;
}

__kernel void heatbath_odd(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd){

  int t, pos;
  int id = get_global_id(0);
  get_odd_site(id, &pos, &t);
  
  hmc_ocl_su3matrix U[SU3SIZE];
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
  hmc_ocl_staplematrix staplematrix[STAPLEMATRIXSIZE];
  
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
  random_1_2_3(order, &rnd[id]);
  //old U
  get_su3matrix(U, gaugefield, pos, t, mu);
     
  hmc_complex det = det_su3matrix(U);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(U);

  calc_staple(gaugefield, staplematrix, pos, t, mu);

  for(int i=0; i<3; i++)
  {
    multiply_staplematrix(W, U, staplematrix); 
    reduction(w, W, order[i]);

    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );

    hmc_float beta_neu =  2.*beta / NC*k;
    SU2Update(r_pauli, beta_neu, &rnd[id]);
    
    w[0].re = (r_pauli[0]*w_pauli[0] + r_pauli[1]*w_pauli[1] + r_pauli[2]*w_pauli[2] + r_pauli[3]*w_pauli[3] )/k;
    w[0].im = (w_pauli[0]*r_pauli[3] - w_pauli[3]*r_pauli[0] + r_pauli[1]*w_pauli[2] - r_pauli[2]*w_pauli[1] )/k;
    w[1].re = (w_pauli[0]*r_pauli[2] - w_pauli[2]*r_pauli[0] + r_pauli[3]*w_pauli[1] - r_pauli[1]*w_pauli[3] )/k;
    w[1].im = (w_pauli[0]*r_pauli[1] - w_pauli[1]*r_pauli[0] + r_pauli[2]*w_pauli[3] - r_pauli[3]*w_pauli[2] )/k;
    w[2].re = -(w_pauli[0]*r_pauli[2] - w_pauli[2]*r_pauli[0] + r_pauli[3]*w_pauli[1] - r_pauli[1]*w_pauli[3] )/k;
    w[2].im = (w_pauli[0]*r_pauli[1] - w_pauli[1]*r_pauli[0] + r_pauli[2]*w_pauli[3] - r_pauli[3]*w_pauli[2] )/k;
    w[3].re = (r_pauli[0]*w_pauli[0] + r_pauli[1]*w_pauli[1] + r_pauli[2]*w_pauli[2] + r_pauli[3]*w_pauli[3] )/k;
    w[3].im = -(w_pauli[0]*r_pauli[3] - w_pauli[3]*r_pauli[0] + r_pauli[1]*w_pauli[2] - r_pauli[2]*w_pauli[1] )/k;
    
    hmc_ocl_su3matrix extW[SU3SIZE]; 
    extend (extW, order[i], w); 
    //perhaps one can define another acc-func here and save one copying step
    accumulate_su3matrix_prod(extW, U);
    copy_su3matrix(U, extW);
  }
  put_su3matrix(gaugefield, U, pos, t, mu);

  return;
}

__kernel void overrelax_even(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd){
  //it is assumed that there are VOL4D/2 threads
  int t, pos;
  int id = get_global_id(0);
  get_even_site(id, &pos, &t);
  
  hmc_ocl_su3matrix U[SU3SIZE];
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
  hmc_ocl_staplematrix staplematrix[STAPLEMATRIXSIZE];
  
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k;
  int order[3]; 
  random_1_2_3(order, &rnd[id]);
  get_su3matrix(U, gaugefield, pos, t, mu);

  hmc_complex det = det_su3matrix(U);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(U);
  
  calc_staple(gaugefield, staplematrix, pos, t, mu);
  hmc_ocl_su3matrix tmp[SU3SIZE]; 
  hmc_ocl_su3matrix extW[SU3SIZE]; 
  for(int i=0; i<3; i++)
  {
    multiply_staplematrix(W, U, staplematrix); 
    reduction(w, W, order[i]);

    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );

    w[0].re = (w_pauli[0]*w_pauli[0] - w_pauli[1]*w_pauli[1] - w_pauli[2]*w_pauli[2] - w_pauli[3]*w_pauli[3])/k/k;
    w[0].im = (-2.*w_pauli[0]*w_pauli[3])/k/k;
    w[1].re = (-2.*w_pauli[0]*w_pauli[2])/k/k;
    w[1].im = (-2.*w_pauli[0]*w_pauli[1])/k/k;
    w[2].re = (2.*w_pauli[0]*w_pauli[2])/k/k;
    w[2].im = (-2.*w_pauli[0]*w_pauli[1])/k/k;
    w[3].re = (w_pauli[0]*w_pauli[0] - w_pauli[1]*w_pauli[1] - w_pauli[2]*w_pauli[2] - w_pauli[3]*w_pauli[3])/k/k;
    w[3].im = (2.*w_pauli[0]*w_pauli[3])/k/k;
    
    extend (extW, order[i], w); 
    multiply_su3matrices(tmp, extW, U);
    copy_su3matrix(U, tmp);
  }
  put_su3matrix(gaugefield, U, pos, t, mu);
  return;
}

__kernel void overrelax_odd(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd){
  int t, pos;
  int id = get_global_id(0);
  get_odd_site(id, &pos, &t);
  
  hmc_ocl_su3matrix U[SU3SIZE];
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
  hmc_ocl_staplematrix staplematrix[STAPLEMATRIXSIZE];
  
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k;
  int order[3]; 
  random_1_2_3(order, &rnd[id]);
  get_su3matrix(U, gaugefield, pos, t, mu);

  hmc_complex det = det_su3matrix(U);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(U);
  
  calc_staple(gaugefield, staplematrix, pos, t, mu);
  hmc_ocl_su3matrix tmp[SU3SIZE]; 
  hmc_ocl_su3matrix extW[SU3SIZE]; 
  for(int i=0; i<3; i++)
  {
    multiply_staplematrix(W, U, staplematrix); 
    reduction(w, W, order[i]);

    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );

    w[0].re = (w_pauli[0]*w_pauli[0] - w_pauli[1]*w_pauli[1] - w_pauli[2]*w_pauli[2] - w_pauli[3]*w_pauli[3])/k/k;
    w[0].im = (-2.*w_pauli[0]*w_pauli[3])/k/k;
    w[1].re = (-2.*w_pauli[0]*w_pauli[2])/k/k;
    w[1].im = (-2.*w_pauli[0]*w_pauli[1])/k/k;
    w[2].re = (2.*w_pauli[0]*w_pauli[2])/k/k;
    w[2].im = (-2.*w_pauli[0]*w_pauli[1])/k/k;
    w[3].re = (w_pauli[0]*w_pauli[0] - w_pauli[1]*w_pauli[1] - w_pauli[2]*w_pauli[2] - w_pauli[3]*w_pauli[3])/k/k;
    w[3].im = (2.*w_pauli[0]*w_pauli[3])/k/k;
    
    extend (extW, order[i], w); 
    //perhaps one can define another acc-func here and save one copying step
    multiply_su3matrices(tmp, extW, U);
    copy_su3matrix(U, tmp);
  }
  put_su3matrix(gaugefield, U, pos, t, mu);
  return;
}
__kernel void plaquette(__global hmc_ocl_gaugefield * field,__global hmc_float * plaq_out, __global hmc_float* tplaq, __global hmc_float* splaq){
  int t, pos;
  int id = get_global_id(0);

  hmc_float plaq=0;
  hmc_float splaq_tmp=0;
  hmc_float tplaq_tmp=0;
  hmc_float tmpfloat = 0;

  hmc_ocl_su3matrix tmp[SU3SIZE]; 
  hmc_ocl_su3matrix prod[SU3SIZE];
  
  //this is copied from the host-code (see comments there)
  //calc even plaquette
  get_even_site(id, &pos, &t);
  for(int mu=0; mu<NDIM; mu++) {
    for(int nu=0;nu<mu; nu++) {
      get_su3matrix(prod,field,pos,t,mu);
      if(mu==0) {
	int newt = (t+1)%NTIME;
	get_su3matrix(tmp,field,pos,newt,nu);
	} else {
	    get_su3matrix(tmp,field,get_neighbor(pos,mu),t,nu);
	}
	accumulate_su3matrix_prod(prod,tmp);
	if(nu==0) {
	  int newt = (t+1)%NTIME;
	  get_su3matrix(tmp,field,pos,newt,mu);
	} else {
	  get_su3matrix(tmp,field,get_neighbor(pos,nu),t,mu);
	}
	adjoin_su3matrix(tmp);
	accumulate_su3matrix_prod(prod,tmp);
	get_su3matrix(tmp,field,pos,t,nu);
	adjoin_su3matrix(tmp);
	accumulate_su3matrix_prod(prod,tmp);
	tmpfloat = trace_su3matrix(prod).re;
	plaq += tmpfloat;
	if(mu==0 || nu==0) {
	  tplaq_tmp+=tmpfloat;
	} else {
	  splaq_tmp+=tmpfloat;
	}
  }}
    
  //calc odd plaquette
  get_odd_site(id, &pos, &t);
  for(int mu=0; mu<NDIM; mu++) {
    for(int nu=0;nu<mu; nu++) {
      get_su3matrix(prod,field,pos,t,mu);
      if(mu==0) {
	int newt = (t+1)%NTIME;
	get_su3matrix(tmp,field,pos,newt,nu);
	} else {
	    get_su3matrix(tmp,field,get_neighbor(pos,mu),t,nu);
	}
	accumulate_su3matrix_prod(prod,tmp);
	if(nu==0) {
	  int newt = (t+1)%NTIME;
	  get_su3matrix(tmp,field,pos,newt,mu);
	} else {
	  get_su3matrix(tmp,field,get_neighbor(pos,nu),t,mu);
	}
	adjoin_su3matrix(tmp);
	accumulate_su3matrix_prod(prod,tmp);
	get_su3matrix(tmp,field,pos,t,nu);
	adjoin_su3matrix(tmp);
	accumulate_su3matrix_prod(prod,tmp);
	tmpfloat = trace_su3matrix(prod).re;
	plaq += tmpfloat;
	if(mu==0 || nu==0) {
	  tplaq_tmp+=tmpfloat;
	} else {
	  splaq_tmp+=tmpfloat;
	}
  }}
  
  (*plaq_out) += plaq;
  (*splaq) += splaq_tmp;
  (*tplaq) += tplaq_tmp;

}


__kernel void polyakov(__global hmc_ocl_gaugefield * field, __global hmc_complex * out){
  int const tdir = 0;
  int pos, t;
  int id = get_global_id(0);
   
  hmc_ocl_su3matrix prod[SU3SIZE];
  hmc_ocl_su3matrix tmp[SU3SIZE];
  unit_su3matrix(prod);
  hmc_complex tmpcomplex;
  
  get_even_site(id, &pos, &t);
  if(t==0){
    for(int t=0; t<NTIME; t++) {
      get_su3matrix(tmp,field,pos,t,tdir);
      accumulate_su3matrix_prod(prod,tmp);
    }
    tmpcomplex = trace_su3matrix(prod);
    (*out).re += tmpcomplex.re;
    (*out).im += tmpcomplex.im;
  }
  else if(t==1){
    get_odd_site(id, &pos, &t);
    for(int t=0; t<NTIME; t++) {
      get_su3matrix(tmp,field,pos,t,tdir);
      accumulate_su3matrix_prod(prod,tmp);
    }
    tmpcomplex = trace_su3matrix(prod);
    (*out).re += tmpcomplex.re;
    (*out).im += tmpcomplex.im;
  }
  else return;
}

//EOF
