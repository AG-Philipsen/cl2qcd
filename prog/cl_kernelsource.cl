//opencl_header.cl

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

/*
#define SPINORSIZE NSPIN*NC
#define HALFSPINORSIZE NSPIN/2*NC
#define SPINORFIELDSIZE SPINORSIZE*NTIME*VOLSPACE
#define EOPREC_SPINORFIELDSIZE SPINORSIZE*NTIME*VOLSPACE/2
*/

#define VOL4D VOLSPACE*NTIME

//opencl_geometry.cl

int inline get_global_pos(int spacepos, int t){
  return spacepos + VOLSPACE * t;
}

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
  t = (int)(idx/(VOLSPACE/2));
  x -= t*VOLSPACE/2;
  z = (int)(x/(NSPACE*NSPACE/2));
  x -= z*NSPACE*NSPACE/2;
  y = (int)(x/NSPACE);
  x -= y*NSPACE;
  (*out_space) =  (int)((z+t)%2)*(1 + 2*x - (int) (2*x/NSPACE)) + (int)((t+z+1)%2)*(2*x + (int) (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
  (*out_t) = t;
}

//it is assumed that idx iterates only over half the number of sites
void inline get_odd_site(int idx, int * out_space, int * out_t){
  int x,y,z,t;
  x = idx;
  t = (int)(idx/(VOLSPACE/2));
  x -= t*VOLSPACE/2;
  z = (int)(x/(NSPACE*NSPACE/2));
  x -= z*NSPACE*NSPACE/2;
  y = (int)(x/NSPACE);
  x -= y*NSPACE;
  (*out_space) =  (int)((z+t+1)%2)*(1 + 2*x - (int) (2*x/NSPACE)) + (int)((t+z)%2)*(2*x + (int) (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
  (*out_t) = t;
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
  get_allspacecoord(nspace, coord);
  coord[dir] = (coord[dir] + 1)%NSPACE;
  return get_nspace(coord);
}

int get_lower_neighbor(const int nspace, int const dir) {
  int coord[NDIM];
  get_allspacecoord(nspace, coord);
  coord[dir] = (coord[dir] - 1 + NSPACE)%NSPACE;
  return get_nspace(coord);
}

int spinor_color(int spinor_element){
  return (int)(spinor_element/NSPIN);
}

int spinor_spin(int spinor_element,int color){
  return spinor_element - NSPIN*color;
}

int spinor_element(int alpha, int color) {
  return alpha + NSPIN*color;
}

int get_n_eoprec(int timepos, int spacepos){
  return (int)((get_global_pos(spacepos, timepos))/2);
}

int eoprec_spinor_field_element(int alpha, int color, int n_eoprec) {
  return alpha + NSPIN*color + NSPIN*NC*n_eoprec;
}

int spinor_field_element(int alpha, int color, int nspace, int t) {
  return alpha + NSPIN*color + NSPIN*NC*(get_global_pos(nspace, t));
}
//opencl_operations_complex.cl


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

hmc_complex complexdivide(hmc_complex* numerator, hmc_complex* denominator){
  hmc_float norm = (*denominator).re*(*denominator).re + (*denominator).im*(*denominator).im;
  hmc_complex res;
  res.re = ((*numerator).re*(*denominator).re+(*numerator).im*(*denominator).im)/norm;
  res.im = ((*numerator).im*(*denominator).re-(*numerator).re*(*denominator).im)/norm;
  return res;
}

//opencl_operations_gaugefield.cl

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

void gaugefield_apply_bc(__private hmc_ocl_su3matrix * in, hmc_float theta){
  hmc_float tmp1,tmp2;
#ifdef _RECONSTRUCT_TWELVE_
  hmc_ocl_su3matrix tmp;
  copy_su3matrix(tmp, mat);
  for(int n=0; n<NC*(NC-1); n++) {
    tmp1 = ((in)[ocl_su3matrix_element(a,b)]).re;
    tmp2 = ((in)[ocl_su3matrix_element(a,b)]).im;
    ((in)[ocl_su3matrix_element(a,b)]).re = cos(theta)*tmp1 - sin(theta)*tmp2;
    ((in)[ocl_su3matrix_element(a,b)]).im = sin(theta)*tmp1 + cos(theta)*tmp2;
  }
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      tmp1 = ((in)[ocl_su3matrix_element(a,b)]).re;
      tmp2 = ((in)[ocl_su3matrix_element(a,b)]).im;
      ((in)[ocl_su3matrix_element(a,b)]).re = cos(theta)*tmp1 - sin(theta)*tmp2;
      ((in)[ocl_su3matrix_element(a,b)]).im = sin(theta)*tmp1 + cos(theta)*tmp2;
    }
  }
#endif
  return;
}

// replace link in with e^(mu.re)*(cos(mu.im) + i*sin(mu.im))
void gaugefield_apply_chem_pot(__private hmc_ocl_su3matrix * u, __private hmc_ocl_su3matrix * udagger, hmc_float chem_pot_re, hmc_float chem_pot_im){
  hmc_float tmp1,tmp2;
#ifdef _RECONSTRUCT_TWELVE_
  hmc_ocl_su3matrix tmp;
  copy_su3matrix(tmp, mat);
  for(int n=0; n<NC*(NC-1); n++) {
    tmp1 = ((u)[ocl_su3matrix_element(a,b)]).re;
    tmp2 = ((u)[ocl_su3matrix_element(a,b)]).im;
    ((u)[ocl_su3matrix_element(a,b)]).re = exp(chem_pot_re)*cos(chem_pot_im)*tmp1;
    ((u)[ocl_su3matrix_element(a,b)]).im = exp(chem_pot_re)*sin(chem_pot_im)*tmp2;
    tmp1 = ((udagger)[ocl_su3matrix_element(a,b)]).re;
    tmp2 = ((udagger)[ocl_su3matrix_element(a,b)]).im;
    ((udagger)[ocl_su3matrix_element(a,b)]).re = exp(-chem_pot_re)*cos(chem_pot_im)*tmp1;
    ((udagger)[ocl_su3matrix_element(a,b)]).im = -exp(-chem_pot_re)*sin(chem_pot_im)*tmp2;
  }
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      tmp1 = ((u)[ocl_su3matrix_element(a,b)]).re;
      tmp2 = ((u)[ocl_su3matrix_element(a,b)]).im;
      ((u)[ocl_su3matrix_element(a,b)]).re = exp(chem_pot_re)*cos(chem_pot_im)*tmp1;
      ((u)[ocl_su3matrix_element(a,b)]).im = exp(chem_pot_re)*sin(chem_pot_im)*tmp2;
      tmp1 = ((udagger)[ocl_su3matrix_element(a,b)]).re;
      tmp2 = ((udagger)[ocl_su3matrix_element(a,b)]).im;
      ((udagger)[ocl_su3matrix_element(a,b)]).re = exp(-chem_pot_re)*cos(chem_pot_im)*tmp1;
      ((udagger)[ocl_su3matrix_element(a,b)]).im = -exp(-chem_pot_re)*sin(chem_pot_im)*tmp2;
    }
  }
#endif
  return;
}
  
  //opencl_operations_spinor

void set_zero_spinor(hmc_spinor_field *field) {
  for (int n=0; n<SPINORFIELDSIZE; n++) {
	field[n].re=0;
	field[n].im=0;
  }
  return;
}

void set_zero_eoprec_spinor(hmc_eoprec_spinor_field *field) {
  for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
	field[n].re=0;
	field[n].im=0;
  }
  return;
}

hmc_float local_squarenorm(hmc_spinor_field *field, int spacepos, int timepos) {
  hmc_float sum=0;
  for (int a=0; a<NSPIN; a++) {
    for(int c=0; c<NC; c++) {
      hmc_float dummy_re = field[spinor_field_element(a,c,spacepos,timepos)].re;
      hmc_float dummy_im = field[spinor_field_element(a,c,spacepos,timepos)].im;
      sum += dummy_re*dummy_re + dummy_im*dummy_im;
    }
  }
  return sum;
}

hmc_float global_squarenorm(hmc_spinor_field *field) {
  hmc_float sum=0;
  for (int t=0; t<NTIME; t++) {
    for (int n=0; n<VOLSPACE; n++) {
      sum += local_squarenorm(field,n,t);
    }
  }
  return sum;
}

void fill_with_one(hmc_spinor_field *field, int spacepos, int timepos, int alpha, int color){
  field[spinor_field_element(alpha,color,spacepos,timepos)].re = hmc_one_f;
  field[spinor_field_element(alpha,color,spacepos,timepos)].im = 0;
  return;
}

void su3matrix_times_colorvector(hmc_ocl_su3matrix* u, hmc_color_vector* in, hmc_color_vector* out){
#ifdef _RECONSTRUCT_TWELVE_
  for(int a=0; a<NC-1; a++) {
    out[a] = hmc_complex_zero;
    for(int b=0; b<NC; b++) {
      hmc_complex tmp = complexmult(&((*u)[ocl_su3matrix_element(a,b)]),&(in[b]));
      complexaccumulate(&out[a],&tmp);
    }
  }
  out[2] = hmc_complex_zero;
  for(int b=0; b<NC; b++) {
    hmc_complex rec = reconstruct_su3(u,b);
    hmc_complex tmp = complexmult(&rec,&(in[b]));
    complexaccumulate(&out[2],&tmp);
  }
#else
  for(int a=0; a<NC; a++) {
    out[a] = hmc_complex_zero;
    for(int b=0; b<NC; b++) {
      hmc_complex tmp = complexmult(&((u)[ocl_su3matrix_element(a,b)]),&(in[b]));
      complexaccumulate(&(out[a]),&tmp);
    }
  }
#endif
  return;
}

void set_local_zero_spinor(hmc_spinor* inout){
  for(int j=0; j<SPINORSIZE; j++) {
    inout[j].re = 0;
    inout[j].im = 0;
  }
  return;
}

hmc_float spinor_squarenorm(hmc_spinor* in){
  hmc_float res=0;
  for(int j=0; j<SPINORSIZE; j++) {
    hmc_complex tmp = complexconj(&in[j]);
    hmc_complex incr= complexmult(&tmp,&in[j]);
    res+=incr.re;
  }
  return res;
}

void real_multiply_spinor(hmc_spinor* inout, hmc_float factor){
  for(int j=0; j<SPINORSIZE; j++) inout[j].re *=factor;
  return;
}

void spinprojectproduct_gamma0(hmc_ocl_su3matrix* u, hmc_spinor* spin,hmc_float sign){
  //out = u*(1+sign*gamma0)*in
  hmc_color_vector vec1[NC];
  hmc_color_vector vec2[NC];
  hmc_color_vector uvec1[NC];
  hmc_color_vector uvec2[NC];
  for(int c=0; c<NC; c++) {
    vec1[c].re = spin[spinor_element(0,c)].re - sign*spin[spinor_element(3,c)].im;
    vec1[c].im = spin[spinor_element(0,c)].im + sign*spin[spinor_element(3,c)].re;
    vec2[c].re = spin[spinor_element(1,c)].re - sign*spin[spinor_element(2,c)].im;
    vec2[c].im = spin[spinor_element(1,c)].im + sign*spin[spinor_element(2,c)].re;
  }
  su3matrix_times_colorvector(u,vec1,uvec1);
  su3matrix_times_colorvector(u,vec2,uvec2);
  for(int c=0; c<NC; c++) {
    spin[spinor_element(0,c)].re = uvec1[c].re;
    spin[spinor_element(0,c)].im = uvec1[c].im;
    spin[spinor_element(1,c)].re = uvec2[c].re;
    spin[spinor_element(1,c)].im = uvec2[c].im;
    spin[spinor_element(2,c)].re = sign*uvec2[c].im;
    spin[spinor_element(2,c)].im = -sign*uvec2[c].re;
    spin[spinor_element(3,c)].re = sign*uvec1[c].im;
    spin[spinor_element(3,c)].im = -sign*uvec1[c].re;
  }
  return;
}

void spinprojectproduct_gamma1(hmc_ocl_su3matrix* u, hmc_spinor* spin, hmc_float sign){
  //out = u*(1+sign*gamma1)*in
  hmc_color_vector vec1[NC];
  hmc_color_vector vec2[NC];
  hmc_color_vector uvec1[NC];
  hmc_color_vector uvec2[NC];
  for(int c=0; c<NC; c++) {
    vec1[c].re = spin[spinor_element(0,c)].re - sign*spin[spinor_element(3,c)].re;
    vec1[c].im = spin[spinor_element(0,c)].im - sign*spin[spinor_element(3,c)].im;
    vec2[c].re = spin[spinor_element(1,c)].re + sign*spin[spinor_element(2,c)].re;
    vec2[c].im = spin[spinor_element(1,c)].im + sign*spin[spinor_element(2,c)].im;
  }
  su3matrix_times_colorvector(u,vec1,uvec1);
  su3matrix_times_colorvector(u,vec2,uvec2);
  for(int c=0; c<NC; c++) {
    spin[spinor_element(0,c)].re = uvec1[c].re;
    spin[spinor_element(0,c)].im = uvec1[c].im;
    spin[spinor_element(1,c)].re = uvec2[c].re;
    spin[spinor_element(1,c)].im = uvec2[c].im;
    spin[spinor_element(2,c)].re = sign*uvec2[c].re;
    spin[spinor_element(2,c)].im = sign*uvec2[c].im;
    spin[spinor_element(3,c)].re = -sign*uvec1[c].re;
    spin[spinor_element(3,c)].im = -sign*uvec1[c].im;
  }
  return;
}

void spinprojectproduct_gamma2(hmc_ocl_su3matrix* u, hmc_spinor* spin, hmc_float sign){
  //out = u*(1+sign*gamma2)*in
  hmc_color_vector vec1[NC];
  hmc_color_vector vec2[NC];
  hmc_color_vector uvec1[NC];
  hmc_color_vector uvec2[NC];
  for(int c=0; c<NC; c++) {
    vec1[c].re = spin[spinor_element(0,c)].re - sign*spin[spinor_element(2,c)].im;
    vec1[c].im = spin[spinor_element(0,c)].im + sign*spin[spinor_element(2,c)].re;
    vec2[c].re = spin[spinor_element(1,c)].re + sign*spin[spinor_element(3,c)].im;
    vec2[c].im = spin[spinor_element(1,c)].im - sign*spin[spinor_element(3,c)].re;
  }
  su3matrix_times_colorvector(u,vec1,uvec1);
  su3matrix_times_colorvector(u,vec2,uvec2);
  for(int c=0; c<NC; c++) {
    spin[spinor_element(0,c)].re = uvec1[c].re;
    spin[spinor_element(0,c)].im = uvec1[c].im;
    spin[spinor_element(1,c)].re = uvec2[c].re;
    spin[spinor_element(1,c)].im = uvec2[c].im;
    spin[spinor_element(2,c)].re = sign*uvec1[c].im;
    spin[spinor_element(2,c)].im = -sign*uvec1[c].re;
    spin[spinor_element(3,c)].re = -sign*uvec2[c].im;
    spin[spinor_element(3,c)].im = sign*uvec2[c].re;
  }
  return;
}

void spinprojectproduct_gamma3(hmc_ocl_su3matrix* u, hmc_spinor* spin, hmc_float sign){
  //out = u*(1+sign*gamma3)*in
  hmc_color_vector vec1[NC];
  hmc_color_vector vec2[NC];
  hmc_color_vector uvec1[NC];
  hmc_color_vector uvec2[NC];
  for(int c=0; c<NC; c++) {
    vec1[c].re = spin[spinor_element(0,c)].re + sign*spin[spinor_element(2,c)].re;
    vec1[c].im = spin[spinor_element(0,c)].im + sign*spin[spinor_element(2,c)].im;
    vec2[c].re = spin[spinor_element(1,c)].re + sign*spin[spinor_element(3,c)].re;
    vec2[c].im = spin[spinor_element(1,c)].im + sign*spin[spinor_element(3,c)].im;
  }
  su3matrix_times_colorvector(u,vec1,uvec1);
  su3matrix_times_colorvector(u,vec2,uvec2);
  for(int c=0; c<NC; c++) {
    spin[spinor_element(0,c)].re = uvec1[c].re;
    spin[spinor_element(0,c)].im = uvec1[c].im;
    spin[spinor_element(1,c)].re = uvec2[c].re;
    spin[spinor_element(1,c)].im = uvec2[c].im;
    spin[spinor_element(2,c)].re = sign*uvec1[c].re;
    spin[spinor_element(2,c)].im = sign*uvec1[c].im;
    spin[spinor_element(3,c)].re = sign*uvec2[c].re;
    spin[spinor_element(3,c)].im = sign*uvec2[c].im;
  }
  return;
}

void spinors_accumulate(hmc_spinor* inout, hmc_spinor* incr){
  for(int j=0; j<SPINORSIZE; j++) {
    inout[j].re += incr[j].re;
    inout[j].im += incr[j].im;
  }
  return;
}

void multiply_spinor_factor_gamma5(hmc_spinor* in, hmc_spinor* out, hmc_float factor){
  for(int c=0; c<NC; c++) {
    out[spinor_element(0,c)].re = -factor*in[spinor_element(0,c)].im;
    out[spinor_element(1,c)].re = -factor*in[spinor_element(1,c)].im;
    out[spinor_element(2,c)].re = factor*in[spinor_element(2,c)].im;
    out[spinor_element(3,c)].re = factor*in[spinor_element(3,c)].im;
    out[spinor_element(0,c)].im = factor*in[spinor_element(0,c)].re;
    out[spinor_element(1,c)].im = factor*in[spinor_element(1,c)].re;
    out[spinor_element(2,c)].im = -factor*in[spinor_element(2,c)].re;
    out[spinor_element(3,c)].im = -factor*in[spinor_element(3,c)].re;
  }
  return;
}

void multiply_spinor_gamma0(hmc_spinor* in,hmc_spinor* out){
  for(int c=0; c<NC; c++) {
    out[spinor_element(0,c)].re = -in[spinor_element(3,c)].im;
    out[spinor_element(1,c)].re = -in[spinor_element(2,c)].im;
    out[spinor_element(2,c)].re = in[spinor_element(1,c)].im;
    out[spinor_element(3,c)].re = in[spinor_element(0,c)].im;
    out[spinor_element(0,c)].im = in[spinor_element(3,c)].re;
    out[spinor_element(1,c)].im = in[spinor_element(2,c)].re;
    out[spinor_element(2,c)].im = -in[spinor_element(1,c)].re;
    out[spinor_element(3,c)].im = -in[spinor_element(0,c)].re;
  }
  return;
}
void multiply_spinor_gamma1(hmc_spinor* in,hmc_spinor* out){
  for(int c=0; c<NC; c++) {
    out[spinor_element(0,c)].re = -in[spinor_element(3,c)].re;
    out[spinor_element(1,c)].re = in[spinor_element(2,c)].re;
    out[spinor_element(2,c)].re = in[spinor_element(1,c)].re;
    out[spinor_element(3,c)].re = -in[spinor_element(0,c)].re;
    out[spinor_element(0,c)].im = -in[spinor_element(3,c)].im;
    out[spinor_element(1,c)].im = in[spinor_element(2,c)].im;
    out[spinor_element(2,c)].im = in[spinor_element(1,c)].im;
    out[spinor_element(3,c)].im = -in[spinor_element(0,c)].im;
  }
  return;
}
void multiply_spinor_gamma2(hmc_spinor* in,hmc_spinor* out){
  for(int c=0; c<NC; c++) {
    out[spinor_element(0,c)].re = -in[spinor_element(2,c)].im;
    out[spinor_element(1,c)].re = in[spinor_element(3,c)].im;
    out[spinor_element(2,c)].re = in[spinor_element(0,c)].im;
    out[spinor_element(3,c)].re = -in[spinor_element(1,c)].im;
    out[spinor_element(0,c)].im = in[spinor_element(2,c)].re;
    out[spinor_element(1,c)].im = -in[spinor_element(3,c)].re;
    out[spinor_element(2,c)].im = -in[spinor_element(0,c)].re;
    out[spinor_element(3,c)].im = in[spinor_element(1,c)].re;
  }
  return;
}
void multiply_spinor_gamma3(hmc_spinor* in,hmc_spinor* out){
  for(int c=0; c<NC; c++) {
    out[spinor_element(0,c)].re = in[spinor_element(2,c)].re;
    out[spinor_element(1,c)].re = in[spinor_element(3,c)].re;
    out[spinor_element(2,c)].re = in[spinor_element(0,c)].re;
    out[spinor_element(3,c)].re = in[spinor_element(1,c)].re;
    out[spinor_element(0,c)].im = in[spinor_element(2,c)].im;
    out[spinor_element(1,c)].im = in[spinor_element(3,c)].im;
    out[spinor_element(2,c)].im = in[spinor_element(0,c)].im;
    out[spinor_element(3,c)].im = in[spinor_element(1,c)].im;
  }
  return;
}

void su3matrix_times_spinor(hmc_ocl_su3matrix* u, hmc_spinor* in, hmc_spinor* out){
  for (int alpha=0; alpha<NSPIN; alpha++) {
    hmc_color_vector vec_in[NC];
    hmc_color_vector vec_out[NC];
    for(int c=0; c<NC; c++) {
      vec_in[c].re = in[spinor_element(alpha,c)].re;
      vec_in[c].im = in[spinor_element(alpha,c)].im;
    }
    su3matrix_times_colorvector(u, vec_in, vec_out);
    for(int c=0; c<NC; c++) {
      in[spinor_element(alpha,c)].re = vec_out[c].re;
      in[spinor_element(alpha,c)].im = vec_out[c].im;
    }
  }
  return;
}

//eoprec operations
void convert_to_eoprec(hmc_eoprec_spinor_field* even, hmc_eoprec_spinor_field* odd, hmc_spinor_field* in){
  int spacepos, timepos;
  for(int n=0; n<VOL4D/2; n++) {
    for(int alpha=0; alpha<NSPIN; alpha++) {
      for(int color=0; color<NC; color++) {
	get_even_site(n, &spacepos, &timepos);
	even[eoprec_spinor_field_element(alpha,color,n)] = in[spinor_field_element(alpha,color,spacepos,timepos)];
	get_odd_site(n, &spacepos, &timepos);
	odd[eoprec_spinor_field_element(alpha,color,n)] = in[spinor_field_element(alpha,color,spacepos,timepos)];
      }
    }
  }
  return;
}

void convert_from_eoprec(hmc_eoprec_spinor_field* even, hmc_eoprec_spinor_field* odd, hmc_spinor_field* out){
  int spacepos, timepos;
  for(int n=0; n<VOL4D/2; n++) {
    for(int alpha=0; alpha<NSPIN; alpha++) {
      for(int color=0; color<NC; color++) {
        get_even_site(n, &spacepos, &timepos);
        out[spinor_field_element(alpha,color,spacepos,timepos)] = even[eoprec_spinor_field_element(alpha,color,n)];
	get_odd_site(n, &spacepos, &timepos);
	out[spinor_field_element(alpha,color,spacepos,timepos)] = odd[eoprec_spinor_field_element(alpha,color,n)];
      }
    }
  }
  return;
}

void convert_to_kappa_format(__global hmc_spinor_field* inout,hmc_float kappa){
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    inout[n].re *= sqrt(2.*kappa);
    inout[n].im *= sqrt(2.*kappa);
  }
  return;
}

void convert_to_kappa_format_eoprec(__global hmc_eoprec_spinor_field* inout,hmc_float kappa){
  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    inout[n].re *= sqrt(2.*kappa);
    inout[n].im *= sqrt(2.*kappa);
  }
  return;
}

void convert_from_kappa_format(__global hmc_spinor_field* in, __global hmc_spinor_field * out,hmc_float kappa){
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    out[n].re = (in[n].re)/sqrt(2.*kappa);
    out[n].im = (in[n].im)/sqrt(2.*kappa);
  }
  return;
}

void convert_from_kappa_format_eoprec(__global hmc_eoprec_spinor_field* in, __global hmc_eoprec_spinor_field * out, hmc_float kappa){
  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    out[n].re = (in[n].re)/sqrt(2.*kappa);
    out[n].im = (in[n].im)/sqrt(2.*kappa);
  }
  return;
}

hmc_float global_squarenorm_eoprec(hmc_eoprec_spinor_field *in) {
  hmc_float sum=0;
  for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    sum += in[n].re*in[n].re + in[n].im*in[n].im;
  }
  return sum;
}

hmc_complex scalar_product_eoprec(hmc_eoprec_spinor_field* a, hmc_eoprec_spinor_field* b){
  // (a,b) = sum_k conj(a_k)*b_k
  hmc_complex res;
  res.re=0;
  res.im=0;
  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    res.re += a[n].re*b[n].re + a[n].im*b[n].im;
    res.im += a[n].re*b[n].im - a[n].im*b[n].re;
  }
  return res;
}

hmc_complex scalar_product(hmc_spinor_field* a, hmc_spinor_field* b){
  // (a,b) = sum_k conj(a_k)*b_k
  hmc_complex res;
  res.re=0;
  res.im=0;
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    res.re += a[n].re*b[n].re + a[n].im*b[n].im;
    res.im += a[n].re*b[n].im - a[n].im*b[n].re;
  }
  return res;
}

void get_spinor_from_eoprec_field(hmc_eoprec_spinor_field* in, hmc_spinor* out, int n_eoprec){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[spinor_element(alpha,color)] = in[eoprec_spinor_field_element(alpha,color,n_eoprec)];
    }
  }
  return;
}

void put_spinor_to_eoprec_field(hmc_spinor* in, hmc_eoprec_spinor_field* out, int n_eoprec){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[eoprec_spinor_field_element(alpha,color,n_eoprec)]=in[spinor_element(alpha,color)];
    }
  }
  return;
}

void get_spinor_from_field(hmc_spinor_field* in, hmc_spinor* out, int n, int t){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[spinor_element(alpha,color)] = in[spinor_field_element(alpha,color,n,t)];
    }
  }
  return;
}

void put_spinor_to_field(hmc_spinor* in, hmc_spinor_field* out, int n, int t){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[spinor_field_element(alpha,color,n,t)]=in[spinor_element(alpha,color)];
    }
  }
  return;
}

void copy_spinor(hmc_complex * in, hmc_complex * out){
  for (int n=0; n<SPINORFIELDSIZE; n++) {
    out[n].re = in[n].re;
    out[n].im = in[n].im;
  }
  return;
}

void copy_spinorfield(__global hmc_complex * in,__global  hmc_complex * out){
  for (int n=0; n<SPINORFIELDSIZE; n++) {
    out[n].re = in[n].re;
    out[n].im = in[n].im;
  }
  return;
}

void get_spinor(__global hmc_complex * in, hmc_complex * out){
  for (int n=0; n<SPINORFIELDSIZE; n++) {
    out[n].re = in[n].re;
    out[n].im = in[n].im;
  }
  return;
}

void put_spinor(hmc_complex * in, __global hmc_complex * out){
  for (int n=0; n<SPINORFIELDSIZE; n++) {
    out[n].re = in[n].re;
    out[n].im = in[n].im;
  }
  return;
}

// -alpha*x + y
//CP: defined with a minus!!!
void saxpy(hmc_spinor_field * x, hmc_spinor_field * y, hmc_complex * alpha, hmc_spinor_field * out){
  for (int n=0; n<SPINORFIELDSIZE; n++) {
    hmc_complex tmp1 = complexmult(alpha,&x[n]);
    ((out)[n]).re = -(tmp1).re + y[n].re;
    ((out)[n]).im = -(tmp1).im + y[n].im;
  }
  return;
}

void saxpy_eoprec(hmc_eoprec_spinor_field * x, hmc_eoprec_spinor_field * y, hmc_complex * alpha, hmc_eoprec_spinor_field * out){
  for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    hmc_complex tmp1 = complexmult(alpha,&x[n]);
    ((out)[n]).re = -(tmp1).re + y[n].re;
    ((out)[n]).im = -(tmp1).im + y[n].im;
  }
  return;
}

//alpha*x + beta*y + z
void saxsbypz(hmc_spinor_field * x, hmc_spinor_field * y,  hmc_spinor_field * z, hmc_complex * alpha, hmc_complex * beta, hmc_spinor_field * out){
  for (int n=0; n<SPINORFIELDSIZE; n++) {
    hmc_complex tmp1 = complexmult(alpha,&x[n]);
    hmc_complex tmp2 = complexmult(beta,&y[n]);
    ((out)[n]).re = (tmp1).re + (tmp2).re + z[n].re;
    ((out)[n]).im = (tmp1).im + (tmp2).im + z[n].im;
  }
  return;
}

void saxsbypz_eoprec(hmc_eoprec_spinor_field * x, hmc_eoprec_spinor_field * y,  hmc_eoprec_spinor_field * z, hmc_complex * alpha, hmc_complex * beta, hmc_eoprec_spinor_field * out){
  for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    hmc_complex tmp1 = complexmult(alpha,&x[n]);
    hmc_complex tmp2 = complexmult(beta,&y[n]);
    ((out)[n]).re = (tmp1).re + (tmp2).re + z[n].re;
    ((out)[n]).im = (tmp1).im + (tmp2).im + z[n].im;
  }
  return;
}

void create_point_source(hmc_spinor_field* b, int i, int spacepos, int timepos, hmc_float kappa, hmc_float mu, __global hmc_ocl_gaugefield * gaugefield){
  set_zero_spinor(b);

  int color = spinor_color(i);
  int spin = spinor_spin(i,color);

  b[spinor_field_element(spin,color,spacepos,timepos)].re = sqrt(2.*kappa);

  return;
}

//!!CP: LZ should update this...
void create_point_source_eoprec(hmc_eoprec_spinor_field* be,hmc_eoprec_spinor_field* bo,int i,int spacepos,int timepos,hmc_float kappa, hmc_float mu, hmc_float theta,hmc_float chem_pot_re, hmc_float chem_pot_im, __global hmc_ocl_gaugefield* gaugefield){
  
  
  return;
}

void spinor_apply_bc(hmc_spinor * in, hmc_float theta){
  for(int n = 0; n<SPINORSIZE; n++){
    hmc_float tmp1 = in[n].re;
    hmc_float tmp2 = in[n].im;
    in[n].re = cos(theta)*tmp1 - sin(theta)*tmp2;
    in[n].im = sin(theta)*tmp1 + cos(theta)*tmp2;
  }
  return; 
}

//opencl_random.cl

typedef ulong4 hmc_ocl_ran;
//PRNG as described in NR3, implemented by MB
inline ulong nr3_int64(__global hmc_ocl_ran * state ) {
	(*state).x = (*state).x * 2862933555777941757L + 7046029254386353087L;
	(*state).y ^= (*state).y >> 17; (*state).y ^= (*state).y << 31; (*state).y ^= (*state).y >> 8;
	(*state).z = 4294957665U*((*state).z & 0xffffffff) + ((*state).z >> 32);
	ulong tmp = (*state).x ^ ((*state).x << 21); tmp ^= tmp >> 35; tmp ^= tmp << 4;
	return (tmp + (*state).y) ^ (*state).z;
}
inline float ocl_new_ran(__global hmc_ocl_ran * state ){
	return 5.42101086242752217E-20f * nr3_int64( state );
}
inline uint nr3_int32(__global hmc_ocl_ran * state ){
	return (uint) nr3_int64( state );
}
int random_int( int range, __global hmc_ocl_ran* taus_state ){
	return convert_int( ocl_new_ran( taus_state ) * range );
}
//returns 1,2,3 in a random way
void random_1_2_3 (int rand[3], __global hmc_ocl_ran * rnd) { 
  rand[ 0 ] = random_int( 3, rnd ) + 1;
  do{
    rand[ 1 ] = random_int( 3, rnd ) + 1;
  } while( rand[ 0 ] == rand[ 1 ] );
  rand[ 2 ] = 6 - rand[ 1 ] - rand[ 0 ];
}
// Construct new SU2 matrix using improved alg by Kennedy Pendleton
void SU2Update(__private hmc_float dst [su2_entries], const hmc_float alpha, __global hmc_ocl_ran * rnd){
  hmc_float delta;
  hmc_float a0 ;
  hmc_float eta ; 
  do {
    delta = -log(ocl_new_ran(rnd))/alpha*pow(cos(2. * PI * ocl_new_ran(rnd)), 2.) -log(ocl_new_ran(rnd))/alpha;
    a0 = 1.-delta;
    eta = ocl_new_ran(rnd);
  } while ( (1.-0.5*delta) < eta*eta); 
  hmc_float phi = 2.*PI*ocl_new_ran(rnd);
  hmc_float theta = asin(2.*ocl_new_ran(rnd) - 1.);
  dst[0] = a0;
  dst[1] = sqrt(1.-a0 * a0)*cos(theta) * cos(phi);
  dst[2] = sqrt(1.-a0 * a0)*cos(theta) * sin(phi);
  dst[3] = sqrt(1.-a0 * a0)*sin(theta);
}
//opencl_update_heatbath.cl

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

void inline perform_heatbath(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd, int pos, int t, int id){
  
    hmc_ocl_su3matrix U[SU3SIZE];
    hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
    hmc_ocl_staplematrix staplematrix[STAPLEMATRIXSIZE];
    int order[3]; 
    hmc_complex w [su2_entries];
    hmc_float w_pauli[su2_entries];
    hmc_float k;
    hmc_float r_pauli[su2_entries];
    hmc_float beta_new;

    random_1_2_3(order, &rnd[id]);
    get_su3matrix(U, gaugefield, pos, t, mu);
  
    hmc_complex det = det_su3matrix(U);
    hmc_complex detadj = complexconj(&det);
    hmc_complex detsqnorm = complexmult(&det, &detadj);
    if( (detsqnorm.re - hmc_one_f) <= projectioneps)
      project_su3(U);

    calc_staple(gaugefield, staplematrix, pos, t, mu);

    for(int i=0; i<NC; i++)
    {
      multiply_staplematrix(W, U, staplematrix); 
      reduction(w, W, order[i]);

      w_pauli[0] = 0.5*(w[0].re + w[3].re);
      w_pauli[1] = 0.5*(w[1].im + w[2].im);
      w_pauli[2] = 0.5*(w[1].re - w[2].re);
      w_pauli[3] = 0.5*(w[0].im - w[3].im);
      k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );

      beta_new =  2.*beta / NC*k;
      SU2Update(r_pauli, beta_new, &rnd[id]);
    
      /*
      w[0].re = (r_pauli[0]*w_pauli[0] + r_pauli[1]*w_pauli[1] + r_pauli[2]*w_pauli[2] + r_pauli[3]*w_pauli[3] )/k;
      w[0].im = (w_pauli[0]*r_pauli[3] - w_pauli[3]*r_pauli[0] + r_pauli[1]*w_pauli[2] - r_pauli[2]*w_pauli[1] )/k;
      w[1].re = (w_pauli[0]*r_pauli[2] - w_pauli[2]*r_pauli[0] + r_pauli[3]*w_pauli[1] - r_pauli[1]*w_pauli[3] )/k;
      w[1].im = (w_pauli[0]*r_pauli[1] - w_pauli[1]*r_pauli[0] + r_pauli[2]*w_pauli[3] - r_pauli[3]*w_pauli[2] )/k;
      w[2].re = -(w_pauli[0]*r_pauli[2] - w_pauli[2]*r_pauli[0] + r_pauli[3]*w_pauli[1] - r_pauli[1]*w_pauli[3] )/k;
      w[2].im = (w_pauli[0]*r_pauli[1] - w_pauli[1]*r_pauli[0] + r_pauli[2]*w_pauli[3] - r_pauli[3]*w_pauli[2] )/k;
      w[3].re = (r_pauli[0]*w_pauli[0] + r_pauli[1]*w_pauli[1] + r_pauli[2]*w_pauli[2] + r_pauli[3]*w_pauli[3] )/k;
      w[3].im = -(w_pauli[0]*r_pauli[3] - w_pauli[3]*r_pauli[0] + r_pauli[1]*w_pauli[2] - r_pauli[2]*w_pauli[1] )/k;
      */
    
      //old:
        w_pauli[0] = w_pauli[0]/k;
        w_pauli[1] = -w_pauli[1]/k;
        w_pauli[2] = -w_pauli[2]/k;
        w_pauli[3] = -w_pauli[3]/k;
    	
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
    
      hmc_ocl_su3matrix extW[SU3SIZE]; 
      extend (extW, order[i], w); 
      accumulate_su3matrix_prod(extW, U);
      copy_su3matrix(U, extW);
    }
    put_su3matrix(gaugefield, U, pos, t, mu);
  
  return;
}

__kernel void heatbath_even(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd){
  int t, pos, id, id_tmp, size;
  id_tmp = get_global_id(0);
  size = get_global_size(0);
  for(id = id_tmp; id<VOLSPACE*NTIME/2; id+=size){
    get_even_site(id, &pos, &t);
    perform_heatbath(gaugefield, beta, mu, rnd, pos, t, id_tmp);
  }
  return;
}

__kernel void heatbath_odd(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd){

  int t, pos, id, id_tmp, size;
  id_tmp = get_global_id(0);
  size = get_global_size(0);
  for(id = id_tmp; id<VOLSPACE*NTIME/2; id+=size){
    get_odd_site(id, &pos, &t);
    perform_heatbath(gaugefield, beta, mu, rnd, pos, t, id_tmp);
  }
  return;
}

void inline perform_overrelaxing(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd, int pos, int t, int id){
  
    hmc_ocl_su3matrix U[SU3SIZE];
    hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
    hmc_ocl_staplematrix staplematrix[STAPLEMATRIXSIZE];
  
    hmc_complex w [su2_entries];
    hmc_float w_pauli[su2_entries];
    hmc_float k;
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
    for(int i=0; i<NC; i++)
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

__kernel void overrelax_even(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd){
  int t, pos, id, id_tmp, size;
  id_tmp = get_global_id(0);
  size = get_global_size(0);
  for(id = id_tmp; id<VOLSPACE*NTIME/2; id+=size){
    get_even_site(id, &pos, &t);
    perform_overrelaxing(gaugefield, beta, mu, rnd, pos, t, id_tmp);  
  }  
  return;
}

__kernel void overrelax_odd(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd){
  int t, pos, id, id_tmp, size;
  id_tmp = get_global_id(0);
  size = get_global_size(0);
  for(id = id_tmp; id<VOLSPACE*NTIME/2; id+=size){
    get_odd_site(id, &pos, &t);
    perform_overrelaxing(gaugefield, beta, mu, rnd, pos, t, id_tmp);  
  }
  return;
}

//opencl_solver.cl
void M_diag(hmc_spinor_field* in, hmc_spinor_field* out, hmc_float kappa, hmc_float mu){
  for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
    for(int timepos=0; timepos<NTIME; timepos++) {
      hmc_spinor spinout[SPINORSIZE];
      hmc_spinor spintmp[SPINORSIZE];
      get_spinor_from_field(in,spinout,spacepos,timepos);
      hmc_float twistfactor = 2*kappa*mu;
      multiply_spinor_factor_gamma5(spinout,spintmp,twistfactor);
      spinors_accumulate(spinout,spintmp);
      put_spinor_to_field(spinout,out,spacepos,timepos);
    }
  }
  return; 
}

void dslash(hmc_spinor_field* in, hmc_spinor_field* out,__global hmc_ocl_gaugefield* gaugefield, hmc_float theta){

  hmc_spinor spinout[SPINORSIZE];

  for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
    for(int timepos=0; timepos<NTIME; timepos++) {

      hmc_spinor spinnext[SPINORSIZE];
      hmc_spinor spinprev[SPINORSIZE];
      hmc_spinor tmp[SPINORSIZE];
      hmc_ocl_su3matrix u;
      hmc_ocl_su3matrix udagger;
      int next;
      int prev;
      hmc_float theta = 0.;

      //like in host_geometry
      int coord[NDIM];
      coord[0]=0;
      for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(spacepos,j);
      
      set_local_zero_spinor(spinout);    
   
      // spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
      next = (timepos+1)%NTIME;
      prev = (timepos-1+NTIME)%NTIME;

      get_spinor_from_field(in, spinnext, spacepos, next);
      get_spinor_from_field(in, spinprev, spacepos, prev);

      if(next == 1) spinor_apply_bc(spinnext, theta);
      else if(prev == NTIME) spinor_apply_bc(spinprev, theta);
      
      get_su3matrix(&u,gaugefield,spacepos,timepos,0);
      get_su3matrix(&udagger,gaugefield,spacepos,prev,0);
      adjoin_su3matrix(&udagger);

      multiply_spinor_gamma0(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma0(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&udagger,spinprev,tmp);
      spinors_accumulate(spinout,tmp);

    // spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev
    
      next = get_neighbor(spacepos,1);
      prev = get_lower_neighbor(spacepos,1);
      get_spinor_from_field(in, spinnext, next, timepos);
      get_spinor_from_field(in, spinprev, prev, timepos);
      get_su3matrix(&u,gaugefield,spacepos,timepos,1);
      get_su3matrix(&udagger,gaugefield,prev,timepos,1);
      adjoin_su3matrix(&udagger);

      if(coord[1] == NSPACE) spinor_apply_bc(spinnext, theta);
      else if(coord[1] == 0) spinor_apply_bc(spinprev, theta);
      
      multiply_spinor_gamma1(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma1(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&u,spinprev,tmp);
      spinors_accumulate(spinout,tmp);


    // spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
      next = get_neighbor(spacepos,2);
      prev = get_lower_neighbor(spacepos,2);
      get_spinor_from_field(in, spinnext, next, timepos);
      get_spinor_from_field(in, spinprev, prev, timepos);
      get_su3matrix(&u,gaugefield,spacepos,timepos,2);
      get_su3matrix(&udagger,gaugefield,prev,timepos,2);
      adjoin_su3matrix(&udagger);
      
      if(coord[2] == NSPACE) spinor_apply_bc(spinnext, theta);
      else if(coord[2] == 0) spinor_apply_bc(spinprev, theta);

      multiply_spinor_gamma2(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma2(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&u,spinprev,tmp);
      spinors_accumulate(spinout,tmp);    


    // spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev
      next = get_neighbor(spacepos,3);
      prev = get_lower_neighbor(spacepos,3);
      get_spinor_from_field(in, spinnext, next, timepos);
      get_spinor_from_field(in, spinprev, prev, timepos);
      get_su3matrix(&u,gaugefield,spacepos,timepos,3);
      get_su3matrix(&udagger,gaugefield,prev,timepos,3);
      adjoin_su3matrix(&udagger);

      if(coord[3] == NSPACE) spinor_apply_bc(spinnext, theta);
      else if(coord[3] == 0) spinor_apply_bc(spinprev, theta);
      
      multiply_spinor_gamma3(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma3(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&u,spinprev,tmp);
      spinors_accumulate(spinout,tmp);

      
      put_spinor_to_field(spinout,out,spacepos,timepos);
    }
  }

  return;
}



void M( hmc_spinor_field* in, hmc_spinor_field* out,__global hmc_ocl_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta){
  
  M_diag(in, out, kappa, mu);    
  hmc_spinor_field tmp[SPINORFIELDSIZE];
  dslash(in,tmp,gaugefield, theta);

  hmc_complex kappa_cmplx = {kappa, 0.};
  saxpy(tmp, out, &kappa_cmplx, out);

  return;
}


void bicgstab(__global hmc_spinor_field* inout, hmc_eoprec_spinor_field* source, __global hmc_ocl_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, int cgmax){

  //BiCGStab according to hep-lat/9404013

  hmc_spinor_field rn[SPINORFIELDSIZE];
  hmc_spinor_field rhat [SPINORFIELDSIZE];
  hmc_spinor_field v [SPINORFIELDSIZE];
  hmc_spinor_field p [SPINORFIELDSIZE];
  hmc_spinor_field s [SPINORFIELDSIZE];
  hmc_spinor_field t [SPINORFIELDSIZE];
  hmc_spinor_field inout_tmp [SPINORFIELDSIZE];
  hmc_complex rho;
  hmc_complex rho_next;
  hmc_complex alpha;
  hmc_complex omega;
  hmc_complex beta;

  hmc_complex tmp1;
  hmc_complex tmp2;
  hmc_complex one = hmc_complex_one;
  hmc_complex minusone = hmc_complex_minusone;
 
  hmc_float resid;

  //CP: copy spinor from global to local mem, this must be avoided in the end
  get_spinor(inout, inout_tmp);

  for(int iter=0; iter<cgmax; iter++){

    if(iter%iter_refresh==0) {
      //fresh start
      M(inout_tmp,rn,gaugefield,kappa,mu, theta);
      saxpy(rn, source, &one, rn);
      copy_spinor(rn, rhat);

      alpha = hmc_complex_one;
      omega = hmc_complex_one;
      rho = hmc_complex_one;
      set_zero_spinor(v);
      set_zero_spinor(p);
    }

    rho_next = scalar_product(rhat,rn);
    tmp1 = complexdivide(&rho_next,&rho);
    rho = rho_next;
    tmp2 = complexdivide(&alpha,&omega);
    beta = complexmult(&tmp1,&tmp2);

    tmp1 = complexmult(&beta,&omega);
    tmp2 = complexmult(&minusone,&tmp1);
    saxsbypz(p, v, rn, &beta, &tmp2, p);

    M(p,v,gaugefield,kappa,mu, theta);

    tmp1 = scalar_product(rhat,v);
    alpha = complexdivide(&rho,&tmp1);

    saxpy(v, rn, &alpha, s);

    M(s,t,gaugefield,kappa,mu, theta);

    tmp1 = scalar_product(t,s);
    tmp2 = scalar_product(t,t);
    omega = complexdivide(&(tmp1),&(tmp2));

    saxpy(t, s, &omega, rn);

    saxsbypz(p, s, inout_tmp, &alpha, &omega, inout_tmp);

    resid = global_squarenorm(rn);

    if(resid<epssquare) {
      hmc_spinor_field aux [SPINORFIELDSIZE];
      M(inout_tmp,aux,gaugefield,kappa,mu, theta);
      saxpy(aux, source, &one, aux);
      hmc_float trueresid = global_squarenorm(aux);

      if(trueresid<epssquare) {
        put_spinor(inout_tmp, inout);
        return;
      }
    }
  }
  

for(int i = 0; i<SPINORFIELDSIZE; i++){
if(i==0)  inout_tmp[i].re = resid;
else inout_tmp[i].re = 0.;
  inout_tmp[i].im = 0.;
}


  //CP: copy spinor back to global mem, this musst be avoided in the end
  put_spinor(inout_tmp, inout);
  //put_spinor(rn, inout);

  return;
}

void solver(__global hmc_spinor_field* in,__global hmc_spinor_field* out,__private hmc_spinor_field* b,__global hmc_ocl_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, int cgmax){

  convert_to_kappa_format(in, kappa);
  bicgstab(in, b, gaugefield, kappa, mu, theta, cgmax);
  convert_from_kappa_format(in, out, kappa);

  return;
}


//opencl_gaugeobservables.cl

__kernel void plaquette(__global hmc_ocl_gaugefield * field,__global hmc_float * plaq_out, __global hmc_float* tplaq, __global hmc_float* splaq){

  int t, pos, id, id_tmp, size;
  id_tmp = get_global_id(0);
  size = get_global_size(0);

  hmc_float plaq=0;
  hmc_float splaq_tmp=0;
  hmc_float tplaq_tmp=0;
  hmc_float tmpfloat = 0;

  hmc_ocl_su3matrix tmp[SU3SIZE]; 
  hmc_ocl_su3matrix prod[SU3SIZE];

  for(id = id_tmp; id<VOLSPACE*NTIME/2; id+=size){
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
  }
  
  (plaq_out)[id_tmp] += plaq;
  (splaq)[id_tmp] += splaq_tmp;
  (tplaq)[id_tmp] += tplaq_tmp;
  
  //wait for all threads to end calculations, does this work in a kernel???
  //cl_finish(queue);
  
  //perform reduction
  int cut1;
  int cut2 = size;
  if(size > 128){
    for(cut1 = 128; cut1>0; cut1/=2){
      for(int i = id_tmp+cut1; i < cut2; i+=cut1){
	(plaq_out)[id_tmp] += (plaq_out)[i];
	(splaq)[id_tmp] += (splaq)[i];
	(tplaq)[id_tmp] += (tplaq)[i];
      }
      cut2 = cut1;
    }
  }
  else if(id_tmp == 0) {
    for(int i = id_tmp; i < size; i++){
      (plaq_out)[id_tmp] += (plaq_out)[i];
      (splaq)[id_tmp] += (splaq)[i];
      (tplaq)[id_tmp] += (tplaq)[i];
    }
  }
    
  return;
}

__kernel void polyakov(__global hmc_ocl_gaugefield * field, __global hmc_complex * out){
  
  int t, pos, id, id_tmp, size;
  id_tmp = get_global_id(0);
  size = get_global_size(0);
  int const tdir = 0;
   
  hmc_ocl_su3matrix prod[SU3SIZE];
  hmc_ocl_su3matrix tmp[SU3SIZE];
  unit_su3matrix(prod);
  hmc_complex tmpcomplex;
  hmc_complex tmp_pol;
  tmp_pol.re = 0.;
  tmp_pol.im = 0.;
  
  for(id = id_tmp; id<VOLSPACE/2; id+=size){
    
    //calc polyakov loop at even site
    get_even_site(id, &pos, &t);
    for(int t=0; t<NTIME; t++) {
      get_su3matrix(tmp,field,pos,t,tdir);
      accumulate_su3matrix_prod(prod,tmp);
    }
    tmpcomplex = trace_su3matrix(prod);
    (tmp_pol).re += tmpcomplex.re;
    (tmp_pol).im += tmpcomplex.im;
    
    //calc polyakov loop at odd site
    get_odd_site(id, &pos, &t);
    for(int t=0; t<NTIME; t++) {
      get_su3matrix(tmp,field,pos,t,tdir);
      accumulate_su3matrix_prod(prod,tmp);
    }
    tmpcomplex = trace_su3matrix(prod);
    (tmp_pol).re += tmpcomplex.re;
    (tmp_pol).im += tmpcomplex.im;
  }
  ((out)[id_tmp]).re += tmp_pol.re;
  ((out)[id_tmp]).im += tmp_pol.im;
  
  //wait for all threads to end calculations, does this work in a kernel???
  //cl_finish(queue);
  
  //perform reduction
  int cut1;
  int cut2 = size;
  if(size > 128){
    for(cut1 = 128; cut1>0; cut1/=2){
      for(int i = id_tmp+cut1; i < cut2; i+=cut1){
	((out)[id_tmp]).re +=  ((out)[i]).re;
	((out)[id_tmp]).im +=  ((out)[i]).im;
      }
      cut2 = cut1;
    }
  }
  else if(id_tmp == 0) {
    for(int i = id_tmp; i < size; i++){
     ((out)[id_tmp]).re +=  ((out)[i]).re;
     ((out)[id_tmp]).im +=  ((out)[i]).im;
    }
  }
  
  return;
}
//opencl_fermionobservables.cl

//!!CP: LT should update this...
void simple_correlator(__global hmc_spinor_field * in, __global hmc_spinor_field * spinor_out, __global hmc_ocl_gaugefield* gaugefield, __global hmc_complex * out, hmc_float kappa, hmc_float mu, hmc_float theta, int cgmax){

  //pseudo scalar, flavour multiplet
  for(int z=0; z<NSPACE; z++) {
    out[z].re = 0;
    out[z].im = 0;
  }

  hmc_spinor_field b[SPINORFIELDSIZE];

  for(int k=0; k<NC*NSPIN; k++) {
    create_point_source(b,k,0,0,kappa,mu,gaugefield);
/*
for(int i = 0; i<SPINORFIELDSIZE; i++){
b[i].re= 0;
b[i].im = 0.;
}
*/
    solver(in, spinor_out, b, gaugefield, kappa, mu, theta, 10000);

    for(int timepos = 0; timepos<NTIME; timepos++) {
      for(int spacepos = 0; spacepos<VOLSPACE; spacepos++) {
	for(int alpha = 0; alpha<NSPIN; alpha++) {
	  for(int c = 0; c<NC; c++) {
	    //	    int j = spinor_element(alpha,c);
	    int n = spinor_field_element(alpha, c, spacepos, timepos);
	    int z = get_spacecoord(spacepos, 3);
	    hmc_complex tmp = spinor_out[n];
	    hmc_complex ctmp = complexconj(&tmp);
	    hmc_complex incr = complexmult(&ctmp,&tmp);
	//    if(timepos == 0 && spacepos == 0){
 	    out[z].re += incr.re;
 	    out[z].im += incr.im;
//}
	  }
	}
      }
    }
  }



  return;
}
//opencl_testing.cl

void testing_heatbath_norandommat_no123(hmc_ocl_su3matrix * in, hmc_ocl_staplematrix * staple_in, hmc_ocl_su3matrix * out, hmc_float beta){
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
      
  hmc_ocl_su3matrix su3_in[SU3SIZE];
  copy_su3matrix(su3_in, in);
  
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
//   random_1_2_3(order);
  order[0] = 1; order[1] = 2; order[2] = 3;
  hmc_complex det = det_su3matrix(su3_in);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(su3_in);
  
  for(int i=0; i<3; i++)
  {
    //Produkt aus U und staple, saved in W
    multiply_staplematrix(W, su3_in, staple_in);
    //Reduktion auf SU(2) submatrix
    reduction (w, W, order[i]);
    //Darstellung von W in Paulibasis
    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);

    //Berechnung der Normierung k, entspricht Determinante in normaler Basis
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
    w_pauli[0] = w_pauli[0]/k;
    w_pauli[1] = -w_pauli[1]/k;
    w_pauli[2] = -w_pauli[2]/k;
    w_pauli[3] = -w_pauli[3]/k;
	
    //beta' = 2Beta/Nc*k
    hmc_float beta_neu =  2.*beta / convert_float(NC)*k;
      
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
    hmc_ocl_su3matrix extW[SU3SIZE], tmp[SU3SIZE];
    extend (extW, order[i], w);
	
    multiply_su3matrices(tmp, extW, su3_in);
    copy_su3matrix(su3_in, tmp);
  }
  copy_su3matrix(out, su3_in);
  
  return;
}

void testing_heatbath_no123(hmc_ocl_su3matrix * in, hmc_ocl_staplematrix * staple_in, hmc_ocl_su3matrix * out, hmc_float beta, hmc_float * rnd_array, int * cter_out){
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
      
  hmc_ocl_su3matrix su3_in[SU3SIZE];
  copy_su3matrix(su3_in, in);
  
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
//   random_1_2_3(order);
  order[0] = 1; order[1] = 2; order[2] = 3;
  hmc_complex det = det_su3matrix(su3_in);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(su3_in);
  
  int cter = 0;
  for(int i=0; i<3; i++)
  {
    //Produkt aus U und staple, saved in W
    multiply_staplematrix(W, su3_in, staple_in);
    //Reduktion auf SU(2) submatrix
    reduction (w, W, order[i]);
    //Darstellung von W in Paulibasis
    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);

    //Berechnung der Normierung k, entspricht Determinante in normaler Basis
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
    w_pauli[0] = w_pauli[0]/k;
    w_pauli[1] = -w_pauli[1]/k;
    w_pauli[2] = -w_pauli[2]/k;
    w_pauli[3] = -w_pauli[3]/k;
	
    //beta' = 2Beta/Nc*k
    hmc_float beta_neu =  2.*beta / convert_float(NC)*k;
      
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
    hmc_ocl_su3matrix extW[SU3SIZE], tmp[SU3SIZE];
    extend (extW, order[i], w);
	
    multiply_su3matrices(tmp, extW, su3_in);
    copy_su3matrix(su3_in, tmp);
  }
  copy_su3matrix(out, su3_in);
  (*cter_out) = cter;
  return;
}

void testing_heatbath(hmc_ocl_su3matrix * in, hmc_ocl_staplematrix * staple_in, hmc_ocl_su3matrix * out, hmc_float beta, hmc_float * rnd_array, int * cter_out){
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
      
  hmc_ocl_su3matrix su3_in[SU3SIZE];
  copy_su3matrix(su3_in, in);
  int cter = 0;
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
  
  order[0] = convert_int( rnd_array[cter] * 3 ) + 1;
  cter ++;
  do
    {order[1] = convert_int( rnd_array[cter] * 3 ) + 1;
     cter ++;}
  while (order[1] == order[0]);
  order[2] = 6 - order[1] - order[0];
    
  hmc_complex det = det_su3matrix(su3_in);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(su3_in);
  
  for(int i=0; i<3; i++)
  {
    //Produkt aus U und staple, saved in W
    multiply_staplematrix(W, su3_in, staple_in);
    //Reduktion auf SU(2) submatrix
    reduction (w, W, order[i]);
    //Darstellung von W in Paulibasis
    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);

    //Berechnung der Normierung k, entspricht Determinante in normaler Basis
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
    w_pauli[0] = w_pauli[0]/k;
    w_pauli[1] = -w_pauli[1]/k;
    w_pauli[2] = -w_pauli[2]/k;
    w_pauli[3] = -w_pauli[3]/k;
	
    //beta' = 2Beta/Nc*k
    hmc_float beta_neu =  2.*beta / convert_float(NC)*k;
      
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
    hmc_ocl_su3matrix extW[SU3SIZE], tmp[SU3SIZE];
    extend (extW, order[i], w);
	
    multiply_su3matrices(tmp, extW, su3_in);
    copy_su3matrix(su3_in, tmp);
  }
  copy_su3matrix(out, su3_in);
  (*cter_out) = cter;
  return;
}

__kernel void test(
__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int nsteps,__global hmc_float* check
//random_test-args
, __global hmc_ocl_ran * rnd, __global int * random_field_int, __global float * random_field_float, __global hmc_float * su2mat, const int size_1, const int size_2
//heatbath_test-args
 ,__global hmc_ocl_su3matrix * heatbath_link_in, __global hmc_ocl_su3matrix * heatbath_staple_in, __global hmc_ocl_su3matrix * heatbath_link_out, __global hmc_float * heatbath_rnd_array_in,__global int * heatbath_cter
//solver_test-args
,__global hmc_ocl_gaugefield* gaugefield2, __global hmc_spinor_field * solver_spinor_in, __global hmc_spinor_field * solver_spinor_out, __global hmc_complex * solver_correlator
)
{
  int id = get_global_id(0);
  if(id >0) return;
  else{
  //test by LZ
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


  //CP: random number test
  
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
  
  //CP: heatbath test
  //take a link and its staplematrix after 200 host-iterations and run the host code and the device code on it

  hmc_ocl_su3matrix out_tmp[SU3SIZE];
  hmc_ocl_su3matrix in_tmp[SU3SIZE];
  hmc_ocl_staplematrix staple_tmp[STAPLEMATRIXSIZE];
  int heatbath_rnd_array_size = 10000;
  hmc_float heatbath_rnd_array[10000];
  for(int i = 0; i<heatbath_rnd_array_size; i++){
    heatbath_rnd_array[i] = heatbath_rnd_array_in[i];
  }
  
  for(int i = 0; i<SU3SIZE; i++){
    (in_tmp[i]).re = (heatbath_link_in[i]).re;
    (in_tmp[i]).im = (heatbath_link_in[i]).im;
  }
  for(int i = 0; i<STAPLEMATRIXSIZE; i++){
    (staple_tmp[i]).re = (heatbath_staple_in[i]).re;
    (staple_tmp[i]).im = (heatbath_staple_in[i]).im;
  }
  
  int cter;
  //CP: three different possibilities to test the heatbath
//   testing_heatbath_norandommat_no123(in_tmp, staple_tmp, out_tmp, 4.2);
//   testing_heatbath_no123(in_tmp, staple_tmp, out_tmp, 4.2, heatbath_rnd_array, &cter);
  testing_heatbath(in_tmp, staple_tmp, out_tmp, 4.2, heatbath_rnd_array, &cter);
  
  for(int i = 0; i<SU3SIZE; i++){
    (heatbath_link_out[i]).re = (out_tmp[i]).re;
    (heatbath_link_out[i]).im = (out_tmp[i]).im;
  }
  heatbath_cter[0] = cter;

  //CP: solver test: invert a small matrix on a cold-gaugeconfiguration and calculate the pion propagator

  simple_correlator(solver_spinor_in, solver_spinor_out, gaugefield2, solver_correlator, 0.125, 0.06, 0., 1000);

  return;
  } //else
}

//EOF
