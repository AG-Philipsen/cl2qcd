#include "host_operations_su3matrix.h"


//operations that contain explicit SU(3) indices!!!

hmc_complex det_su3matrix(hmc_su3matrix * U){
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
  tmp1 = complexmult( &(*U)[0], &(*U)[3] );
  tmp2 = reconstruct_su3(U,2);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&det,&tmp3);

  tmp1 = complexmult( &(*U)[2], &(*U)[5] );
  tmp2 = reconstruct_su3(U,0);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&det,&tmp3);

  tmp1 = complexmult( &(*U)[4], &(*U)[1] );
  tmp2 = reconstruct_su3(U,1);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&det,&tmp3);

  tmp1 = complexmult( &(*U)[3], &(*U)[4] );
  tmp2 = reconstruct_su3(U,0);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&subdet,&tmp3);

  tmp1 = complexmult( &(*U)[5], &(*U)[0] );
  tmp2 = reconstruct_su3(U,1);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&subdet,&tmp3);

  tmp1 = complexmult( &(*U)[1], &(*U)[2] );
  tmp2 = reconstruct_su3(U,2);
  tmp3 = complexmult( &tmp1, &tmp2 );
  complexaccumulate(&subdet,&tmp3);

  det.re -= subdet.re;
  det.im -= subdet.im;

#else
  hmc_complex det, det1, det2, det3, det4, det5, det6, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
  det.re=0;
  det.im=0;
  tmp1 = complexmult( &(*U)[1][1], &(*U)[2][2] );
  det1 = complexmult( &(*U)[0][0] , &tmp1);
  tmp2 = complexmult( &(*U)[1][2], &(*U)[2][0] );
  det2 = complexmult( &(*U)[0][1] , &tmp2);
  tmp3 = complexmult( &(*U)[1][0], &(*U)[2][1] );
  det3 = complexmult( &(*U)[0][2] , &tmp3);
  tmp4 = complexmult( &(*U)[1][1], &(*U)[2][0] );
  det4 = complexmult( &(*U)[0][2] , &tmp4);
  tmp5 = complexmult( &(*U)[1][0], &(*U)[2][2] );
  det5 = complexmult( &(*U)[0][1] , &tmp5);
  tmp6 = complexmult( &(*U)[1][2], &(*U)[2][1] );
  det6 = complexmult( &(*U)[0][0] , &tmp6);

  det.re = det1.re + det2.re + det3.re - det4.re - det5.re - det6.re;
  det.im = det1.im + det2.im + det3.im - det4.im - det5.im - det6.im;

#endif
  return det;
}

/** @todo memcpy ... */
void copy_su3matrix(hmc_su3matrix *out, hmc_su3matrix *in){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
    (*out)[n] = (*in)[n];
  }
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*out)[a][b] = (*in)[a][b];
    }
  }
#endif
  return;
}

/** @todo memcpy ... */
void copy_staplematrix(hmc_staplematrix *out, hmc_staplematrix *in){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*NC; n++) {
    (*out)[n] = (*in)[n];
  }
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*out)[a][b] = (*in)[a][b];
    }
  }
#endif
  return;
}

/** @todo memset ... */
void zero_su3matrix(hmc_su3matrix * u){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
    (*u)[n].re = 0;
    (*u)[n].im = 0;
  }
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*u)[a][b].re = 0;
      (*u)[a][b].im = 0;
    }
  }
#endif
  return;
}

/** @todo memset ... */
void zero_staplematrix(hmc_staplematrix * u){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*NC; n++) {
    (*u)[n].re = 0;
    (*u)[n].im = 0;
  }
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*u)[a][b].re = 0;
      (*u)[a][b].im = 0;
    }
  }
#endif
  return;
}

Matrixsu3 unit_matrixsu3()
{
	Matrixsu3 out;
	out.e00.re = 1.;
	out.e00.im = 0.;
	out.e01.re = 0.;
	out.e01.im = 0.;
	out.e02.re = 0.;
	out.e02.im = 0.;

	out.e10.re = 0.;
	out.e10.im = 0.;
	out.e11.re = 1.;
	out.e11.im = 0.;
	out.e12.re = 0.;
	out.e12.im = 0.;

#ifndef _RECONSTRUCT_TWELVE_
	out.e20.re = 0.;
	out.e20.im = 0.;
	out.e21.re = 0.;
	out.e21.im = 0.;
	out.e22.re = 1.;
	out.e22.im = 0.;
#endif
	return out;
}

void unit_su3matrix(hmc_su3matrix * u){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
    if( n%NC == 0) {
      (*u)[n].re = hmc_one_f;
    } else {
      (*u)[n].re = 0;
    }
    (*u)[n].im = 0;
  }
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      if(a!=b) {
	(*u)[a][b].re = 0;
      } else {
	(*u)[a][b].re = hmc_one_f;
      }
      (*u)[a][b].im = 0;
    }
  }
#endif
  return;
}

void random_su3matrix(hmc_su3matrix *){
  throw Print_Error_Message("random su3matrix needs to be implemented...");
  return;
}

#ifdef _RECONSTRUCT_TWELVE_
hmc_complex reconstruct_su3(hmc_su3matrix *in, int ncomp){
  int jplusone = (ncomp+1)%NC;
  int jplustwo = (ncomp+2)%NC;
  hmc_complex first = complexmult(&((*in)[(NC-1)*jplusone]),&((*in)[1+(NC-1)*jplustwo]));
  hmc_complex second = complexmult(&((*in)[(NC-1)*jplustwo]),&((*in)[1+(NC-1)*jplusone]));
  hmc_complex result = complexsubtract(&first,&second);
  return complexconj(&result);
}
#endif


void multiply_su3matrices(hmc_su3matrix *out, hmc_su3matrix *p, hmc_su3matrix *q){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
      (*out)[n].re=0;
      (*out)[n].im=0;
      for(int j=0;j<NC;j++) {
	int k = (int)(n/(NC-1));
	int i = n - (NC-1)*k;
	int np = i + (NC-1)*j;
	hmc_complex qcomponent;
	if(j==2) {
	  qcomponent = reconstruct_su3(q,k);
	} else {
	  int nq = j + (NC-1)*k;
	  qcomponent = (*q)[nq];
	}
	hmc_complex tmp = complexmult(&(*p)[np],&qcomponent);
	complexaccumulate(&(*out)[n],&tmp);
      }
    }
#else
  for(int i=0; i<NC; i++) {
    for(int k=0; k<NC; k++) {
      (*out)[i][k].re=0;
      (*out)[i][k].im=0;
      for(int j=0;j<NC;j++) {
	hmc_complex tmp = complexmult(&(*p)[i][j],&(*q)[j][k]);
	complexaccumulate(&(*out)[i][k],&tmp);
      }
    }
  }
#endif
  return;
}

void multiply_staplematrix(hmc_staplematrix *out, hmc_su3matrix *p, hmc_staplematrix *q){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
      (*out)[n].re=0;
      (*out)[n].im=0;
      for(int j=0;j<NC;j++) {
	int k = (int)(n/(NC-1));
	int i = n - (NC-1)*k;
	int np = i + (NC-1)*j;
	hmc_complex qcomponent;
	if(j==2) {
// 	  qcomponent = reconstruct_su3(q,k);
          qcomponent = (*q)[NC*(NC-1)+k];
	} else {
	  int nq = j + (NC-1)*k;
	  qcomponent = (*q)[nq];
	}
	hmc_complex tmp = complexmult(&(*p)[np],&qcomponent);
	complexaccumulate(&(*out)[n],&tmp);
      }
    }
    //the left components:
    hmc_complex X = reconstruct_su3(p,0);
    hmc_complex Y = reconstruct_su3(p,1);
    hmc_complex Z = reconstruct_su3(p,2);
    hmc_complex tmp;
    ((*out)[6]).re=0;
    ((*out)[6]).im=0;
    ((*out)[7]).re=0;
    ((*out)[7]).im=0;
    ((*out)[8]).re=0;
    ((*out)[8]).im=0;
    
    tmp = complexmult(&X,&(*q)[0]);
    complexaccumulate(&(*out)[6],&tmp);
    tmp = complexmult(&Y,&(*q)[1]);
    complexaccumulate(&(*out)[6],&tmp);
    tmp = complexmult(&Z,&(*q)[6]);
    complexaccumulate(&(*out)[6],&tmp);

    tmp = complexmult(&X,&(*q)[2]);
    complexaccumulate(&(*out)[7],&tmp);
    tmp = complexmult(&Y,&(*q)[3]);
    complexaccumulate(&(*out)[7],&tmp);
    tmp = complexmult(&Z,&(*q)[7]);
    complexaccumulate(&(*out)[7],&tmp);

    tmp = complexmult(&X,&(*q)[4]);
    complexaccumulate(&(*out)[8],&tmp);
    tmp = complexmult(&Y,&(*q)[5]);
    complexaccumulate(&(*out)[8],&tmp);
    tmp = complexmult(&Z,&(*q)[8]);
    complexaccumulate(&(*out)[8],&tmp);
    
#else
  for(int i=0; i<NC; i++) {
    for(int k=0; k<NC; k++) {
      (*out)[i][k].re=0;
      (*out)[i][k].im=0;
      for(int j=0;j<NC;j++) {
	hmc_complex tmp = complexmult(&(*p)[i][j],&(*q)[j][k]);
	complexaccumulate(&(*out)[i][k],&tmp);
      }
    }
  }
#endif
  return;
}

void accumulate_su3matrices_add(hmc_staplematrix *p, hmc_su3matrix *q){
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
    complexaccumulate(&(*p)[n], &(*q)[n]);
  }
  for(int n=NC*(NC-1);  n<NC*NC; n++) {
    hmc_complex tmp = reconstruct_su3(q, n-NC*(NC-1)); 
    complexaccumulate(&(*p)[n], &(tmp));
  }  
#else

  for(int i=0; i<NC; i++) {
    for(int k=0; k<NC; k++) {
      complexaccumulate(&(*p)[i][k],&(*q)[i][k]);
    }
  }
#endif
  return;
}

void accumulate_su3matrix_prod(hmc_su3matrix *acc, hmc_su3matrix *mult){
  hmc_su3matrix tmp;
  multiply_su3matrices(&tmp,acc,mult);
  copy_su3matrix(acc,&tmp);
  return;
}

void adjoin_su3matrix(hmc_su3matrix * mat){
#ifdef _RECONSTRUCT_TWELVE_
  hmc_su3matrix tmp;
  copy_su3matrix(&tmp, mat);
  for(int n=0; n<NC*(NC-1); n++) {
    int j = (int)(n/(NC-1));
    int i = n - (NC-1)*j;
    hmc_complex element;
    if ( j==2 ) {
      element = reconstruct_su3(&tmp,i);
    } else {
      int nnew = j + (NC-1)*i;
      element = tmp[nnew];
    }
    (*mat)[n] = complexconj(&element);
  }
#else
  hmc_su3matrix tmp;
  copy_su3matrix(&tmp, mat);
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*mat)[a][b] = complexconj(&(tmp[b][a]));
    }
  }
#endif
  return;
}

hmc_complex trace_su3matrix(hmc_su3matrix * mat){
  hmc_complex trace;
#ifdef _RECONSTRUCT_TWELVE_
  trace = reconstruct_su3(mat,NC-1);
  for(int n=0; n<(NC-1); n++) complexaccumulate(&trace,&((*mat)[NC*n]));
#else
  trace.re=0;
  trace.im=0;
  for(int a=0; a<NC; a++) complexaccumulate(&trace,&((*mat)[a][a]));;
#endif
  return trace;
}

void gaugefield_apply_bc(hmc_su3matrix * in, hmc_float theta){
  hmc_float tmp1,tmp2;
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
    tmp1 = ((*in)[n]).re;
    tmp2 = ((*in)[n]).im;
    ((*in)[n]).re = cos(theta)*tmp1 - sin(theta)*tmp2;
    ((*in)[n]).im = sin(theta)*tmp1 + cos(theta)*tmp2;
  }
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      tmp1 = ((*in)[a][b]).re;
      tmp2 = ((*in)[a][b]).im;
      ((*in)[a][b]).re = cos(theta)*tmp1 - sin(theta)*tmp2;
      ((*in)[a][b]).im = sin(theta)*tmp1 + cos(theta)*tmp2;
    }
  }
#endif
  return;
}

// replace link in with e^(mu.re)*(cos(mu.im) + i*sin(mu.im))
void gaugefield_apply_chem_pot(hmc_su3matrix * u, hmc_su3matrix * udagger, hmc_float chem_pot_re, hmc_float chem_pot_im){
  hmc_float tmp1,tmp2;
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) {
    tmp1 = ((*u)[n]).re;
    tmp2 = ((*u)[n]).im;
    ((*u)[n]).re = exp(chem_pot_re)*( cos(chem_pot_im)*tmp1 - sin(chem_pot_im)*tmp2 );
    ((*u)[n]).im = exp(chem_pot_re)*( sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
    tmp1 = ((*udagger)[n]).re;
    tmp2 = ((*udagger)[n]).im;
    ((*udagger)[n]).re = exp(-chem_pot_re)*( cos(chem_pot_im)*tmp1 + sin(chem_pot_im)*tmp2 );
    ((*udagger)[n]).im = exp(-chem_pot_re)*( -sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
  }
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      tmp1 = ((*u)[a][b]).re;
      tmp2 = ((*u)[a][b]).im;
      ((*u)[a][b]).re = exp(chem_pot_re)*( cos(chem_pot_im)*tmp1 - sin(chem_pot_im)*tmp2 );
      ((*u)[a][b]).im = exp(chem_pot_re)*( sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
      tmp1 = ((*udagger)[a][b]).re;
      tmp2 = ((*udagger)[a][b]).im;
      ((*udagger)[a][b]).re = exp(-chem_pot_re)*( cos(chem_pot_im)*tmp1 + sin(chem_pot_im)*tmp2 );
      ((*udagger)[a][b]).im = exp(-chem_pot_re)*( -sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
    }
  }
#endif
  return;
}

