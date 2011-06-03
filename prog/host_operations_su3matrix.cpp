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

using namespace std;

hmc_error project_su3(hmc_su3matrix *U){

// 	cout << "0" << endl;
	
  //Extract initial vectors
  hmc_complex a[NC];
  hmc_complex b[NC];
  hmc_complex c[NC];
#ifdef _RECONSTRUCT_TWELVE_
  a[0] = (*U)[0];
  a[1] = (*U)[2];
  a[2] = (*U)[4];
  b[0] = (*U)[1];
  b[1] = (*U)[3];
  b[2] = (*U)[5];
  c[0] = reconstruct_su3(U,0);
  c[1] = reconstruct_su3(U,1);
  c[2] = reconstruct_su3(U,2);
#else
    for (int i = 0; i<NC; i++){
     a[i] = (*U)[0][i];
     b[i] = (*U)[1][i];
     c[i] = (*U)[2][i];
    }
#endif
//   	cout << "1" << endl;
  //New SU3-Matrix
  //first vector
  //norm
  hmc_float norm = 0.;
  for (int i=0; i<NC; i++){
    hmc_complex tmp = complexconj(&(a[i]));
    tmp = complexmult (& a[i], & tmp);
    norm += tmp.re;
  }
  norm = 1./sqrt(norm);
  //rescale
  for (int i=0; i<NC; i++){
    //perhaps define a new complex-function for multiplying with a real number
    a[i].re *= norm;
    a[i].im *= norm;
  }
//   	cout << "2" << endl;
  //second vector
  //orthogonal vector
  hmc_complex factor;
  factor.re = 0.0;
  factor.im = 0.0;
  for (int i=0; i<NC; i++){
    hmc_complex tmp;
    tmp = complexconj (&(b[i]));
    tmp = complexmult (&(a[i]), &tmp);
    factor = complexadd (&factor, &tmp);
  }
  for (int i=0; i<NC; i++){
    hmc_complex tmp;
    tmp = complexmult(&factor, &(a[i]));
    b[i] = complexsubtract(&(b[i]), &tmp); 
  }
//   	cout << "3" << endl;
//norm
  norm = 0.;
  for (int i=0; i<NC; i++)
  {
    hmc_complex tmp;
    tmp = complexconj(&(b[i]));
    tmp = complexmult (&(b[i]), &tmp);
    norm +=  tmp.re;
  }
  norm = 1./sqrt(norm);
  //rescale
  for  (int i=0; i<NC; i++){
    b[i].re *= norm;
    b[i].im *= norm;
  }
// 	cout << "4" << endl;
#ifdef _RECONSTRUCT_TWELVE_
  //third vector 
  //orthogonal vector
  hmc_complex tmp;
  hmc_complex tmp2;
  tmp = complexmult(&(a[1]), &(b[2]));
  tmp = complexconj(&tmp);
  tmp2 = complexmult(&(a[2]), &(b[1]));
  tmp2 = complexconj(&tmp2);
  c[0] = complexsubtract(&tmp, &tmp2);
  tmp = complexmult(&(a[2]), &(b[0]));
  tmp = complexconj(&tmp);
  tmp2 = complexmult(&(a[0]), &(b[2]));
  tmp2 = complexconj(&tmp2);
  c[1] = complexsubtract(&tmp, &tmp2);
  
  //Set new values to matrix
  (*U)[0] = a[0];
  (*U)[1] = b[0];
  (*U)[2] = a[1];
  (*U)[3] = b[1];
  (*U)[4] = a[2];
  (*U)[5] = b[2];
#else
  //third vector 
  //orthogonal vector
  hmc_complex tmp;
  hmc_complex tmp2;
  tmp = complexmult(&(a[1]), &(b[2]));
  tmp = complexconj(&tmp);
  tmp2 = complexmult(&(a[2]), &(b[1]));
  tmp2 = complexconj(&tmp2);
  c[0] = complexsubtract(&tmp, &tmp2);
  tmp = complexmult(&(a[2]), &(b[0]));
  tmp = complexconj(&tmp);
  tmp2 = complexmult(&(a[0]), &(b[2]));
  tmp2 = complexconj(&tmp2);
  c[1] = complexsubtract(&tmp, &tmp2);
  tmp = complexmult(&(a[0]), &(b[1]));
  tmp = complexconj(&tmp);
  tmp2 = complexmult(&(a[1]), &(b[0]));
  tmp2 = complexconj(&tmp2);
  c[2] = complexsubtract(&tmp, &tmp2);
  
  //Set new values to matrix
  for(int i=0; i<NC; i++) {
     (*U)[0][i] = a[i];
     (*U)[1][i] = b[i];
     (*U)[2][i] = c[i];
  }
#endif
// 	cout << "5" << endl;
  return HMC_SUCCESS;
}


hmc_error project_su3_old(hmc_su3matrix *U){
//old code
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
    hmc_complex tmp = (*U)[n];
    (*U)[n].re = (tmp.re*norm.re+tmp.im*norm.im)/normsqunorm;
    (*U)[n].im = (tmp.im*norm.re-tmp.re*norm.im)/normsqunorm;
  }
  #else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      hmc_complex tmp = (*U)[a][b]; 
      (*U)[a][b].re = (tmp.re*norm.re+tmp.im*norm.im)/normsqunorm;
      (*U)[a][b].im = (tmp.im*norm.re-tmp.re*norm.im)/normsqunorm;
    }
  }
  #endif
  return HMC_SUCCESS;
}

/** @todo memcpy ... */
hmc_error copy_su3matrix(hmc_su3matrix *out, hmc_su3matrix *in){
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
  return HMC_SUCCESS;
}

/** @todo memcpy ... */
hmc_error copy_staplematrix(hmc_staplematrix *out, hmc_staplematrix *in){
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
  return HMC_SUCCESS;
}

/** @todo memset ... */
hmc_error zero_su3matrix(hmc_su3matrix * u){
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
  return HMC_SUCCESS;
}

/** @todo memset ... */
hmc_error zero_staplematrix(hmc_staplematrix * u){
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
  return HMC_SUCCESS;
}

hmc_error unit_su3matrix(hmc_su3matrix * u){
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
  return HMC_SUCCESS;
}

hmc_error random_su3matrix(hmc_su3matrix * u){
  printf("random su3matrix needs to be implemented...\n");
  exit(HMC_UNDEFINEDERROR);
  return HMC_SUCCESS;
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


hmc_error multiply_su3matrices(hmc_su3matrix *out, hmc_su3matrix *p, hmc_su3matrix *q){
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
  return HMC_SUCCESS;
}

// wird wohl doch nicht gebraucht...
// hmc_error add_su3matrices(hmc_su3matrix *p, hmc_su3matrix *q){
// #ifdef _RECONSTRUCT_TWELVE
//   for(int n=0; n<NC*(NC-1); n++) {
//       (*p)[n].re = (*p)[n].re + (*q)[n].re;
//       (*p)[n].im = (*p)[n].im + (*q)[n].im;
// #else
//   for(int i=0; i<NC; i++) {
//     for(int k=0; k<NC; k++) {
//       (*p)[i][k].re = (*p)[i][k].re + (*q)[i][k].re;
//       (*p)[i][k].im = (*p)[i][k].im + (*q)[i][k].im;
//     }
//   }
// #endif
// }


hmc_error multiply_staplematrix(hmc_staplematrix *out, hmc_su3matrix *p, hmc_staplematrix *q){
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
  return HMC_SUCCESS;
}

hmc_error accumulate_su3matrices_add(hmc_staplematrix *p, hmc_su3matrix *q){
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
  return HMC_SUCCESS;
}

hmc_error accumulate_su3matrix_prod(hmc_su3matrix *acc, hmc_su3matrix *mult){
  hmc_su3matrix tmp;
  multiply_su3matrices(&tmp,acc,mult);
  copy_su3matrix(acc,&tmp);
  return HMC_SUCCESS;
}

hmc_error adjoin_su3matrix(hmc_su3matrix * mat){
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
  return HMC_SUCCESS;
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

void reduction (hmc_complex dest[su2_entries], hmc_staplematrix src, const int rand){
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
  else std::cout<<"error at reduction, rand not 1,2,3"<<std::endl; 
#else
  if(rand == 1)
  {
    dest[0] = src[0][0];
    dest[1] = src[0][1];
    dest[2] = src[1][0];
    dest[3] = src[1][1];
  }
  else if (rand==2)
  {
    dest[0] = src[1][1];
    dest[1] = src[1][2];
    dest[2] = src[2][1];
    dest[3] = src[2][2];
  }
  else if (rand==3)
  {
    dest[0] = src[0][0];
    dest[1] = src[0][2];
    dest[2] = src[2][0];
    dest[3] = src[2][2];
  }
  else
    std::cout<<"error at reduction, rand not 1,2,3"<<std::endl;
#endif
}

// return an SU2 matrix (std basis) extended to SU3 (std basis)
void extend (hmc_su3matrix * dest, const int random, hmc_complex src[su2_entries]){
#ifdef _RECONSTRUCT_TWELVE_
  if (random == 1){
    (*dest)[0] = src[0];
    (*dest)[2] = src[1];
    (*dest)[4] = hmc_complex_zero;
    (*dest)[1] = src[2];
    (*dest)[3] = src[3];
    (*dest)[5] = hmc_complex_zero;
  }
  else if (random == 2){
    (*dest)[0] = hmc_complex_one;
    (*dest)[2] = hmc_complex_zero;
    (*dest)[4] = hmc_complex_zero;
    (*dest)[1] = hmc_complex_zero;
    (*dest)[3] = src[0];
    (*dest)[5] = src[1];
  }
  else if (random == 3){
    (*dest)[0] = src[0];
    (*dest)[2] = hmc_complex_zero;
    (*dest)[4] = src[1];
    (*dest)[1] = hmc_complex_zero;
    (*dest)[3] = hmc_complex_one;
    (*dest)[5] = hmc_complex_zero;
  }
  else
    std::cout<<"error at extend, random not 1,2,3"<<std::endl;

#else
  if (random == 1){
    (*dest)[0][0] = src[0];
    (*dest)[0][1] = src[1];
    (*dest)[0][2] = hmc_complex_zero;
    (*dest)[1][0] = src[2];
    (*dest)[1][1] = src[3];
    (*dest)[1][2] = hmc_complex_zero;
    (*dest)[2][0] = hmc_complex_zero;
    (*dest)[2][1] = hmc_complex_zero;
    (*dest)[2][2] = hmc_complex_one;
  }
  else if (random == 2){
    (*dest)[0][0] = hmc_complex_one;
    (*dest)[0][1] = hmc_complex_zero;
    (*dest)[0][2] = hmc_complex_zero;
    (*dest)[1][0] = hmc_complex_zero;
    (*dest)[1][1] = src[0];
    (*dest)[1][2] = src[1];
    (*dest)[2][0] = hmc_complex_zero;
    (*dest)[2][1] = src[2];
    (*dest)[2][2] = src[3];
  }
  else if (random == 3){
    (*dest)[0][0] = src[0];
    (*dest)[0][1] = hmc_complex_zero;
    (*dest)[0][2] = src[1];
    (*dest)[1][0] = hmc_complex_zero;
    (*dest)[1][1] = hmc_complex_one;
    (*dest)[1][2] = hmc_complex_zero;
    (*dest)[2][0] = src[2];
    (*dest)[2][1] = hmc_complex_zero;
    (*dest)[2][2] = src[3];
  }
  else
    std::cout<<"error at extend, random not 1,2,3"<<std::endl;
#endif
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

//CP: tested version that recreates the matrices made by tmlqcd
/** @todo recheck the factor 0.5 (or F_1_2) that has been deleted here */
/** @todo use of hmc_algebraelement!!!! */
hmc_error build_su3matrix_by_exponentiation(hmc_algebraelement2 inn, hmc_su3matrix* out, hmc_float epsilon){
		//CP: workaround for struct hmc_algebraelement
		hmc_float in [8];
		in[0] = inn.e0;
		in[1] = inn.e1;
		in[2] = inn.e2;
		in[3] = inn.e3;
		in[4] = inn.e4;
		in[5] = inn.e5;
		in[6] = inn.e6;
		in[7] = inn.e7;
		
		// this is the case where one actually evaluates -as matrices- many orders of exp(i*e*P)=1+i*e*P+(1/2)(i*e*P)^2 + ...
		// 1. create mat = (i epsilon)/2 p_i lambda_i    with lambda=gellmann matrix
		hmc_3x3matrix eMat;
		hmc_float halfeps = epsilon;//*F_1_2;
		eMat[0][0].re = 0.0;
		eMat[0][0].im = halfeps*(in[7]*F_1_S3+in[2]);
		eMat[0][1].re = halfeps*in[1];
		eMat[0][1].im = halfeps*in[0];
		eMat[0][2].re = halfeps*in[4];
		eMat[0][2].im = halfeps*in[3];
		eMat[1][0].re = -halfeps*in[1];;
		eMat[1][0].im = halfeps*in[0];;
		eMat[1][1].re = 0.0;
		eMat[1][1].im = halfeps*(in[7]*F_1_S3-in[2]);
		eMat[1][2].re = halfeps*in[6];
		eMat[1][2].im = halfeps*in[5];
		eMat[2][0].re = -halfeps*in[4];
		eMat[2][0].im = halfeps*in[3];
		eMat[2][1].re = -halfeps*in[6];
		eMat[2][1].im = halfeps*in[5];
		eMat[2][2].re = 0.0;
		eMat[2][2].im = -halfeps*2.*in[7]*F_1_S3;
		// 2. start with the exp(...) expansion by using standard 3x3 sum and multiplication
		hmc_3x3matrix eRes, eCurPower, eNextPower, eLastResult;
		hmc_float eAccuracyCheck;
		set_to_3x3_identity(&eRes);
		set_to_3x3_identity(&eCurPower);
// 		hmc_float eCoefficient = 1.0;
		int factorial = 1;
		for(int power=1;power<_EXACT_EXPONENTIATION_MAX_POWER_;power++){
			multiply_3x3matrix(&eNextPower, &eMat, &eCurPower);
			copy_3x3_matrix(&eCurPower, &eNextPower);
			//CP: the calc of factorial can be simplified by just doing factorial*=power!!
			for(int i = power; i>1; i--) factorial *= i;
			multiply_3x3matrix_by_real(&eCurPower, (1./(1.*factorial)));
			factorial = 1;
			copy_3x3_matrix(&eLastResult, &eRes);
			add_3x3matrix(&eRes, &eRes, &eCurPower);
			absoluteDifference_3x3_matrix(&eAccuracyCheck, &eRes, &eLastResult);
			if(eAccuracyCheck < _EXACT_EXPONENTIATION_ACCURACY_){
				break;
			}
		}
		// 3. here I have the exponentiated matrix in 3x3 generic form (eRes), project it
		#ifdef _RECONSTRUCT_TWELVE_
			(*out)[0] = eRes[0][0];
			(*out)[1] = eRes[0][1];
			(*out)[2] = eRes[0][2];
			(*out)[3] = eRes[1][0];
			(*out)[4] = eRes[1][1];
			(*out)[5] = eRes[1][2];
		#else
			(*out)[0][0] = eRes[0][0];
			(*out)[0][1] = eRes[0][1];
			(*out)[0][2] = eRes[0][2];
			(*out)[1][0] = eRes[1][0];
			(*out)[1][1] = eRes[1][1];
			(*out)[1][2] = eRes[1][2];
			(*out)[2][0] = eRes[2][0];
			(*out)[2][1] = eRes[2][1];
			(*out)[2][2] = eRes[2][2];
		#endif // _RECONSTRUCT_TWELVE_
		project_su3(out);
	
	return HMC_SUCCESS;
}


/*
hmc_error build_su3matrix_by_exponentiation(hmc_algebraelement in, hmc_su3matrix* out, hmc_float epsilon){
	// SL: this takes 8 real numbers and builds the su3 matrix exp(i*epsilon*p_i*T_i)
	//     either by truncated "smart" series expansion to order eps^2 or eps^3, or by direct evaluation of the series
	// NOTE: the code here is strictly SU(3)-specific! (and mostly generated by other code...)
	
	#ifdef _USE_MORNGINGSTAR_PEARDON_
	//TODO CP: this has to be implemented
	
	#else

	#ifndef _EXPONENTIATE_ALGEBRA_ALL_ORDERS_

	// Phase 1: calculate (number) coefficients for the reconstruction
	hmc_float beta_0, gamma_0;
	hmc_float beta[8], gamma[8];
	// the above are: coefficients for identity (beta_0+i*gamma_0), and coefficients for the T_l, namely T_L*(beta_l + i*gamma_l).
	#ifdef _EXPONENTIATE_ALGEBRA_ORDER_2_
	
		hmc_float pR, pQ[8], pPtilde[8];

		pR = in[0]*in[0]+in[1]*in[1]+in[2]*in[2]+in[3]*in[3]+in[4]*in[4]+in[5]*in[5]+in[6]*in[6]+in[7]*in[7];

		pPtilde[0] = (F_1_S3*in[7]*in[0])+(F_1_2*(in[5]*in[3]+in[6]*in[4]));
		pPtilde[1] = (F_1_S3*in[7]*in[1])+(F_1_2*(in[5]*in[4]-in[6]*in[3]));
		pPtilde[2] = (F_1_S3*in[7]*in[2]);
		pPtilde[3] = (F_1_2*(in[5]*in[0]-in[6]*in[1]+in[3]*in[2]))-(F_1_2S3*in[7]*in[3]);
		pPtilde[4] = (F_1_2*(in[6]*in[0]+in[5]*in[1]+in[4]*in[2]))-(F_1_2S3*in[7]*in[4]);
		pPtilde[5] = (F_1_2*(in[3]*in[0]+in[4]*in[1]-in[5]*in[2]))-(F_1_2S3*in[7]*in[5]);
		pPtilde[6] = (F_1_2*(in[4]*in[0]-in[3]*in[1]-in[6]*in[2]))-(F_1_2S3*in[7]*in[6]);
		pPtilde[7] = 0.0;

		pQ[2] = (F_1_2*(in[3]*in[3]+in[4]*in[4]-in[5]*in[5]-in[6]*in[6]));
		pQ[7] = (F_1_S3*(in[0]*in[0]+in[1]*in[1]+in[2]*in[2]-in[7]*in[7]))-(F_1_2S3*(in[3]*in[3]+in[4]*in[4]+in[5]*in[5]+in[6]*in[6]));
		pQ[0] = pQ[1] = pQ[3] = pQ[4] = pQ[5] = pQ[6] = 0.0;

		hmc_float eps_squared = epsilon*epsilon;
		beta_0  = 1.0 - (eps_squared*pR/12.);
		gamma_0 = 0.0;
		for(int il=0;il<8;il++){
			beta[il] = -0.25*eps_squared*(pQ[il]+2.0*pPtilde[il]);
			gamma[il]= epsilon*in[il];
		}
	#endif // EXPONENTIATE_ALGEBRA_ORDER_2
	#ifdef _EXPONENTIATE_ALGEBRA_ORDER_3_
		hmc_float pR, pQ[8], pPtilde[8];
		hmc_float pT, pS[8], pU[8];

		pR = in[0]*in[0]+in[1]*in[1]+in[2]*in[2]+in[3]*in[3]+in[4]*in[4]+in[5]*in[5]+in[6]*in[6]+in[7]*in[7];
		
		pPtilde[0] = (F_1_S3*in[7]*in[0])+(F_1_2*(in[5]*in[3]+in[6]*in[4]));
		pPtilde[1] = (F_1_S3*in[7]*in[1])+(F_1_2*(in[5]*in[4]-in[6]*in[3]));
		pPtilde[2] = (F_1_S3*in[7]*in[2]);
		pPtilde[3] = (F_1_2*(in[5]*in[0]-in[6]*in[1]+in[3]*in[2]))-(F_1_2S3*in[7]*in[3]);
		pPtilde[4] = (F_1_2*(in[6]*in[0]+in[5]*in[1]+in[4]*in[2]))-(F_1_2S3*in[7]*in[4]);
		pPtilde[5] = (F_1_2*(in[3]*in[0]+in[4]*in[1]-in[5]*in[2]))-(F_1_2S3*in[7]*in[5]);
		pPtilde[6] = (F_1_2*(in[4]*in[0]-in[3]*in[1]-in[6]*in[2]))-(F_1_2S3*in[7]*in[6]);
		pPtilde[7] = 0.0;

		pQ[2] = (F_1_2*(in[3]*in[3]+in[4]*in[4]-in[5]*in[5]-in[6]*in[6]));
		pQ[7] = (F_1_S3*(in[0]*in[0]+in[1]*in[1]+in[2]*in[2]-in[7]*in[7]))-(F_1_2S3*(in[3]*in[3]+in[4]*in[4]+in[5]*in[5]+in[6]*in[6]));
		pQ[0] = pQ[1] = pQ[3] = pQ[4] = pQ[5] = pQ[6] = 0.0;
		
		pT = (in[0]*(pQ[0]+2*pPtilde[0]))+(in[1]*(pQ[1]+2*pPtilde[1]))+(in[2]*(pQ[2]+2*pPtilde[2]))+(in[3]*(pQ[3]+2*pPtilde[3]))
			+(in[4]*(pQ[4]+2*pPtilde[4]))+(in[5]*(pQ[5]+2*pPtilde[5]))+(in[6]*(pQ[6]+2*pPtilde[6]))+(in[7]*(pQ[7]+2*pPtilde[7]));
			
		pS[0] = (F_1_2*(in[5]*(pQ[3]+2*pPtilde[3])+in[6]*(pQ[4]+2*pPtilde[4])))+(F_1_S3*in[7]*(pQ[0]+2*pPtilde[0]));
		pS[1] = (F_1_2*(in[5]*(pQ[4]+2*pPtilde[4])-in[6]*(pQ[3]+2*pPtilde[3])))+(F_1_S3*in[7]*(pQ[1]+2*pPtilde[1]));
		pS[2] = (F_1_S3*in[7]*(pQ[2]+2*pPtilde[2]));
		pS[3] = (F_1_2*(in[3]*(pQ[2]+2*pPtilde[2])+in[5]*(pQ[0]+2*pPtilde[0])-in[6]*(pQ[1]+2*pPtilde[1])))-(F_1_2S3*in[7]*(pQ[3]+2*pPtilde[3]));
		pS[4] = (F_1_2*(in[4]*(pQ[2]+2*pPtilde[2])+in[5]*(pQ[1]+2*pPtilde[1])+in[6]*(pQ[0]+2*pPtilde[0])))-(F_1_2S3*in[7]*(pQ[4]+2*pPtilde[4]));
		pS[5] = (F_1_2*(in[3]*(pQ[0]+2*pPtilde[0])+in[4]*(pQ[1]+2*pPtilde[1])-in[5]*(pQ[2]+2*pPtilde[2])))-(F_1_2S3*in[7]*(pQ[5]+2*pPtilde[5]));
		pS[6] = (F_1_2*(in[4]*(pQ[0]+2*pPtilde[0])-in[3]*(pQ[1]+2*pPtilde[1])-in[6]*(pQ[2]+2*pPtilde[2])))-(F_1_2S3*in[7]*(pQ[6]+2*pPtilde[6]));
		pS[7] = 0.0;

		pU[2] = (F_1_2*(in[3]*(pQ[3]+2*pPtilde[3])+in[4]*(pQ[4]+2*pPtilde[4])-in[5]*(pQ[5]+2*pPtilde[5])-in[6]*(pQ[6]+2*pPtilde[6])));
		pU[7] = (F_1_S3*(in[0]*(pQ[0]+2*pPtilde[0])+in[1]*(pQ[1]+2*pPtilde[1])+in[2]*(pQ[2]+2*pPtilde[2])-in[7]*(pQ[7]+2*pPtilde[7])))
			-(F_1_2S3*(in[3]*(pQ[3]+2*pPtilde[3])+in[4]*(pQ[4]+2*pPtilde[4])+in[5]*(pQ[5]+2*pPtilde[5])+in[6]*(pQ[6]+2*pPtilde[6])));
		pU[0] = pU[1] = pU[3] = pU[4] = pU[5] = pU[6] = 0.0;
		
		hmc_float eps_squared, eps_cubed;
		eps_squared = epsilon*epsilon;
		eps_cubed = eps_squared*epsilon;
		beta_0  = 1.0 - (eps_squared*pR/12.);
		gamma_0 = -eps_cubed*pT/72.;
		for(int il=0;il<8;il++){
			beta[il]  = -0.25*eps_squared*(pQ[il]+2.0*pPtilde[il]);
			gamma[il] = epsilon*in[il] - (eps_cubed/72.)*(2.*pR*in[il] + 3.*pU[il] + 6.*pS[il]);
		}
		#endif // EXPONENTIATE_ALGEBRA_ORDER_3
		
		// Phase 2: build a (generic 3x3) matrix from beta's and gamma's
		hmc_3x3matrix combination;
		construct_3x3_combination(beta_0, gamma_0, beta, gamma, combination);
		
		// Phase 3: project back onto SU(3). This requires knowledge whether one is using _RECONSTRUCT_TWELVE_ when filling "out"
		#ifdef _RECONSTRUCT_TWELVE_
			(*out)[0] = combination[0][0];
			(*out)[1] = combination[0][1];
			(*out)[2] = combination[0][2];
			(*out)[3] = combination[1][0];
			(*out)[4] = combination[1][1];
			(*out)[5] = combination[1][2];
		#else
			(*out)[0][0] = combination[0][0];
			(*out)[0][1] = combination[0][1];
			(*out)[0][2] = combination[0][2];
			(*out)[1][0] = combination[1][0];
			(*out)[1][1] = combination[1][1];
			(*out)[1][2] = combination[1][2];
			(*out)[2][0] = combination[2][0];
			(*out)[2][1] = combination[2][1];
			(*out)[2][2] = combination[2][2];
		#endif // _RECONSTRUCT_TWELVE_
				cout << "perform projection" << endl;
		project_su3(out);
		
	#else  // (ifndef) _EXPONENTIATE_ALGEBRA_ALL_ORDERS_
		// this is the case where one actually evaluates -as matrices- many orders of exp(i*e*P)=1+i*e*P+(1/2)(i*e*P) + ...
		// 1. create mat = (i epsilon)/2 p_i lambda_i    with lambda=gellmann matrix
		hmc_3x3matrix eMat;
		hmc_float halfeps = epsilon*F_1_2;
		eMat[0][0].re = 0.0;
		eMat[0][0].im = halfeps*(in[7]*F_1_S3+in[2]);
		eMat[0][1].re = halfeps*in[1];
		eMat[0][1].im = halfeps*in[0];
		eMat[0][2].re = halfeps*in[4];
		eMat[0][2].im = halfeps*in[3];
		eMat[1][0].re = -halfeps*in[1];;
		eMat[1][0].im = halfeps*in[0];;
		eMat[1][1].re = 0.0;
		eMat[1][1].im = halfeps*(in[7]*F_1_S3-in[2]);
		eMat[1][2].re = halfeps*in[6];
		eMat[1][2].im = halfeps*in[5];
		eMat[2][0].re = -halfeps*in[4];
		eMat[2][0].im = halfeps*in[3];
		eMat[2][1].re = -halfeps*in[6];
		eMat[2][1].im = halfeps*in[5];
		eMat[2][2].re = 0.0;
		eMat[2][2].im = -epsilon*in[7]*F_1_S3;
		// 2. start with the exp(...) expansion by using standard 3x3 sum and multiplication
		hmc_3x3matrix eRes, eCurPower, eNextPower, eLastResult;
		hmc_float eAccuracyCheck;
		set_to_3x3_identity(&eRes);
		set_to_3x3_identity(&eCurPower);
		hmc_float eCoefficient = 1.0;
		for(int power=1;power<_EXACT_EXPONENTIATION_MAX_POWER_;power++){
			multiply_3x3matrix(&eNextPower, &eMat, &eCurPower);
			copy_3x3_matrix(&eCurPower, &eNextPower);
			multiply_3x3matrix_by_real(&eCurPower, (1./(1.*power)));
			copy_3x3_matrix(&eLastResult, &eRes);
			add_3x3matrix(&eRes, &eRes, &eCurPower);
			absoluteDifference_3x3_matrix(&eAccuracyCheck, &eRes, &eLastResult);
			if(eAccuracyCheck < _EXACT_EXPONENTIATION_ACCURACY_){
				break;
			}
		}
		// 3. here I have the exponentiated matrix in 3x3 generic form (eRes), project it
		#ifdef _RECONSTRUCT_TWELVE_
			(*out)[0] = eRes[0][0];
			(*out)[1] = eRes[0][1];
			(*out)[2] = eRes[0][2];
			(*out)[3] = eRes[1][0];
			(*out)[4] = eRes[1][1];
			(*out)[5] = eRes[1][2];
		#else
			(*out)[0][0] = eRes[0][0];
			(*out)[0][1] = eRes[0][1];
			(*out)[0][2] = eRes[0][2];
			(*out)[1][0] = eRes[1][0];
			(*out)[1][1] = eRes[1][1];
			(*out)[1][2] = eRes[1][2];
			(*out)[2][0] = eRes[2][0];
			(*out)[2][1] = eRes[2][1];
			(*out)[2][2] = eRes[2][2];
		#endif // _RECONSTRUCT_TWELVE_
		project_su3(out);
	#endif // EXPONENTIATE_ALGEBRA_ALL_ORDERS_
	
	#endif // _USE_MORNGINGSTAR_PEARDON_

	
	return HMC_SUCCESS;
}
*/