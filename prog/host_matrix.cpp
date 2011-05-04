#include "host_matrix.h"

hmc_error multiply_3x3matrix (hmc_3x3matrix *out, hmc_3x3matrix *p, hmc_3x3matrix *q){
  for(int i=0; i<3; i++) {
    for(int k=0; k<3; k++) {
      (*out)[i][k].re=0;
      (*out)[i][k].im=0;
      for(int j=0;j<3;j++) {
	hmc_complex tmp = complexmult(&(*p)[i][j],&(*q)[j][k]);
	complexaccumulate(&(*out)[i][k],&tmp);
      }
    }
  }
  return HMC_SUCCESS;
}


hmc_error add_3x3matrix (hmc_3x3matrix *out, hmc_3x3matrix *p, hmc_3x3matrix *q){
  for(int i=0; i<3; i++) {
    for(int k=0; k<3; k++) {
      (*out)[i][k] = complexadd(&(*p)[i][k],&(*q)[i][k]);
    }
  }
  return HMC_SUCCESS;
}


hmc_error substract_3x3matrix (hmc_3x3matrix *out, hmc_3x3matrix *p, hmc_3x3matrix *q){
  for(int i=0; i<3; i++) {
    for(int k=0; k<3; k++) {
      (*out)[i][k] = complexsubtract(&(*p)[i][k],&(*q)[i][k]);
    }
  }
  return HMC_SUCCESS;
}
