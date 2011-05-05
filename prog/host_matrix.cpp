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


hmc_error subtract_3x3matrix (hmc_3x3matrix *out, hmc_3x3matrix *p, hmc_3x3matrix *q){
  for(int i=0; i<3; i++) {
    for(int k=0; k<3; k++) {
      (*out)[i][k] = complexsubtract(&(*p)[i][k],&(*q)[i][k]);
    }
  }
  return HMC_SUCCESS;
}


hmc_error set_to_3x3_identity(hmc_3x3matrix mat){
	// simply sets to identity a generic 3x3 complex matrix
	mat[0][0].re=1.0;
	mat[0][0].im=0.0;
	mat[0][1].re=mat[0][1].im=0.0;
	mat[0][2].re=mat[0][2].im=0.0;
	mat[1][0].re=mat[1][0].im=0.0;
	mat[1][1].re=1.0;
	mat[1][1].im=0.0;
	mat[1][2].re=mat[1][2].im=0.0;
	mat[2][0].re=mat[2][0].im=0.0;
	mat[2][1].re=mat[2][1].im=0.0;
	mat[2][2].re=1.0;
	mat[2][2].im=0.0;
	return HMC_SUCCESS;
}

hmc_error copy_3x3_matrix(hmc_3x3matrix *dest, hmc_3x3matrix *src){
	// copies the values in src into dest
	dest[0][0]=src[0][0];
	dest[0][1]=src[0][1];
	dest[0][2]=src[0][2];
	dest[1][0]=src[1][0];
	dest[1][1]=src[1][1];
	dest[1][2]=src[1][2];
	dest[2][0]=src[2][0];
	dest[2][1]=src[2][1];
	dest[2][2]=src[2][2];
	return HMC_SUCCESS;
}

hmc_error multiply_3x3matrix_by_real(hmc_3x3matrix *mat, hmc_float factor){
	// applies a real factor to all elements of matrix mat
	complexmult_real(*mat[0][0], *factor);
	complexmult_real(*mat[0][1], *factor);
	complexmult_real(*mat[0][2], *factor);
	complexmult_real(*mat[1][0], *factor);
	complexmult_real(*mat[1][1], *factor);
	complexmult_real(*mat[1][2], *factor);
	complexmult_real(*mat[2][0], *factor);
	complexmult_real(*mat[2][1], *factor);
	complexmult_real(*mat[2][2], *factor);
	return HMC_SUCCESS;
}

absoluteDifference_3x3_matrix(hmc_float *result, hmc_3x3matrix *mat1, hmc_3x3matrix *mat2){
	result=0.0;
	for(int i=;li<3;i++)for(j=0;j<3;j++){
		result += fabs(mat1[i][j].re-mat2[i][j].re);
		result += fabs(mat1[i][j].im-mat2[i][j].im);
	}
	return HMC_SUCCESS;
}