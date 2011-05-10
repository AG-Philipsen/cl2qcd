#include "host_matrix.h"

hmc_error multiply_3x3matrix (hmc_3x3matrix *out, hmc_3x3matrix *p, hmc_3x3matrix *q)
{
	for(int i=0; i<3; i++) {
		for(int k=0; k<3; k++) {
			(*out)[i][k].re=0;
			(*out)[i][k].im=0;
			for(int j=0; j<3; j++) {
				hmc_complex tmp = complexmult(&(*p)[i][j],&(*q)[j][k]);
				complexaccumulate(&(*out)[i][k],&tmp);
			}
		}
	}
	return HMC_SUCCESS;
}


hmc_error add_3x3matrix (hmc_3x3matrix *out, hmc_3x3matrix *p, hmc_3x3matrix *q)
{
	for(int i=0; i<3; i++) {
		for(int k=0; k<3; k++) {
			(*out)[i][k] = complexadd(&(*p)[i][k],&(*q)[i][k]);
		}
	}
	return HMC_SUCCESS;
}


hmc_error subtract_3x3matrix (hmc_3x3matrix *out, hmc_3x3matrix *p, hmc_3x3matrix *q)
{
	for(int i=0; i<3; i++) {
		for(int k=0; k<3; k++) {
			(*out)[i][k] = complexsubtract(&(*p)[i][k],&(*q)[i][k]);
		}
	}
	return HMC_SUCCESS;
}


hmc_error set_to_3x3_identity(hmc_3x3matrix mat)
{
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

hmc_error copy_3x3_matrix(hmc_3x3matrix *dest, hmc_3x3matrix *src)
{
	// copies the values in src into dest
	(*dest)[0][0] = (*src)[0][0];
	(*dest)[0][1] = (*src)[0][1];
	(*dest)[0][2] = (*src)[0][2];
	(*dest)[1][0] = (*src)[1][0];
	(*dest)[1][1] = (*src)[1][1];
	(*dest)[1][2] = (*src)[1][2];
	(*dest)[2][0] = (*src)[2][0];
	(*dest)[2][1] = (*src)[2][1];
	(*dest)[2][2] = (*src)[2][2];
	return HMC_SUCCESS;
}

hmc_error multiply_3x3matrix_by_real(hmc_3x3matrix *mat, hmc_float factor)
{
	// applies a real factor to all elements of matrix mat
	complexmult_real(&(*mat)[0][0], &factor);
	complexmult_real(&(*mat)[0][1], &factor);
	complexmult_real(&(*mat)[0][2], &factor);
	complexmult_real(&(*mat)[1][0], &factor);
	complexmult_real(&(*mat)[1][1], &factor);
	complexmult_real(&(*mat)[1][2], &factor);
	complexmult_real(&(*mat)[2][0], &factor);
	complexmult_real(&(*mat)[2][1], &factor);
	complexmult_real(&(*mat)[2][2], &factor);
	return HMC_SUCCESS;
}

int absoluteDifference_3x3_matrix(hmc_float *result, hmc_3x3matrix *mat1, hmc_3x3matrix *mat2)
{
	*result=0.0;
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++) {
			*result += fabs((*mat1)[i][j].re - (*mat2)[i][j].re);
			*result += fabs((*mat1)[i][j].im - (*mat2)[i][j].im);
		}
	return HMC_SUCCESS;
}

hmc_error accumulate_su3matrix_3x3_add(hmc_3x3matrix *out, hmc_su3matrix *q){
#ifdef _RECONSTRUCT_TWELVE_
  for (int i=0; i< NC-1; i++)
  {
    for (int k=0; k<NC; k++)
    {
      complexaccumulate(&(*out)[i][k], &(*q)[i+(NC-1)*k]);
    }
    for (int k=0; k<NC; k++)
      {
	hmc_complex tmp = reconstruct_su3(q, k); 
	complexaccumulate(&(*out)[2][k], &(tmp));
      }
  }
#else

  for(int i=0; i<NC; i++) {
    for(int k=0; k<NC; k++) {
      complexaccumulate(&(*out)[i][k],&(*q)[i][k]);
    }
  }
#endif
  return HMC_SUCCESS;
}

hmc_error trace_3x3matrix (hmc_complex * out, hmc_3x3matrix *q){
  (*out).re = (*q[0][0]).re;
  (*out).im = (*q[0][0]).im;
  (*out).re += (*q[1][1]).re;
  (*out).im += (*q[1][1]).im;
  (*out).re += (*q[2][2]).re;
  (*out).im += (*q[2][2]).im;
  
  return HMC_SUCCESS;
}

hmc_error adjoint_3x3matrix (hmc_3x3matrix * out, hmc_3x3matrix *q){
  for(int a=0; a<3; a++) {
    for(int b=0; b<3; b++) {
      (*out)[a][b] = complexconj(&(*q)[b][a]);
    }
  }  
  return HMC_SUCCESS;
}

hmc_error multiply_generator_3x3matrix (hmc_3x3matrix * out, int gen_index, hmc_3x3matrix *in){
	// useful constants:
	// F_1_2   = 1/2
	// F_1_2S3 = 1/(2*sqrt(3))
	// F_1_S3  = 1/sqrt(3)

	// SL: not yet tested!

	switch(gen_index){
		case 1:
			(*out)[0][0].re = +F_1_2*(*in)[1][0].re;
			(*out)[0][0].im = +F_1_2*(*in)[1][0].im;
			(*out)[0][1].re = +F_1_2*(*in)[1][1].re;
			(*out)[0][1].im = +F_1_2*(*in)[1][1].im;
			(*out)[0][2].re = +F_1_2*(*in)[1][2].re;
			(*out)[0][2].im = +F_1_2*(*in)[1][2].im;
			(*out)[1][0].re = +F_1_2*(*in)[0][0].re;
			(*out)[1][0].im = +F_1_2*(*in)[0][0].im;
			(*out)[1][1].re = +F_1_2*(*in)[0][1].re;
			(*out)[1][1].im = +F_1_2*(*in)[0][1].im;
			(*out)[1][2].re = +F_1_2*(*in)[0][2].re;
			(*out)[1][2].im = +F_1_2*(*in)[0][2].im;
			(*out)[2][0].re = 0.0;
			(*out)[2][0].im = 0.0;
			(*out)[2][1].re = 0.0;
			(*out)[2][1].im = 0.0;
			(*out)[2][2].re = 0.0;
			(*out)[2][2].im = 0.0;
			break;
		case 2:
			(*out)[0][0].re = +F_1_2*(*in)[1][0].im;
			(*out)[0][0].im = -F_1_2*(*in)[1][0].re;
			(*out)[0][1].re = +F_1_2*(*in)[1][1].im;
			(*out)[0][1].im = -F_1_2*(*in)[1][1].re;
			(*out)[0][2].re = +F_1_2*(*in)[1][2].im;
			(*out)[0][2].im = -F_1_2*(*in)[1][2].re;
			(*out)[1][0].re = -F_1_2*(*in)[0][0].im;
			(*out)[1][0].im = +F_1_2*(*in)[0][0].re;
			(*out)[1][1].re = -F_1_2*(*in)[0][1].im;
			(*out)[1][1].im = +F_1_2*(*in)[0][1].re;
			(*out)[1][2].re = -F_1_2*(*in)[0][2].im;
			(*out)[1][2].im = +F_1_2*(*in)[0][2].re;
			(*out)[2][0].re = 0.0;
			(*out)[2][0].im = 0.0;
			(*out)[2][1].re = 0.0;
			(*out)[2][1].im = 0.0;
			(*out)[2][2].re = 0.0;
			(*out)[2][2].im = 0.0;
			break;
		case 3:
			(*out)[0][0].re = +F_1_2*(*in)[0][0].re;
			(*out)[0][0].im = +F_1_2*(*in)[0][0].im;
			(*out)[0][1].re = +F_1_2*(*in)[0][1].re;
			(*out)[0][1].im = +F_1_2*(*in)[0][1].im;
			(*out)[0][2].re = +F_1_2*(*in)[0][2].re;
			(*out)[0][2].im = +F_1_2*(*in)[0][2].im;
			(*out)[1][0].re = -F_1_2*(*in)[1][0].re;
			(*out)[1][0].im = -F_1_2*(*in)[1][0].im;
			(*out)[1][1].re = -F_1_2*(*in)[1][1].re;
			(*out)[1][1].im = -F_1_2*(*in)[1][1].im;
			(*out)[1][2].re = -F_1_2*(*in)[1][2].re;
			(*out)[1][2].im = -F_1_2*(*in)[1][2].im;
			(*out)[2][0].re = 0.0;
			(*out)[2][0].im = 0.0;
			(*out)[2][1].re = 0.0;
			(*out)[2][1].im = 0.0;
			(*out)[2][2].re = 0.0;
			(*out)[2][2].im = 0.0;
			break;
		case 4:
			(*out)[0][0].re = +F_1_2*(*in)[2][0].re;
			(*out)[0][0].im = +F_1_2*(*in)[2][0].im;
			(*out)[0][1].re = +F_1_2*(*in)[2][1].re;
			(*out)[0][1].im = +F_1_2*(*in)[2][1].im;
			(*out)[0][2].re = +F_1_2*(*in)[2][2].re;
			(*out)[0][2].im = +F_1_2*(*in)[2][2].im;
			(*out)[1][0].re = 0.0;
			(*out)[1][0].im = 0.0;
			(*out)[1][1].re = 0.0;
			(*out)[1][1].im = 0.0;
			(*out)[1][2].re = 0.0;
			(*out)[1][2].im = 0.0;
			(*out)[2][0].re = +F_1_2*(*in)[0][0].re;
			(*out)[2][0].im = +F_1_2*(*in)[0][0].im;
			(*out)[2][1].re = +F_1_2*(*in)[0][1].re;
			(*out)[2][1].im = +F_1_2*(*in)[0][1].im;
			(*out)[2][2].re = +F_1_2*(*in)[0][2].re;
			(*out)[2][2].im = +F_1_2*(*in)[0][2].im;
			break;
		case 5:
			(*out)[0][0].re = +F_1_2*(*in)[2][0].im;
			(*out)[0][0].im = -F_1_2*(*in)[2][0].re;
			(*out)[0][1].re = +F_1_2*(*in)[2][1].im;
			(*out)[0][1].im = -F_1_2*(*in)[2][1].re;
			(*out)[0][2].re = +F_1_2*(*in)[2][2].im;
			(*out)[0][2].im = -F_1_2*(*in)[2][2].re;
			(*out)[1][0].re = 0.0;
			(*out)[1][0].im = 0.0;
			(*out)[1][1].re = 0.0;
			(*out)[1][1].im = 0.0;
			(*out)[1][2].re = 0.0;
			(*out)[1][2].im = 0.0;
			(*out)[2][0].re = -F_1_2*(*in)[0][0].im;
			(*out)[2][0].im = +F_1_2*(*in)[0][0].re;
			(*out)[2][1].re = -F_1_2*(*in)[0][1].im;
			(*out)[2][1].im = +F_1_2*(*in)[0][1].re;
			(*out)[2][2].re = -F_1_2*(*in)[0][2].im;
			(*out)[2][2].im = +F_1_2*(*in)[0][2].re;
			break;
		case 6:
			(*out)[0][0].re = 0.0;
			(*out)[0][0].im = 0.0;
			(*out)[0][1].re = 0.0;
			(*out)[0][1].im = 0.0;
			(*out)[0][2].re = 0.0;
			(*out)[0][2].im = 0.0;
			(*out)[1][0].re = +F_1_2*(*in)[2][0].re;
			(*out)[1][0].im = +F_1_2*(*in)[2][0].im;
			(*out)[1][1].re = +F_1_2*(*in)[2][1].re;
			(*out)[1][1].im = +F_1_2*(*in)[2][1].im;
			(*out)[1][2].re = +F_1_2*(*in)[2][2].re;
			(*out)[1][2].im = +F_1_2*(*in)[2][2].im;
			(*out)[2][0].re = +F_1_2*(*in)[1][0].re;
			(*out)[2][0].im = +F_1_2*(*in)[1][0].im;
			(*out)[2][1].re = +F_1_2*(*in)[1][1].re;
			(*out)[2][1].im = +F_1_2*(*in)[1][1].im;
			(*out)[2][2].re = +F_1_2*(*in)[1][2].re;
			(*out)[2][2].im = +F_1_2*(*in)[1][2].im;
			break;
		case 7:
			(*out)[0][0].re = 0.0;
			(*out)[0][0].im = 0.0;
			(*out)[0][1].re = 0.0;
			(*out)[0][1].im = 0.0;
			(*out)[0][2].re = 0.0;
			(*out)[0][2].im = 0.0;
			(*out)[1][0].re = +F_1_2*(*in)[2][0].im;
			(*out)[1][0].im = -F_1_2*(*in)[2][0].re;
			(*out)[1][1].re = +F_1_2*(*in)[2][1].im;
			(*out)[1][1].im = -F_1_2*(*in)[2][1].re;
			(*out)[1][2].re = +F_1_2*(*in)[2][2].im;
			(*out)[1][2].im = -F_1_2*(*in)[2][2].re;
			(*out)[2][0].re = -F_1_2*(*in)[1][0].im;
			(*out)[2][0].im = +F_1_2*(*in)[1][0].re;
			(*out)[2][1].re = -F_1_2*(*in)[1][1].im;
			(*out)[2][1].im = +F_1_2*(*in)[1][1].re;
			(*out)[2][2].re = -F_1_2*(*in)[1][2].im;
			(*out)[2][2].im = +F_1_2*(*in)[1][2].re;
			break;
		case 8:
			(*out)[0][0].re = +F_1_2S3*(*in)[0][0].re;
			(*out)[0][0].im = +F_1_2S3*(*in)[0][0].im;
			(*out)[0][1].re = +F_1_2S3*(*in)[0][1].re;
			(*out)[0][1].im = +F_1_2S3*(*in)[0][1].im;
			(*out)[0][2].re = +F_1_2S3*(*in)[0][2].re;
			(*out)[0][2].im = +F_1_2S3*(*in)[0][2].im;
			(*out)[1][0].re = +F_1_2S3*(*in)[1][0].re;
			(*out)[1][0].im = +F_1_2S3*(*in)[1][0].im;
			(*out)[1][1].re = +F_1_2S3*(*in)[1][1].re;
			(*out)[1][1].im = +F_1_2S3*(*in)[1][1].im;
			(*out)[1][2].re = +F_1_2S3*(*in)[1][2].re;
			(*out)[1][2].im = +F_1_2S3*(*in)[1][2].im;
			(*out)[2][0].re = -F_1_S3*(*in)[2][0].re;
			(*out)[2][0].im = -F_1_S3*(*in)[2][0].im;
			(*out)[2][1].re = -F_1_S3*(*in)[2][1].re;
			(*out)[2][1].im = -F_1_S3*(*in)[2][1].im;
			(*out)[2][2].re = -F_1_S3*(*in)[2][2].re;
			(*out)[2][2].im = -F_1_S3*(*in)[2][2].im;
			break;
		default:
			return HMC_INVALID_GENERATOR_INDEX;
	}
	return HMC_SUCCESS;
}

hmc_error multiply_3x3matrix_generator (hmc_3x3matrix * out, hmc_3x3matrix *in, int gen_index){
	// useful constants:
	// F_1_2   = 1/2
	// F_1_2S3 = 1/(2*sqrt(3))
	// F_1_S3  = 1/sqrt(3)

	// SL: not yet tested!

	switch(gen_index){
		case 1:
			(*out)[0][0].re = +F_1_2*(*in)[0][1].re;
			(*out)[0][0].im = +F_1_2*(*in)[0][1].im;
			(*out)[0][1].re = +F_1_2*(*in)[0][0].re;
			(*out)[0][1].im = +F_1_2*(*in)[0][0].im;
			(*out)[0][2].re = 0.0;
			(*out)[0][2].im = 0.0;
			(*out)[1][0].re = +F_1_2*(*in)[1][1].re;
			(*out)[1][0].im = +F_1_2*(*in)[1][1].im;
			(*out)[1][1].re = +F_1_2*(*in)[1][0].re;
			(*out)[1][1].im = +F_1_2*(*in)[1][0].im;
			(*out)[1][2].re = 0.0;
			(*out)[1][2].im = 0.0;
			(*out)[2][0].re = +F_1_2*(*in)[2][1].re;
			(*out)[2][0].im = +F_1_2*(*in)[2][1].im;
			(*out)[2][1].re = +F_1_2*(*in)[2][0].re;
			(*out)[2][1].im = +F_1_2*(*in)[2][0].im;
			(*out)[2][2].re = 0.0;
			(*out)[2][2].im = 0.0;
			break;
		case 2:
			(*out)[0][0].re = -F_1_2*(*in)[0][1].im;
			(*out)[0][0].im = +F_1_2*(*in)[0][1].re;
			(*out)[0][1].re = +F_1_2*(*in)[0][0].im;
			(*out)[0][1].im = -F_1_2*(*in)[0][0].re;
			(*out)[0][2].re = 0.0;
			(*out)[0][2].im = 0.0;
			(*out)[1][0].re = -F_1_2*(*in)[1][1].im;
			(*out)[1][0].im = +F_1_2*(*in)[1][1].re;
			(*out)[1][1].re = +F_1_2*(*in)[1][0].im;
			(*out)[1][1].im = -F_1_2*(*in)[1][0].re;
			(*out)[1][2].re = 0.0;
			(*out)[1][2].im = 0.0;
			(*out)[2][0].re = -F_1_2*(*in)[2][1].im;
			(*out)[2][0].im = +F_1_2*(*in)[2][1].re;
			(*out)[2][1].re = +F_1_2*(*in)[2][0].im;
			(*out)[2][1].im = -F_1_2*(*in)[2][0].re;
			(*out)[2][2].re = 0.0;
			(*out)[2][2].im = 0.0;
			break;
		case 3:
			(*out)[0][0].re = +F_1_2*(*in)[0][0].re;
			(*out)[0][0].im = +F_1_2*(*in)[0][0].im;
			(*out)[0][1].re = -F_1_2*(*in)[0][1].re;
			(*out)[0][1].im = -F_1_2*(*in)[0][1].im;
			(*out)[0][2].re = 0.0;
			(*out)[0][2].im = 0.0;
			(*out)[1][0].re = +F_1_2*(*in)[1][0].re;
			(*out)[1][0].im = +F_1_2*(*in)[1][0].im;
			(*out)[1][1].re = -F_1_2*(*in)[1][1].re;
			(*out)[1][1].im = -F_1_2*(*in)[1][1].im;
			(*out)[1][2].re = 0.0;
			(*out)[1][2].im = 0.0;
			(*out)[2][0].re = +F_1_2*(*in)[2][0].re;
			(*out)[2][0].im = +F_1_2*(*in)[2][0].im;
			(*out)[2][1].re = -F_1_2*(*in)[2][1].re;
			(*out)[2][1].im = -F_1_2*(*in)[2][1].im;
			(*out)[2][2].re = 0.0;
			(*out)[2][2].im = 0.0;
			break;
		case 4:
			(*out)[0][0].re = +F_1_2*(*in)[0][2].re;
			(*out)[0][0].im = +F_1_2*(*in)[0][2].im;
			(*out)[0][1].re = 0.0;
			(*out)[0][1].im = 0.0;
			(*out)[0][2].re = +F_1_2*(*in)[0][0].re;
			(*out)[0][2].im = +F_1_2*(*in)[0][0].im;
			(*out)[1][0].re = +F_1_2*(*in)[1][2].re;
			(*out)[1][0].im = +F_1_2*(*in)[1][2].im;
			(*out)[1][1].re = 0.0;
			(*out)[1][1].im = 0.0;
			(*out)[1][2].re = +F_1_2*(*in)[1][0].re;
			(*out)[1][2].im = +F_1_2*(*in)[1][0].im;
			(*out)[2][0].re = +F_1_2*(*in)[2][2].re;
			(*out)[2][0].im = +F_1_2*(*in)[2][2].im;
			(*out)[2][1].re = 0.0;
			(*out)[2][1].im = 0.0;
			(*out)[2][2].re = +F_1_2*(*in)[2][0].re;
			(*out)[2][2].im = +F_1_2*(*in)[2][0].im;
			break;
		case 5:
			(*out)[0][0].re = -F_1_2*(*in)[0][2].im;
			(*out)[0][0].im = +F_1_2*(*in)[0][2].re;
			(*out)[0][1].re = 0.0;
			(*out)[0][1].im = 0.0;
			(*out)[0][2].re = +F_1_2*(*in)[0][0].im;
			(*out)[0][2].im = -F_1_2*(*in)[0][0].re;
			(*out)[1][0].re = -F_1_2*(*in)[1][2].im;
			(*out)[1][0].im = +F_1_2*(*in)[1][2].re;
			(*out)[1][1].re = 0.0;
			(*out)[1][1].im = 0.0;
			(*out)[1][2].re = +F_1_2*(*in)[1][0].im;
			(*out)[1][2].im = -F_1_2*(*in)[1][0].re;
			(*out)[2][0].re = -F_1_2*(*in)[2][2].im;
			(*out)[2][0].im = +F_1_2*(*in)[2][2].re;
			(*out)[2][1].re = 0.0;
			(*out)[2][1].im = 0.0;
			(*out)[2][2].re = +F_1_2*(*in)[2][0].im;
			(*out)[2][2].im = -F_1_2*(*in)[2][0].re;
			break;
		case 6:
			(*out)[0][0].re = 0.0;
			(*out)[0][0].im = 0.0;
			(*out)[0][1].re = F_1_2*(*in)[0][2].re;
			(*out)[0][1].im = F_1_2*(*in)[0][2].im;
			(*out)[0][2].re = F_1_2*(*in)[0][1].re;
			(*out)[0][2].im = F_1_2*(*in)[0][1].im;
			(*out)[1][0].re = 0.0;
			(*out)[1][0].im = 0.0;
			(*out)[1][1].re = F_1_2*(*in)[1][2].re;
			(*out)[1][1].im = F_1_2*(*in)[1][2].im;
			(*out)[1][2].re = F_1_2*(*in)[1][1].re;
			(*out)[1][2].im = F_1_2*(*in)[1][1].im;
			(*out)[2][0].re = 0.0;
			(*out)[2][0].im = 0.0;
			(*out)[2][1].re = F_1_2*(*in)[2][2].re;
			(*out)[2][1].im = F_1_2*(*in)[2][2].im;
			(*out)[2][2].re = F_1_2*(*in)[2][1].re;
			(*out)[2][2].im = F_1_2*(*in)[2][1].im;
			break;
		case 7:
			(*out)[0][0].re = 0.0;
			(*out)[0][0].im = 0.0;
			(*out)[0][1].re = -F_1_2*(*in)[0][2].im;
			(*out)[0][1].im = +F_1_2*(*in)[0][2].re;
			(*out)[0][2].re = +F_1_2*(*in)[0][1].im;
			(*out)[0][2].im = -F_1_2*(*in)[0][1].re;
			(*out)[1][0].re = 0.0;
			(*out)[1][0].im = 0.0;
			(*out)[1][1].re = -F_1_2*(*in)[1][2].im;
			(*out)[1][1].im = +F_1_2*(*in)[1][2].re;
			(*out)[1][2].re = +F_1_2*(*in)[1][1].im;
			(*out)[1][2].im = -F_1_2*(*in)[1][1].re;
			(*out)[2][0].re = 0.0;
			(*out)[2][0].im = 0.0;
			(*out)[2][1].re = -F_1_2*(*in)[2][2].im;
			(*out)[2][1].im = +F_1_2*(*in)[2][2].re;
			(*out)[2][2].re = +F_1_2*(*in)[2][1].im;
			(*out)[2][2].im = -F_1_2*(*in)[2][1].re;
			break;
		case 8:
			(*out)[0][0].re = +F_1_2S3*(*in)[0][0].re;
			(*out)[0][0].im = +F_1_2S3*(*in)[0][0].im;
			(*out)[0][1].re = +F_1_2S3*(*in)[0][1].re;
			(*out)[0][1].im = +F_1_2S3*(*in)[0][1].im;
			(*out)[0][2].re = -F_1_S3*(*in)[0][2].re;
			(*out)[0][2].im = -F_1_S3*(*in)[0][2].im;
			(*out)[1][0].re = +F_1_2S3*(*in)[1][0].re;
			(*out)[1][0].im = +F_1_2S3*(*in)[1][0].im;
			(*out)[1][1].re = +F_1_2S3*(*in)[1][1].re;
			(*out)[1][1].im = +F_1_2S3*(*in)[1][1].im;
			(*out)[1][2].re = -F_1_S3*(*in)[1][2].re;
			(*out)[1][2].im = -F_1_S3*(*in)[1][2].im;
			(*out)[2][0].re = +F_1_2S3*(*in)[2][0].re;
			(*out)[2][0].im = +F_1_2S3*(*in)[2][0].im;
			(*out)[2][1].re = +F_1_2S3*(*in)[2][1].re;
			(*out)[2][1].im = +F_1_2S3*(*in)[2][1].im;
			(*out)[2][2].re = -F_1_S3*(*in)[2][2].re;
			(*out)[2][2].im = -F_1_S3*(*in)[2][2].im;
			break;
			default:
				return HMC_INVALID_GENERATOR_INDEX;
	}
	return HMC_SUCCESS;
}

hmc_error multiply_generator_su3matrix (hmc_3x3matrix * out, int gen_index, hmc_su3matrix *in){
	// if needed, construct the full 3x3 matrix and then invoke the general_3x3 version of this
	// SL: not yet tested!
	#ifdef _RECONSTRUCT_TWELVE_
		hmc_3x3matrix complete_reconstructed;
		complete_reconstructed[0][0]=*in[0];
		complete_reconstructed[0][1]=*in[1];
		complete_reconstructed[0][2]=*in[2];
		complete_reconstructed[1][0]=*in[3];
		complete_reconstructed[1][1]=*in[4];
		complete_reconstructed[1][2]=*in[5];
		complete_reconstructed[2][0]=reconstruct_su3(in,0);
		complete_reconstructed[2][1]=reconstruct_su3(in,1);
		complete_reconstructed[2][2]=reconstruct_su3(in,2);
		#define _MULTIPLY_SU3MATRIX_GENERATOR_TARGET_ (&complete_reconstructed)
	#else
		#define _MULTIPLY_SU3MATRIX_GENERATOR_TARGET_ (in)
	#endif
	return multiply_generator_3x3matrix(out, gen_index, _MULTIPLY_SU3MATRIX_GENERATOR_TARGET_);
}

hmc_error multiply_su3matrix_generator (hmc_3x3matrix * out, hmc_su3matrix *in, int gen_index){
	// if needed, construct the full 3x3 matrix and then invoke the general_3x3 version of this
	// SL: not yet tested!
	#ifdef _RECONSTRUCT_TWELVE_
		hmc_3x3matrix complete_reconstructed;
		complete_reconstructed[0][0]=*in[0];
		complete_reconstructed[0][1]=*in[1];
		complete_reconstructed[0][2]=*in[2];
		complete_reconstructed[1][0]=*in[3];
		complete_reconstructed[1][1]=*in[4];
		complete_reconstructed[1][2]=*in[5];
		complete_reconstructed[2][0]=reconstruct_su3(in,0);
		complete_reconstructed[2][1]=reconstruct_su3(in,1);
		complete_reconstructed[2][2]=reconstruct_su3(in,2);
		#define _MULTIPLY_SU3MATRIX_GENERATOR_TARGET_ (&complete_reconstructed)
	#else
		#define _MULTIPLY_SU3MATRIX_GENERATOR_TARGET_ (in)
	#endif
	return multiply_3x3matrix_generator(out, _MULTIPLY_SU3MATRIX_GENERATOR_TARGET_, gen_index);
}


hmc_error construct_3x3_combination(hmc_float beta_0, hmc_float gamma_0, hmc_float beta[], hmc_float gamma[], hmc_3x3matrix out){
	// called by build_su3matrix_by_exponentiation in case of "smart" approach
	// takes the 2*(8+1) real parameters beta_0, gamma_0, beta[8], gamma[8] and compiles
	// all components of the generic 3x3 complex matrix that is the linear combination of identity+generators
	hmc_float redb8 = beta[7] * F_1_2S3;
	hmc_float redg8 = gamma[7]* F_1_2S3;
	out[0][0].re = beta_0 + 0.5*beta[2] + redb8;
	out[0][0].im = gamma_0 + 0.5*gamma[2] + redg8;
	out[0][1].re = 0.5*(beta[0]+gamma[1]);
	out[0][1].im = 0.5*(gamma[0]-beta[1]);
	out[0][2].re = 0.5*(beta[3]+gamma[4]);
	out[0][2].im = 0.5*(gamma[3]-beta[4]);
	out[1][0].re = 0.5*(beta[0]-gamma[1]);
	out[1][0].im = 0.5*(gamma[0]+beta[1]);
	out[1][1].re = beta_0 - 0.5*beta[2] + redb8;
	out[1][1].im = gamma_0 - 0.5*gamma[2] + redg8;
	out[1][2].re = 0.5*(beta[5]+gamma[6]);
	out[1][2].im = 0.5*(gamma[5]-beta[6]);
	out[2][0].re = 0.5*(beta[3]-gamma[4]);
	out[2][0].im = 0.5*(gamma[3]+beta[4]);
	out[2][1].re = 0.5*(beta[5]-gamma[6]);
	out[2][1].im = 0.5*(gamma[5]+beta[6]);
	out[2][2].re = beta_0 + 2*redb8;
	out[2][2].im = gamma_0 - 2*redg8;
	return HMC_SUCCESS;
}