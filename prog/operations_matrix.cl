/** @file
 * Device code implementing 3x3 matrices
 */

//operations_matrix.cl

Matrix3x3 zero_matrix3x3 (){
    Matrix3x3 out;
    out.e00.re = 0.;
    out.e00.im = 0.;
    out.e01.re = 0.;
    out.e01.im = 0.;
    out.e02.re = 0.;
    out.e02.im = 0.;
    out.e10.re = 0.;
    out.e10.im = 0.;
    out.e11.re = 0.;
    out.e11.im = 0.;
    out.e12.re = 0.;
    out.e12.im = 0.;
    out.e20.re = 0.;
    out.e20.im = 0.;
    out.e21.re = 0.;
    out.e21.im = 0.;
    out.e22.re = 0.;
    out.e22.im = 0.;

    return out;
}

Matrix3x3 multiply_matrix3x3 (const Matrix3x3 p, const Matrix3x3 q)
{
    Matrix3x3 out;
    out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e10.re + p.e20.re * q.e02.re
	       - p.e00.im * q.e00.im - p.e01.im * q.e10.im - p.e20.im * q.e02.im;
    out.e00.im = p.e00.re * q.e00.im + p.e01.re * q.e10.im + p.e20.re * q.e02.im
	       + p.e00.im * q.e00.re + p.e01.im * q.e10.re + p.e20.im * q.e02.re;

    out.e01.re = p.e00.re * q.e01.re + p.e01.re * q.e11.re + p.e02.re * q.e21.re
	       - p.e00.im * q.e01.im - p.e01.im * q.e11.im - p.e02.im * q.e21.im;
    out.e01.im = p.e00.re * q.e01.im + p.e01.re * q.e11.im + p.e02.re * q.e21.im
	       + p.e00.im * q.e01.re + p.e01.im * q.e11.re + p.e02.im * q.e21.re;
    
    out.e02.re = p.e00.re * q.e02.re + p.e01.re * q.e12.re + p.e02.re * q.e22.re
	       - p.e00.im * q.e02.im - p.e01.im * q.e12.im - p.e02.im * q.e22.im;
    out.e02.im = p.e00.re * q.e02.im + p.e01.re * q.e12.im + p.e02.re * q.e22.im
	       + p.e00.im * q.e02.re + p.e01.im * q.e12.re + p.e02.im * q.e22.re;
	       
    out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e10.re + p.e12.re * q.e20.re
	       - p.e10.im * q.e00.im - p.e11.im * q.e10.im - p.e12.im * q.e20.im;
    out.e10.im = p.e10.re * q.e00.im + p.e11.re * q.e10.im + p.e12.re * q.e20.im
	       + p.e10.im * q.e00.re + p.e11.im * q.e10.re + p.e12.im * q.e20.re;
	       
    out.e11.re = p.e10.re * q.e01.re + p.e11.re * q.e11.re + p.e12.re * q.e21.re
	       - p.e10.im * q.e01.im - p.e11.im * q.e11.im - p.e12.im * q.e21.im;
    out.e11.im = p.e10.re * q.e01.im + p.e11.re * q.e11.im + p.e12.re * q.e21.im
	       + p.e10.im * q.e01.re + p.e11.im * q.e11.re + p.e12.im * q.e21.re;
	       
    out.e12.re = p.e10.re * q.e02.re + p.e11.re * q.e12.re + p.e12.re * q.e22.re
	       - p.e10.im * q.e02.im - p.e11.im * q.e12.im - p.e12.im * q.e22.im;
    out.e12.im = p.e10.re * q.e02.im + p.e11.re * q.e12.im + p.e12.re * q.e22.im
	       + p.e10.im * q.e02.re + p.e11.im * q.e12.re + p.e12.im * q.e22.re;	       

    out.e20.re = p.e20.re * q.e00.re + p.e21.re * q.e10.re + p.e22.re * q.e20.re
	       - p.e20.im * q.e00.im - p.e21.im * q.e10.im - p.e22.im * q.e20.im;
    out.e20.im = p.e20.re * q.e00.im + p.e21.re * q.e10.im + p.e22.re * q.e20.im
	       + p.e20.im * q.e00.re + p.e21.im * q.e10.re + p.e22.im * q.e20.re;
	       
    out.e21.re = p.e20.re * q.e01.re + p.e21.re * q.e11.re + p.e22.re * q.e21.re
	       - p.e20.im * q.e01.im - p.e21.im * q.e11.im - p.e22.im * q.e21.im;
    out.e21.im = p.e20.re * q.e01.im + p.e21.re * q.e11.im + p.e22.re * q.e21.im
	       + p.e20.im * q.e01.re + p.e21.im * q.e11.re + p.e22.im * q.e21.re;
	       
    out.e22.re = p.e20.re * q.e02.re + p.e21.re * q.e12.re + p.e22.re * q.e22.re
	       - p.e20.im * q.e02.im - p.e21.im * q.e12.im - p.e22.im * q.e22.im;
    out.e22.im = p.e20.re * q.e02.im + p.e21.re * q.e12.im + p.e22.re * q.e22.im
	       + p.e20.im * q.e02.re + p.e21.im * q.e12.re + p.e22.im * q.e22.re;
    
    return out;
}

Matrix3x3 add_matrix3x3 ( const Matrix3x3 p, const Matrix3x3 q)
{
	Matrix3x3 out;
	out.e00.re = p.e00.re + q.e00.re;
	out.e00.im = p.e00.im + q.e00.im;
	out.e01.re = p.e01.re + q.e01.re;
	out.e01.im = p.e01.im + q.e01.im;
	out.e02.re = p.e02.re + q.e02.re;
	out.e02.im = p.e02.im + q.e02.im;

	out.e10.re = p.e10.re + q.e10.re;
	out.e10.im = p.e10.im + q.e10.im;
	out.e11.re = p.e11.re + q.e11.re;
	out.e11.im = p.e11.im + q.e11.im;
	out.e12.re = p.e12.re + q.e12.re;
	out.e12.im = p.e12.im + q.e12.im;

	out.e20.re = p.e20.re + q.e20.re;
	out.e20.im = p.e20.im + q.e20.im;
	out.e21.re = p.e21.re + q.e21.re;
	out.e21.im = p.e21.im + q.e21.im;
	out.e22.re = p.e22.re + q.e22.re;
	out.e22.im = p.e22.im + q.e22.im;

	return out;
}

Matrix3x3 subtract_matrix3x3 ( const Matrix3x3 p, const Matrix3x3 q)
{
	Matrix3x3 out;
	out.e00.re = p.e00.re - q.e00.re;
	out.e00.im = p.e00.im - q.e00.im;
	out.e01.re = p.e01.re - q.e01.re;
	out.e01.im = p.e01.im - q.e01.im;
	out.e02.re = p.e02.re - q.e02.re;
	out.e02.im = p.e02.im - q.e02.im;

	out.e10.re = p.e10.re - q.e10.re;
	out.e10.im = p.e10.im - q.e10.im;
	out.e11.re = p.e11.re - q.e11.re;
	out.e11.im = p.e11.im - q.e11.im;
	out.e12.re = p.e12.re - q.e12.re;
	out.e12.im = p.e12.im - q.e12.im;

	out.e20.re = p.e20.re - q.e20.re;
	out.e20.im = p.e20.im - q.e20.im;
	out.e21.re = p.e21.re - q.e21.re;
	out.e21.im = p.e21.im - q.e21.im;
	out.e22.re = p.e22.re - q.e22.re;
	out.e22.im = p.e22.im - q.e22.im;

	return out;
}

Matrix3x3 copy_matrix3x3(const Matrix3x3 p)
{
	Matrix3x3 out;
	out.e00.re = p.e00.re; 
	out.e00.im = p.e00.im; 
	out.e01.re = p.e01.re ;
	out.e01.im = p.e01.im; 
	out.e02.re = p.e02.re; 
	out.e02.im = p.e02.im; 

	out.e10.re = p.e10.re; 
	out.e10.im = p.e10.im; 
	out.e11.re = p.e11.re; 
	out.e11.im = p.e11.im; 
	out.e12.re = p.e12.re; 
	out.e12.im = p.e12.im; 

	out.e20.re = p.e20.re; 
	out.e20.im = p.e20.im; 
	out.e21.re = p.e21.re; 
	out.e21.im = p.e21.im; 
	out.e22.re = p.e22.re; 
	out.e22.im = p.e22.im;
  
	return out;
}

hmc_complex trace_matrix3x3 (const Matrix3x3 p)
{
	hmc_complex out;
	out.re = p.e00.re;
	out.im = p.e00.im;
	out.re += p.e11.re;
	out.im += p.e11.im;
	out.re += p.e22.re;
	out.im += p.e22.im;
  
  return out;
}

Matrix3x3 adjoint_matrix3x3 (const Matrix3x3 p){

	Matrix3x3 out;
        out.e00.re = p.e00.re;
	out.e00.im = - p.e00.im;
        out.e01.re = p.e10.re;
	out.e01.im = - p.e10.im;
        out.e02.re = p.e20.re;
	out.e02.im = - p.e20.im;


	out.e10.re = p.e01.re;
	out.e10.im = - p.e01.im;
	out.e11.re = p.e11.re;
	out.e11.im = - p.e11.im;
	out.e12.re = p.e21.re;
	out.e12.im = - p.e21.im;

	out.e20.re = p.e02.re;
	out.e20.im = - p.e02.im;
	out.e21.re = p.e12.re;
	out.e21.im = - p.e12.im;
	out.e22.re = p.e22.re;
	out.e22.im = - p.e22.im;

	return out;

}

Matrix3x3 matrix_su3to3x3 (const Matrixsu3 p){

	Matrix3x3 out;
	out.e00.re = p.e00.re;
	out.e00.im = p.e00.im;
	out.e01.re = p.e01.re;
	out.e01.im = p.e01.im;
	out.e02.re = p.e02.re;
	out.e02.im = p.e02.im;

	out.e10.re = p.e10.re;
	out.e10.im = p.e10.im;
	out.e11.re = p.e11.re;
	out.e11.im = p.e11.im;
	out.e12.re = p.e12.re;
	out.e12.im = p.e12.im;

#ifdef _RECONSTRUCT_TWELVE_
	out.e20.re = reconstruct_su3(p, 0).re;
	out.e20.im = reconstruct_su3(p, 0).im;
	out.e21.re = reconstruct_su3(p, 1).re;
	out.e21.im = reconstruct_su3(p, 1).im;
	out.e22.re = reconstruct_su3(p, 2).re;
	out.e22.im = reconstruct_su3(p, 2).im;
#else
	out.e20.re = p.e20.re;
	out.e20.im = p.e20.im;
	out.e21.re = p.e21.re;
	out.e21.im = p.e21.im;
	out.e22.re = p.e22.re;
	out.e22.im = p.e22.im;
#endif

	return out;
}


