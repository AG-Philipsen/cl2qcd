/** @file
 * Device code implementing SU(3)matrices with and without reconstruct 12
 */
//operations_matrix_su3.cl



//not valid for REC12
void print_matrixsu3(Matrixsu3 in){
     printf("(%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f)\n", in.e00.re, in.e00.im, in.e01.re, in.e01.im, in.e02.re, in.e02.im, in.e10.re, in.e10.im, in.e11.re, in.e11.im, in.e12.re, in.e12.im, in.e02.re, in.e02.im, in.e21.re, in.e21.im, in.e22.re, in.e22.im);
     printf("\n");
}


#ifdef _RECONSTRUCT_TWELVE_
hmc_complex reconstruct_su3(const Matrixsu3 p, const int ncomp)
{
	hmc_complex out;
	hmc_complex tmp;
	
	switch (ncomp){

	case 0:
	    tmp = complexmult (p.e01, p.e12);
	    out = complexmult (p.e02, p.e11);
	    break;

	case 1:
	    tmp = complexmult (p.e02, p.e10);
	    out = complexmult (p.e00, p.e12);
	    break;
	case 2:
	    tmp = complexmult (p.e00, p.e11);
	    out = complexmult (p.e01, p.e10);
	    break;
	}

	out = complexsubtract (tmp, out);
	return complexconj (out);
}
#endif

Matrixsu3 get_matrixsu3( __global ocl_s_gaugefield * field, const int spacepos, const int timepos, const int mu)
{
	Matrixsu3  out;
	out = field [get_global_link_pos(mu, spacepos, timepos)];


	return out;
}

void put_matrixsu3(__global ocl_s_gaugefield * field, const Matrixsu3 in, const int spacepos, const int timepos, const int mu)
{
	    field [get_global_link_pos(mu, spacepos, timepos)] = in;
}

Matrixsu3 copy_matrixsu3(const Matrixsu3 in)
{
	Matrixsu3 out;

	out.e00.re = in.e00.re;
	out.e00.im = in.e00.im;
	out.e01.re = in.e01.re;
	out.e01.im = in.e01.im;
	out.e02.re = in.e02.re;
	out.e02.im = in.e02.im;

	out.e10.re = in.e10.re;
	out.e10.im = in.e10.im;
	out.e11.re = in.e11.re;
	out.e11.im = in.e11.im;
	out.e12.re = in.e12.re;
	out.e12.im = in.e12.im;

#ifndef _RECONSTRUCT_TWELVE_
	out.e20.re = in.e20.re;
	out.e20.im = in.e20.im;
	out.e21.re = in.e21.re;
	out.e21.im = in.e21.im;
	out.e22.re = in.e22.re;
	out.e22.im = in.e22.im;
#endif
	return out;
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


Matrixsu3 zero_matrixsu3()
{
	Matrixsu3 out;
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

#ifndef _RECONSTRUCT_TWELVE_
	out.e20.re = 0.;
	out.e20.im = 0.;
	out.e21.re = 0.;
	out.e21.im = 0.;
	out.e22.re = 0.;
	out.e22.im = 0.;
#endif 

	return out;
}

Matrixsu3 multiply_matrixsu3(const Matrixsu3 p, const Matrixsu3 q)
{
    Matrixsu3 out;
#ifdef _RECONSTRUCT_TWELVE_
    hmc_float p20_re = reconstruct_su3(p, 0).re;
    hmc_float p20_im = reconstruct_su3(p, 0).im;
    hmc_float p21_re = reconstruct_su3(p, 1).re;
    hmc_float p21_im = reconstruct_su3(p, 1).im;
    hmc_float p22_re = reconstruct_su3(p, 2).re;
    hmc_float p22_im = reconstruct_su3(p, 2).im;

    hmc_float q20_re = reconstruct_su3(q, 0).re;
    hmc_float q20_im = reconstruct_su3(q, 0).im;
    hmc_float q21_re = reconstruct_su3(q, 1).re;
    hmc_float q21_im = reconstruct_su3(q, 1).im;
    hmc_float q22_re = reconstruct_su3(q, 2).re;
    hmc_float q22_im = reconstruct_su3(q, 2).im;

    out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e10.re + p.e02.re * q20_re
							 - p.e00.im * q.e00.im - p.e01.im * q.e10.im - p.e02.im * q20_im;
    out.e00.im = p.e00.re * q.e00.im + p.e01.re * q.e10.im + p.e02.re * q20_im
							 + p.e00.im * q.e00.re + p.e01.im * q.e10.re + p.e02.im * q20_re;
	      
    out.e01.re = p.e00.re * q.e01.re + p.e01.re * q.e11.re + p.e02.re * q21_re
							 - p.e00.im * q.e01.im - p.e01.im * q.e11.im - p.e02.im * q21_im;
    out.e01.im = p.e00.re * q.e01.im + p.e01.re * q.e11.im + p.e02.re * q21_im
	       			 + p.e00.im * q.e01.re + p.e01.im * q.e11.re + p.e02.im * q21_re;
	       
    out.e02.re = p.e00.re * q.e02.re + p.e01.re * q.e12.re + p.e02.re * q22_re
	       			 - p.e00.im * q.e02.im - p.e01.im * q.e12.im - p.e02.im * q22_im;
    out.e02.im = p.e00.re * q.e02.im + p.e01.re * q.e12.im + p.e02.re * q22_im
	       			 + p.e00.im * q.e02.re + p.e01.im * q.e12.re + p.e02.im * q22_re;
	       
    out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e10.re + p.e12.re * q20_re
	       			 - p.e10.im * q.e00.im - p.e11.im * q.e10.im - p.e12.im * q20_im;
    out.e10.im = p.e10.re * q.e00.im + p.e11.re * q.e10.im + p.e12.re * q20_im
	       			 + p.e10.im * q.e00.re + p.e11.im * q.e10.re + p.e12.im * q20_re;
	       
    out.e11.re = p.e10.re * q.e01.re + p.e11.re * q.e11.re + p.e12.re * q21_re
	       			 - p.e10.im * q.e01.im - p.e11.im * q.e11.im - p.e12.im * q21_im;
    out.e11.im = p.e10.re * q.e01.im + p.e11.re * q.e11.im + p.e12.re * q21_im
	       			 + p.e10.im * q.e01.re + p.e11.im * q.e11.re + p.e12.im * q21_re;
	       
    out.e12.re = p.e10.re * q.e02.re + p.e11.re * q.e12.re + p.e12.re * q22_re
	       			 - p.e10.im * q.e02.im - p.e11.im * q.e12.im - p.e12.im * q22_im;
    out.e12.im = p.e10.re * q.e02.im + p.e11.re * q.e12.im + p.e12.re * q22_im
	       			 + p.e10.im * q.e02.re + p.e11.im * q.e12.re + p.e12.im * q22_re;	       
#else
    
    out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e10.re + p.e02.re * q.e20.re
	       			 - p.e00.im * q.e00.im - p.e01.im * q.e10.im - p.e02.im * q.e20.im;
    out.e00.im = p.e00.re * q.e00.im + p.e01.re * q.e10.im + p.e02.re * q.e20.im
	       			 + p.e00.im * q.e00.re + p.e01.im * q.e10.re + p.e02.im * q.e20.re;

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

#endif
    
    return out;
}

Matrixsu3 multiply_matrixsu3_dagger(const Matrixsu3 p, const Matrixsu3 q)
{
    Matrixsu3 out;
#ifdef _RECONSTRUCT_TWELVE_
    hmc_float p20_re = reconstruct_su3(p, 0).re;
    hmc_float p20_im = reconstruct_su3(p, 0).im;
    hmc_float p21_re = reconstruct_su3(p, 1).re;
    hmc_float p21_im = reconstruct_su3(p, 1).im;
    hmc_float p22_re = reconstruct_su3(p, 2).re;
    hmc_float p22_im = reconstruct_su3(p, 2).im;

    hmc_float q20_re = reconstruct_su3(q, 0).re;
    hmc_float q20_im = reconstruct_su3(q, 0).im;
    hmc_float q21_re = reconstruct_su3(q, 1).re;
    hmc_float q21_im = reconstruct_su3(q, 1).im;
    hmc_float q22_re = reconstruct_su3(q, 2).re;
    hmc_float q22_im = reconstruct_su3(q, 2).im;

		out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e01.re + p.e02.re * q02_re
		           + p.e00.im * q.e00.im + p.e01.im * q.e01.im + p.e02.im * q02_im;
		out.e00.im =-p.e00.re * q.e00.im - p.e01.re * q.e01.im - p.e02.re * q02_im
		           + p.e00.im * q.e00.re + p.e01.im * q.e01.re + p.e02.im * q02_re;
	      
		out.e01.re = p.e00.re * q.e10.re + p.e01.re * q.e11.re + p.e02.re * q12_re
		           + p.e00.im * q.e10.im + p.e01.im * q.e11.im + p.e02.im * q12_im;
		out.e01.im =-p.e00.re * q.e10.im - p.e01.re * q.e11.im - p.e02.re * q12_im
		           + p.e00.im * q.e10.re + p.e01.im * q.e11.re + p.e02.im * q12_re;
	       
		out.e02.re = p.e00.re * q.e20.re + p.e01.re * q.e21.re + p.e02.re * q22_re
		           + p.e00.im * q.e20.im + p.e01.im * q.e21.im + p.e02.im * q22_im;
		out.e02.im =-p.e00.re * q.e20.im - p.e01.re * q.e21.im - p.e02.re * q22_im
		           + p.e00.im * q.e20.re + p.e01.im * q.e21.re + p.e02.im * q22_re;
	       
		out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e01.re + p.e12.re * q02_re
		           + p.e10.im * q.e00.im + p.e11.im * q.e01.im + p.e12.im * q02_im;
		out.e10.im =-p.e10.re * q.e00.im - p.e11.re * q.e01.im - p.e12.re * q02_im
		           + p.e10.im * q.e00.re + p.e11.im * q.e01.re + p.e12.im * q02_re;
	       
		out.e11.re = p.e10.re * q.e10.re + p.e11.re * q.e11.re + p.e12.re * q12_re
		           + p.e10.im * q.e10.im + p.e11.im * q.e11.im + p.e12.im * q12_im;
		out.e11.im =-p.e10.re * q.e10.im - p.e11.re * q.e11.im - p.e12.re * q12_im
		           + p.e10.im * q.e10.re + p.e11.im * q.e11.re + p.e12.im * q12_re;
	       
		out.e12.re = p.e10.re * q.e20.re + p.e11.re * q.e21.re + p.e12.re * q22_re
		           + p.e10.im * q.e20.im + p.e11.im * q.e21.im + p.e12.im * q22_im;
		out.e12.im =-p.e10.re * q.e20.im - p.e11.re * q.e21.im - p.e12.re * q22_im
		           + p.e10.im * q.e20.re + p.e11.im * q.e21.re + p.e12.im * q22_re;	       
#else
    
		out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e01.re + p.e02.re * q.e02.re
		           + p.e00.im * q.e00.im + p.e01.im * q.e01.im + p.e02.im * q.e02.im;
		out.e00.im =-p.e00.re * q.e00.im - p.e01.re * q.e01.im - p.e02.re * q.e02.im
		           + p.e00.im * q.e00.re + p.e01.im * q.e01.re + p.e02.im * q.e02.re;

		out.e01.re = p.e00.re * q.e10.re + p.e01.re * q.e11.re + p.e02.re * q.e12.re
		           + p.e00.im * q.e10.im + p.e01.im * q.e11.im + p.e02.im * q.e12.im;
		out.e01.im =-p.e00.re * q.e10.im - p.e01.re * q.e11.im - p.e02.re * q.e12.im
		           + p.e00.im * q.e10.re + p.e01.im * q.e11.re + p.e02.im * q.e12.re;
    
		out.e02.re = p.e00.re * q.e20.re + p.e01.re * q.e21.re + p.e02.re * q.e22.re
		           + p.e00.im * q.e20.im + p.e01.im * q.e21.im + p.e02.im * q.e22.im;
		out.e02.im =-p.e00.re * q.e20.im - p.e01.re * q.e21.im - p.e02.re * q.e22.im
		           + p.e00.im * q.e20.re + p.e01.im * q.e21.re + p.e02.im * q.e22.re;
	       
		out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e01.re + p.e12.re * q.e02.re
		           + p.e10.im * q.e00.im + p.e11.im * q.e01.im + p.e12.im * q.e02.im;
		out.e10.im =-p.e10.re * q.e00.im - p.e11.re * q.e01.im - p.e12.re * q.e02.im
		           + p.e10.im * q.e00.re + p.e11.im * q.e01.re + p.e12.im * q.e02.re;
	       
		out.e11.re = p.e10.re * q.e10.re + p.e11.re * q.e11.re + p.e12.re * q.e12.re
		           + p.e10.im * q.e10.im + p.e11.im * q.e11.im + p.e12.im * q.e12.im;
		out.e11.im =-p.e10.re * q.e10.im - p.e11.re * q.e11.im - p.e12.re * q.e12.im
		           + p.e10.im * q.e10.re + p.e11.im * q.e11.re + p.e12.im * q.e12.re;
	       
		out.e12.re = p.e10.re * q.e20.re + p.e11.re * q.e21.re + p.e12.re * q.e22.re
		           + p.e10.im * q.e20.im + p.e11.im * q.e21.im + p.e12.im * q.e22.im;
		out.e12.im =-p.e10.re * q.e20.im - p.e11.re * q.e21.im - p.e12.re * q.e22.im
		           + p.e10.im * q.e20.re + p.e11.im * q.e21.re + p.e12.im * q.e22.re;	       

		out.e20.re = p.e20.re * q.e00.re + p.e21.re * q.e01.re + p.e22.re * q.e02.re
		           + p.e20.im * q.e00.im + p.e21.im * q.e01.im + p.e22.im * q.e02.im;
		out.e20.im =-p.e20.re * q.e00.im - p.e21.re * q.e01.im - p.e22.re * q.e02.im
		           + p.e20.im * q.e00.re + p.e21.im * q.e01.re + p.e22.im * q.e02.re;
	       
		out.e21.re = p.e20.re * q.e10.re + p.e21.re * q.e11.re + p.e22.re * q.e12.re
		           + p.e20.im * q.e10.im + p.e21.im * q.e11.im + p.e22.im * q.e12.im;
		out.e21.im =-p.e20.re * q.e10.im - p.e21.re * q.e11.im - p.e22.re * q.e12.im
		           + p.e20.im * q.e10.re + p.e21.im * q.e11.re + p.e22.im * q.e12.re;
	       
		out.e22.re = p.e20.re * q.e20.re + p.e21.re * q.e21.re + p.e22.re * q.e22.re
		           + p.e20.im * q.e20.im + p.e21.im * q.e21.im + p.e22.im * q.e22.im;
		out.e22.im =-p.e20.re * q.e20.im - p.e21.re * q.e21.im - p.e22.re * q.e22.im
		           + p.e20.im * q.e20.re + p.e21.im * q.e21.re + p.e22.im * q.e22.re;
    
    return out;
}

Matrixsu3 adjoint_matrixsu3(const Matrixsu3 p)
{
	Matrixsu3 out;
	out.e00.re = p.e00.re;
	out.e00.im = - p.e00.im;
        out.e01.re = p.e10.re;
	out.e01.im = - p.e10.im;
        
	out.e10.re = p.e01.re;
	out.e10.im = - p.e01.im;
	out.e11.re = p.e11.re;
	out.e11.im = - p.e11.im;
	

#ifdef _RECONSTRUCT_TWELVE_
	out.e02.re = reconstruct_su3(p, 0).re;
	out.e02.im = - reconstruct_su3(p, 0).im;
	
	out.e12.re = reconstruct_su3(p, 1).re;
	out.e12.im = - reconstruct_su3(p, 1).im;
#else
	out.e02.re = p.e20.re;
	out.e02.im = - p.e20.im;
	
	out.e12.re = p.e21.re;
	out.e12.im = - p.e21.im;

	out.e20.re = p.e02.re;
	out.e20.im = - p.e02.im;
	out.e21.re = p.e12.re;
	out.e21.im = - p.e12.im;
	out.e22.re = p.e22.re;
	out.e22.im = - p.e22.im;
#endif
	return out;
}

hmc_complex trace_matrixsu3(const Matrixsu3 p)
{
	hmc_complex out;
	out.re = p.e00.re;
	out.im = p.e00.im;
	out.re += p.e11.re;
	out.im += p.e11.im;
#ifdef _RECONSTRUCT_TWELVE_
	out.re += reconstruct_su3(p, 2).re;
	out.re += reconstruct_su3(p, 2).im;
#else
	out.re += p.e22.re;
	out.im += p.e22.im;
#endif
	return out;
}

hmc_complex det_matrixsu3(const Matrixsu3 p)
{

	hmc_complex out;

#ifdef _RECONSTRUCT_TWELVE_
	hmc_complex c0 = reconstruct_su3(p, 0);
	hmc_complex c1 = reconstruct_su3(p, 1);
	hmc_complex c2 = reconstruct_su3(p, 2);


	hmc_complex det1, det2, det3, det4, det5, det6, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
	tmp1 = complexmult( p.e11, c2);
	det1 = complexmult( p.e00, tmp1);
	tmp2 = complexmult( p.e12, c0);
	det2 = complexmult( p.e01, tmp2);
	tmp3 = complexmult( p.e10, c1);
	det3 = complexmult( p.e02, tmp3);
	tmp4 = complexmult( p.e11, c0 );
	det4 = complexmult( p.e02, tmp4);
	tmp5 = complexmult( p.e10, c2 );
	det5 = complexmult( p.e01, tmp5);
	tmp6 = complexmult( p.e12, c1 );
	det6 = complexmult( p.e00, tmp6);
#else
	hmc_complex det1, det2, det3, det4, det5, det6, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
	tmp1 = complexmult( p.e11, p.e22);
	det1 = complexmult( p.e00, tmp1);
	tmp2 = complexmult( p.e12, p.e20);
	det2 = complexmult( p.e01, tmp2);
	tmp3 = complexmult( p.e10, p.e21);
	det3 = complexmult( p.e02, tmp3);
	tmp4 = complexmult( p.e11, p.e20 );
	det4 = complexmult( p.e02, tmp4);
	tmp5 = complexmult( p.e10, p.e22 );
	det5 = complexmult( p.e01, tmp5);
	tmp6 = complexmult( p.e12, p.e21 );
	det6 = complexmult( p.e00, tmp6);
#endif
	out.re = det1.re + det2.re + det3.re - det4.re - det5.re - det6.re;
	out.im = det1.im + det2.im + det3.im - det4.im - det5.im - det6.im;

	return out;
}

