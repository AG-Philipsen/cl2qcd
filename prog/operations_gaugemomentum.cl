/**
 * @file operations on algebraelements
 */

ae set_zero_ae(){
	ae tmp;
	tmp.e0 = 0.;
	tmp.e1 = 0.;
	tmp.e2 = 0.;
	tmp.e3 = 0.;
	tmp.e4 = 0.;
	tmp.e5 = 0.;
	tmp.e6 = 0.;
	tmp.e7 = 0.;
	return tmp;
}

hmc_float ae_squarenorm(ae in){
	hmc_float result = 
		(in).e0*(in).e0 +
		(in).e1*(in).e1 +
		(in).e2*(in).e2 +
		(in).e3*(in).e3 +
		(in).e4*(in).e4 +
		(in).e5*(in).e5 +
		(in).e6*(in).e6 +
		(in).e7*(in).e7;
	return result;
}

/** @todo this can propably be defined with a minus...*/
ae acc_factor_times_algebraelement(ae in, hmc_float factor, ae force_in){
	ae tmp;
	tmp.e0 = (in).e0+factor*(force_in).e0; 
	tmp.e1 = (in).e1+factor*(force_in).e1;
	tmp.e2 = (in).e2+factor*(force_in).e2;
	tmp.e3 = (in).e3+factor*(force_in).e3;
	tmp.e4 = (in).e4+factor*(force_in).e4;
	tmp.e5 = (in).e5+factor*(force_in).e5;
	tmp.e6 = (in).e6+factor*(force_in).e6;
	tmp.e7 = (in).e7+factor*(force_in).e7;
	return tmp;
}

//calculates the trace of i times generator times 3x3-matrix and stores this in a su3-algebraelement
//now using structs
ae tr_lambda_u(Matrix3x3 in){
	ae tmp;
	tmp.e0 = ( -in.e10.im - in.e01.im);
	tmp.e1 = (+in.e10.re-in.e01.re);
	tmp.e2 = (-in.e00.im+in.e11.im);
	tmp.e3 = (-in.e20.im-in.e02.im);
	tmp.e4 = (+in.e20.re-in.e02.re);
	tmp.e5 = (-in.e21.im-in.e12.im);
	tmp.e6 = (+in.e21.re-in.e12.re);
	tmp.e7 = (-in.e00.im-in.e11.im + 2.0*in.e22.im)*0.577350269189625;
	return tmp;
}

/**
 * @file operations used by gaugemomentum
 */


/** @todo memcpy ... */
/** @todo get rid of the workaround.. */
// hmc_error copy_gaugemomenta(hmc_gauge_momentum * source, hmc_gauge_momentum * dest){
// 	// copies source to destination within cpu memory, layer for momentum array
// 	return hmc_floatcopy((hmc_float *)source, (hmc_float *)dest, GAUGEMOMENTASIZE); // SL: not tested
// }

//deprecated, this can be done by in = out in the code...
/*
hmc_algebraelement copy_gaugemomenta(hmc_algebraelement2 * source, hmc_algebraelement2 * dest){
	hmc_algebraelement2 tmp;
	for(int i = 0; i<GAUGEMOMENTASIZE; i++){	
		(dest[i]).e0 = (source[i]).e0;
		(dest[i]).e1 = (source[i]).e1;
		(dest[i]).e2 = (source[i]).e2;
		(dest[i]).e3 = (source[i]).e3;
		(dest[i]).e4 = (source[i]).e4;
		(dest[i]).e5 = (source[i]).e5;
		(dest[i]).e6 = (source[i]).e6;
		(dest[i]).e7 = (source[i]).e7;
	}
	return HMC_SUCCESS;
}
*/
	
/** @todo add args for reduction... */
__kernel void gaugemomenta_squarenorm(__global ae * in){
	int id = get_global_id(0);
	if(id == 0){
		hmc_float result = 0.;
		for(int i = 0; i<GAUGEMOMENTASIZE; i++){
			result += ae_squarenorm(in[i]);
		}
		
		/** @todo add reduction.. */
		
	}
}

/** @todo memset... */
__kernel void set_zero_gaugemomentum(__global ae * in){
	int id = get_global_id(0);
	if(id == 0){
		for(int i = 0; i<GAUGEMOMENTASIZE; i++){
			in[i] = set_zero_ae();
		}
	}
}

/** @todo rewrite gaussianComplexVector to handle structs??*/
__kernel void generate_gaussian_gaugemomenta(__global ae * out, __global hmc_ocl_ran * rnd){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	
	hmc_complex tmp;
	
	for(int id_tmp = id; id_tmp < GAUGEMOMENTASIZE; id_tmp += global_size) {	
		//CP: there are 8 elements in ae
		tmp = gaussianNormalPair(&rnd[id]);
		out[id_tmp].e0 = tmp.re;
		out[id_tmp].e1 = tmp.im;
		tmp = gaussianNormalPair(&rnd[id]);
		out[id_tmp].e2 = tmp.re;
		out[id_tmp].e3 = tmp.im;
		tmp = gaussianNormalPair(&rnd[id]);
		out[id_tmp].e4 = tmp.re;
		out[id_tmp].e5 = tmp.im;
		tmp = gaussianNormalPair(&rnd[id]);
		out[id_tmp].e6 = tmp.re;
		out[id_tmp].e7 = tmp.im;
	}
}



void update_gaugemomentum(ae in, hmc_float factor, int global_link_pos, __global ae * out){
	ae tmp = out[global_link_pos];
	acc_factor_times_algebraelement(tmp, factor, in);
	out[global_link_pos] = tmp;
}


//CP: this is the method from tmlqcd
/** @todo decide wheter to use this one or our old one */
Matrixsu3 build_su3matrix_by_exponentiation(ae in, hmc_float epsilon){
	int i;
  Matrixsu3 v,v2,vr;
  hmc_float fac,r;
  hmc_float a,b;
  hmc_complex a0,a1,a2,a1p;

  //make in an su3matrix
	hmc_float halfeps = epsilon;//*F_1_2;
	v.e00.re = 0.0;
	v.e00.im = halfeps*(in.e7*F_1_S3+in.e2);
	v.e01.re = halfeps*in.e1;
	v.e01.im = halfeps*in.e0;
	v.e02.re = halfeps*in.e4;
	v.e02.im = halfeps*in.e3;
	v.e10.re = -halfeps*in.e1;;
	v.e10.im = halfeps*in.e0;;
	v.e11.re = 0.0;
	v.e11.im = halfeps*(in.e7*F_1_S3-in.e2);
	v.e12.re = halfeps*in.e6;
	v.e12.im = halfeps*in.e5;
	v.e20.re = -halfeps*in.e4;
	v.e20.im = halfeps*in.e3;
	v.e21.re = -halfeps*in.e6;
	v.e21.im = halfeps*in.e5;
	v.e22.re = 0.0;
	v.e22.im = -halfeps*2.*in.e7*F_1_S3;
	
  /* calculates v^2 */
  v2 = multiply_matrixsu3 (v,v);
  /* */
  a=0.5*(v2.e00.re+v2.e11.re+v2.e22.re);
  /* 1/3 imaginary part of tr v*v2 */
  b = 0.33333333333333333*
    (v.e00.re*v2.e00.im+v.e00.im*v2.e00.re
     +v.e01.re*v2.e10.im+v.e01.im*v2.e10.re
     +v.e02.re*v2.e20.im+v.e02.im*v2.e20.re
     +v.e10.re*v2.e01.im+v.e10.im*v2.e01.re
     +v.e11.re*v2.e11.im+v.e11.im*v2.e11.re
     +v.e12.re*v2.e21.im+v.e12.im*v2.e21.re
     +v.e20.re*v2.e02.im+v.e20.im*v2.e02.re
     +v.e21.re*v2.e12.im+v.e21.im*v2.e12.re
     +v.e22.re*v2.e22.im+v.e22.im*v2.e22.re  );
  a0.re=0.16059043836821615e-9;    /*  1/13! */
  a0.im=0.0;
  a1.re=0.11470745597729725e-10;   /*  1/14! */
  a1.im=0.0;
  a2.re=0.76471637318198165e-12;   /*  1/15! */
  a2.im=0.0;
  fac=0.20876756987868099e-8;      /*  1/12! */
  r=12.0;
  for(i = 3; i <= 15; i++) {
    a1p.re = a0.re + a * a2.re;
    a1p.im = a0.im + a * a2.im;
    a0.re = fac - b * a2.im;
    a0.im =     + b * a2.re;
    a2.re = a1.re; 
    a2.im = a1.im;
    a1.re = a1p.re; 
    a1.im = a1p.im;
    fac *= r;  
    r -= 1.0;
  }
  /* vr = a0 + a1*v + a2*v2 */
  vr.e00.re = a0.re + a1.re*v.e00.re - a1.im*v.e00.im + a2.re*v2.e00.re - a2.im*v2.e00.im;
  vr.e00.im = a0.im + a1.re*v.e00.im + a1.im*v.e00.re + a2.re*v2.e00.im + a2.im*v2.e00.re;
  vr.e01.re =         a1.re*v.e01.re - a1.im*v.e01.im + a2.re*v2.e01.re - a2.im*v2.e01.im;
  vr.e01.im =         a1.re*v.e01.im + a1.im*v.e01.re + a2.re*v2.e01.im + a2.im*v2.e01.re;
  vr.e02.re =         a1.re*v.e02.re - a1.im*v.e02.im + a2.re*v2.e02.re - a2.im*v2.e02.im;
  vr.e02.im =         a1.re*v.e02.im + a1.im*v.e02.re + a2.re*v2.e02.im + a2.im*v2.e02.re;
  vr.e10.re =         a1.re*v.e10.re - a1.im*v.e10.im + a2.re*v2.e10.re - a2.im*v2.e10.im;
  vr.e10.im =         a1.re*v.e10.im + a1.im*v.e10.re + a2.re*v2.e10.im + a2.im*v2.e10.re;
  vr.e11.re = a0.re + a1.re*v.e11.re - a1.im*v.e11.im + a2.re*v2.e11.re - a2.im*v2.e11.im;
  vr.e11.im = a0.im + a1.re*v.e11.im + a1.im*v.e11.re + a2.re*v2.e11.im + a2.im*v2.e11.re;
  vr.e12.re =         a1.re*v.e12.re - a1.im*v.e12.im + a2.re*v2.e12.re - a2.im*v2.e12.im;
  vr.e12.im =         a1.re*v.e12.im + a1.im*v.e12.re + a2.re*v2.e12.im + a2.im*v2.e12.re;
  vr.e20.re =         a1.re*v.e20.re - a1.im*v.e20.im + a2.re*v2.e20.re - a2.im*v2.e20.im;
  vr.e20.im =         a1.re*v.e20.im + a1.im*v.e20.re + a2.re*v2.e20.im + a2.im*v2.e20.re;
  vr.e21.re =         a1.re*v.e21.re - a1.im*v.e21.im + a2.re*v2.e21.re - a2.im*v2.e21.im;
  vr.e21.im =         a1.re*v.e21.im + a1.im*v.e21.re + a2.re*v2.e21.im + a2.im*v2.e21.re;
  vr.e22.re = a0.re + a1.re*v.e22.re - a1.im*v.e22.im + a2.re*v2.e22.re - a2.im*v2.e22.im;
  vr.e22.im = a0.im + a1.re*v.e22.im + a1.im*v.e22.re + a2.re*v2.e22.im + a2.im*v2.e22.re;
  return vr;
}

//CP: tested version that recreates the matrices made by tmlqcd
/** @todo recheck the factor 0.5 (or F_1_2) that has been deleted here */
/** @todo use of hmc_algebraelement!!!! */
/*
typedef hmc_3x3matrix hmc_complex[3][3]
Matrix_su3 build_su3matrix_by_exponentiation(ae inn, hmc_float epsilon){
		Matrixsu3 tmp;
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
	
	return tmp;
}
*/
