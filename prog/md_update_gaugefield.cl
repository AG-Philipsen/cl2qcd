/**
 * @file operations used for the molecular dynamics update
 */

//CP: tested version that recreates the matrices made by tmlqcd
/** @todo recheck the factor 0.5 (or F_1_2) that has been deleted here */

//this build a su3-matrix from an algebraelement!!
Matrixsu3 build_su3_from_ae(ae in){
	Matrixsu3 v;
	
	v.e00.re = 0.0;
	v.e00.im = (in.e7 * F_1_S3 + in.e2);
	v.e01.re = in.e1;
	v.e01.im = in.e0;
	v.e02.re = in.e4;
	v.e02.im = in.e3;
	v.e10.re = -in.e1;;
	v.e10.im = in.e0;;
	v.e11.re = 0.0;
	v.e11.im = (in.e7 * F_1_S3 - in.e2);
	v.e12.re = in.e6;
	v.e12.im = in.e5;
	v.e20.re = -in.e4;
	v.e20.im = in.e3;
	v.e21.re = -in.e6;
	v.e21.im = in.e5;
	v.e22.re = 0.0;
	v.e22.im = -2.*in.e7 * F_1_S3;
	
	return v;
}

//this build a su3-matrix from an algebraelement and in addition multiplies it by a real number!!
Matrixsu3 build_su3_from_ae_times_real(ae in, hmc_float eps){
	Matrixsu3 v;
	
	v.e00.re = 0.0;
	v.e00.im = eps * (in.e7 * F_1_S3 + in.e2);
	v.e01.re = eps * in.e1;
	v.e01.im = eps * in.e0;
	v.e02.re = eps * in.e4;
	v.e02.im = eps * in.e3;
	v.e10.re = -eps * in.e1;;
	v.e10.im = eps * in.e0;;
	v.e11.re = 0.0;
	v.e11.im = eps * (in.e7 * F_1_S3 - in.e2);
	v.e12.re = eps * in.e6;
	v.e12.im = eps * in.e5;
	v.e20.re = -eps * in.e4;
	v.e20.im = eps * in.e3;
	v.e21.re = -eps * in.e6;
	v.e21.im = eps * in.e5;
	v.e22.re = 0.0;
	v.e22.im = -eps * 2.*in.e7 * F_1_S3;
	
	return v;
}


Matrixsu3 build_su3matrix_by_exponentiation(ae inn, hmc_float epsilon)
{
	//CP: this is the old version
/*
	Matrixsu3 tmp;
	hmc_float in [8];
	in[0] = inn.e0;
	in[1] = inn.e1;
	in[2] = inn.e2;
	in[3] = inn.e3;
	in[4] = inn.e4;
	in[5] = inn.e5;
	in[6] = inn.e6;
	in[7] = inn.e7;

	//this is the case where one actually evaluates -as matrices- many orders of exp(i*e*P)=1+i*e*P+(1/2)(i*e*P)^2 + ...
	//1. create mat = (i epsilon)/2 p_i lambda_i    with lambda=gellmann matrix
	Matrix3x3 eMat;
	hmc_float halfeps = epsilon;//*F_1_2;
	eMat.e00.re = 0.0;
	eMat.e00.im = halfeps * (in[7] * F_1_S3 + in[2]);
	eMat.e01.re = halfeps * in[1];
	eMat.e01.im = halfeps * in[0];
	eMat.e02.re = halfeps * in[4];
	eMat.e02.im = halfeps * in[3];
	eMat.e10.re = -halfeps * in[1];;
	eMat.e10.im = halfeps * in[0];;
	eMat.e11.re = 0.0;
	eMat.e11.im = halfeps * (in[7] * F_1_S3 - in[2]);
	eMat.e12.re = halfeps * in[6];
	eMat.e12.im = halfeps * in[5];
	eMat.e20.re = -halfeps * in[4];
	eMat.e20.im = halfeps * in[3];
	eMat.e21.re = -halfeps * in[6];
	eMat.e21.im = halfeps * in[5];
	eMat.e22.re = 0.0;
	eMat.e22.im = -halfeps * 2.*in[7] * F_1_S3;

	//2. start with the exp(...) expansion by using standard 3x3 sum and multiplication
	Matrix3x3 eRes, eCurPower, eNextPower, eLastResult;
	hmc_float eAccuracyCheck;
	eRes = identity_matrix3x3();
	eCurPower = identity_matrix3x3();
	hmc_float eCoefficient = 1.0;
	int factorial = 1;
	for(int power = 1; power < 50; power++) {
		eNextPower = multiply_matrix3x3(eMat, eCurPower);
		eCurPower = eNextPower;
		//CP: the calc of factorial can be simplified by just doing factorial*=power!!
		for(int i = power; i > 1; i--) factorial *= i;
		eCurPower = multiply_matrix3x3_by_real(eCurPower, (1. / (1.*factorial)));
		factorial = 1;
		eLastResult = eRes;
		eRes = add_matrix3x3(eRes, eCurPower);
		eAccuracyCheck = absoluteDifference_matrix3x3(eRes, eLastResult);
		if(eAccuracyCheck < 1e-15) {
			break;
		}
	}
	//3. here I have the exponentiated matrix in 3x3 generic form (eRes), project it
#ifdef _RECONSTRUCT_TWELVE_
	tmp.e0 = eRes.e00;
	tmp.e1 = eRes.e01;
	tmp.e2 = eRes.e02;
	tmp.e3 = eRes.e10;
	tmp.e4 = eRes.e11;
	tmp.e5 = eRes.e12;
#else
	tmp.e00 = eRes.e00;
	tmp.e01 = eRes.e01;
	tmp.e02 = eRes.e02;
	tmp.e10 = eRes.e10;
	tmp.e11 = eRes.e11;
	tmp.e12 = eRes.e12;
	tmp.e20 = eRes.e20;
	tmp.e21 = eRes.e21;
	tmp.e22 = eRes.e22;
#endif // _RECONSTRUCT_TWELVE_
	project_su3(tmp);

	return tmp;
	*/
	
	//this is the method taken from tmqlcd. It is a cut-off series.
	int i;
  Matrixsu3 v,v2,vr;
  hmc_float fac,r;
  hmc_float a,b;
  hmc_complex a0,a1,a2,a1p;

  //original in tmlqcd: _make_su3(v,p);
	//Here, the stepsize factor is multiplied in right away!!
	//Also, a factor of 0.5 is taken out to fit the different trace-definition from tmqlcd
	hmc_float halfeps = epsilon;// *F_1_2;
	v = build_su3_from_ae_times_real(inn, halfeps);
	
  // calculates v^2 
  v2 = multiply_matrixsu3(v,v);
  a=0.5*(v2.e00.re+v2.e11.re+v2.e22.re);
  // 1/3 imaginary part of tr v*v2 
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
  a0.re=0.16059043836821615e-9;    //  1/13! 
  a0.im=0.0;
  a1.re=0.11470745597729725e-10;   //  1/14! 
  a1.im=0.0;
  a2.re=0.76471637318198165e-12;   //  1/15! 
  a2.im=0.0;
  fac=0.20876756987868099e-8;      //  1/12! 
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
  // vr = a0 + a1*v + a2*v2 
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

//molecular dynamics update for the gaugefield:
//u_out = exp(i eps p_in) u_in
__kernel void md_update_gaugefield(hmc_float eps, __global ae * p_in, __global ocl_s_gaugefield * u_inout)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int n, t, index;
	Matrixsu3 tmp;
	Matrixsu3 tmp2;

	//CP: it is GAUGEMOMENTASIZE = NDIM * SPINORFIELDSIZE
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {
		/** @todo this must be done more efficient */
		if(id_tmp < VOLSPACE * NTIME / 2)
			get_even_site(id_tmp, &n, &t);
		else
			get_odd_site(id_tmp-(VOLSPACE*NTIME/2), &n, &t);
		
		for(int mu = 0; mu < NDIM; mu++) {
			index = get_global_link_pos(mu, n, t);
			// an su3 algebra element has NC*NC-1 = 8 hmc_float entries
			// &(p_in[index*8]) should point to the right position for the pos-th element of the long gaugemomentum vector p_in
			tmp2 = build_su3matrix_by_exponentiation((p_in[index]), eps);

			tmp = get_matrixsu3(u_inout, n, t, mu);
			tmp2 = multiply_matrixsu3( tmp2, tmp);
			put_matrixsu3(u_inout, tmp2, n, t, mu);
		}
	}
}
