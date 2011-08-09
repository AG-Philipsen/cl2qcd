/**
 * @file operations used for the molecular dynamics update
 */

//CP: tested version that recreates the matrices made by tmlqcd
/** @todo recheck the factor 0.5 (or F_1_2) that has been deleted here */
/** @todo use of hmc_algebraelement!!!! */

Matrixsu3 build_su3matrix_by_exponentiation(ae inn, hmc_float epsilon)
{

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
