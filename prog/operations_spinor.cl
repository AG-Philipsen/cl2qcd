/** @file
 * Device code for operations on spinors
 */

//opencl_operations_spinor

#ifdef _FERMIONS_
spinor set_spinor_zero()
{
	spinor tmp;
	tmp.e0 = set_su3vec_zero();
	tmp.e1 = set_su3vec_zero();
	tmp.e2 = set_su3vec_zero();
	tmp.e3 = set_su3vec_zero();
	return tmp;
}

hmc_float spinor_squarenorm(spinor in)
{
	hmc_float res=0;
	res += su3vec_squarenorm(in.e0);
	res += su3vec_squarenorm(in.e1);
	res += su3vec_squarenorm(in.e2);
	res += su3vec_squarenorm(in.e3);
	return res;
}

hmc_complex spinor_scalarproduct(spinor in1, spinor in2){
	hmc_complex res = hmc_complex_zero;
	hmc_complex tmp;
	tmp = su3vec_scalarproduct(in1.e0, in2.e0);
	res.re += tmp.re;
	res.im += tmp.im;
	tmp = su3vec_scalarproduct(in1.e1, in2.e1);
	res.re += tmp.re;
	res.im += tmp.im;
	tmp = su3vec_scalarproduct(in1.e2, in2.e2);
	res.re += tmp.re;
	res.im += tmp.im;
	tmp = su3vec_scalarproduct(in1.e3, in2.e3);
	res.re += tmp.re;
	res.im += tmp.im;
	return res;
}

spinor real_multiply_spinor(spinor in, hmc_float factor)
{
	spinor tmp;
	tmp.e0 = su3vec_times_real(in.e0, factor);
	tmp.e1 = su3vec_times_real(in.e1, factor);
	tmp.e2 = su3vec_times_real(in.e2, factor);
	tmp.e3 = su3vec_times_real(in.e3, factor);
	return tmp;
}

spinor spinor_times_complex(spinor in, hmc_complex factor)
{
	spinor tmp;
	tmp.e0 = su3vec_times_complex(in.e0, factor);
	tmp.e1 = su3vec_times_complex(in.e1, factor);
	tmp.e2 = su3vec_times_complex(in.e2, factor);
	tmp.e3 = su3vec_times_complex(in.e3, factor);
	return tmp;
}

spinor spinor_dim(spinor in1, spinor in2)
{
	spinor tmp;
	tmp.e0 = su3vec_dim(in1.e0, in2.e0);
	tmp.e1 = su3vec_dim(in1.e1, in2.e1);
	tmp.e2 = su3vec_dim(in1.e2, in2.e2);
	tmp.e3 = su3vec_dim(in1.e3, in2.e3);
	return tmp;
}

spinor spinor_acc_acc(spinor in1, spinor in2, spinor in3)
{
	spinor tmp;
	tmp.e0 = su3vec_acc_acc(in1.e0, in2.e0, in3.e0);
	tmp.e1 = su3vec_acc_acc(in1.e1, in2.e1, in3.e1);
	tmp.e2 = su3vec_acc_acc(in1.e2, in2.e2, in3.e2);
	tmp.e3 = su3vec_acc_acc(in1.e3, in2.e3, in3.e3);
	return tmp;
}



//CP: I dont think this function is needed anywhere?!?!
// void spinors_accumulate(hmc_spinor* inout, hmc_spinor* incr)
// {
// 	for(int j=0; j<SPINORSIZE; j++) {
// 		inout[j].re += incr[j].re;
// 		inout[j].im += incr[j].im;
// 	}
// 	return;
// }

//CP: I think this can be done by using su3vec_times_complex...
// void spinor_apply_bc(hmc_spinor * in, hmc_float theta)
// {
// 	for(int n = 0; n<SPINORSIZE; n++) {
// 		hmc_float tmp1 = in[n].re;
// 		hmc_float tmp2 = in[n].im;
// 		in[n].re = cos(theta)*tmp1 - sin(theta)*tmp2;
// 		in[n].im = sin(theta)*tmp1 + cos(theta)*tmp2;
// 	}
// 	return;
// }

//spinout =  (1 + 2*i*gamma_5*kappa*mu)spin_in
spinor M_diag_local(spinor y, hmc_float mubar)
{
	spinor out_tmp;

	//CP: how is this called now??
	#ifdef _TWISTEDMASS_
	hmc_complex twistfactor = {1., mubar};
	hmc_complex twistfactor_minus = {1., -1.*mubar};
	//Diagonalpart:
	//	(1+i*mubar*gamma_5)psi = (1, mubar)psi.0,1 (1,-mubar)psi.2,3
        out_tmp.e0 = su3vec_times_complex(y.e0, twistfactor);
	out_tmp.e1 = su3vec_times_complex(y.e1, twistfactor);
	out_tmp.e2 = su3vec_times_complex(y.e2, twistfactor_minus);
	out_tmp.e3 = su3vec_times_complex(y.e3, twistfactor_minus);
	#else
	//Pure Wilson:
	//Diagonalpart:	(1)psi -> spinor not changed!!
	#endif
	return out_tmp;
}


//CP: this is all deprecated!!!
/*
//spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
void dslash_0(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_ocl_su3matrix* u, hmc_ocl_su3matrix* udagger)
{

	spinprojectproduct_gamma0(u,spinnext,-hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma0(udagger,spinprev,hmc_one_f);
	spinors_accumulate(spinout,spinprev);

	return;
}

// spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev
void dslash_1(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_ocl_su3matrix* u, hmc_ocl_su3matrix* udagger)
{

	spinprojectproduct_gamma1(u,spinnext,-hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma1(udagger,spinprev,hmc_one_f);
	spinors_accumulate(spinout,spinprev);

	return;
}

// spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
void dslash_2(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_ocl_su3matrix* u, hmc_ocl_su3matrix* udagger)
{

	spinprojectproduct_gamma2(u,spinnext,-hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma2(udagger,spinprev,hmc_one_f);
	spinors_accumulate(spinout,spinprev);

	return;
}

// spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev
void dslash_3(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_ocl_su3matrix* u, hmc_ocl_su3matrix* udagger)
{

	spinprojectproduct_gamma3(u,spinnext,-hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma3(udagger,spinprev,hmc_one_f);
	spinors_accumulate(spinout,spinprev);

	return;
}
*/
#endif //_FERMIONS_
