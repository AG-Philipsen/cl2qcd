/** @file
 * Device code for operations on the fermion matrix
 */

//operations_gaugefield.cl


// hmc_complex global_trace_su3(__global hmc_ocl_gaugefield * field, int mu)
// {
// 	hmc_complex sum;
// 	sum.re = 0;
// 	sum.im = 0;
// 	for(int t=0; t<NTIME; t++) {
// 		for(int n=0; n<VOLSPACE; n++) {
// 			hmc_ocl_su3matrix tmp[SU3SIZE];
// 			get_su3matrix(tmp, field, n, t, mu);
// 			sum.re += trace_su3matrix(tmp).re;
// 			sum.im += trace_su3matrix(tmp).im;
// 		}
// 	}
// 	return sum;
// }

Matrixsu3 project_su3(const Matrixsu3 U){

  Matrixsu3 out;
  
  //Extract initial vectors
  hmc_complex a[NC];
  hmc_complex b[NC];
  
#ifdef _RECONSTRUCT_TWELVE_
  a[0] = U.e00;
  a[1] = U.e01;
  a[2] = U.e02;
  b[0] = U.e10;
  b[1] = U.e11;
  b[2] = U.e12;
#else
  hmc_complex c[NC];
  a[0] = U.e00;
  a[1] = U.e01;
  a[2] = U.e02;
  b[0] = U.e10;
  b[1] = U.e11;
  b[2] = U.e12;
  c[0] = U.e20;
  c[1] = U.e21;
  c[2] = U.e22;
#endif
  
  //New SU3-Matrix
  //first vector
  //norm
  hmc_float norm = 0.;
  for (int i=0; i<NC; i++){
    hmc_complex tmp = complexconj(a[i]);
    tmp = complexmult (a[i], tmp);
    norm += tmp.re;
  }
  norm = 1./sqrt(norm);
  //rescale
  for (int i=0; i<NC; i++){
    a[i].re *= norm;
    a[i].im *= norm;
  }
  
  //second vector
  //orthogonal vector
  hmc_complex factor;
  factor.re = 0.0;
  factor.im = 0.0;
  for (int i=0; i<NC; i++){
    hmc_complex tmp;
    tmp = complexconj (b[i]);
    tmp = complexmult (a[i], tmp);
    factor = complexadd (factor, tmp);
  }
  for (int i=0; i<NC; i++){
    hmc_complex tmp;
    tmp = complexmult(factor, a[i]);
    b[i] = complexsubtract(b[i], tmp); 
  }
  
//norm
  norm = 0.;
  for (int i=0; i<NC; i++)
  {
    hmc_complex tmp;
    tmp = complexconj(b[i]);
    tmp = complexmult (b[i], tmp);
    norm +=  tmp.re;
  }
  norm = 1./sqrt(norm);
  //rescale
  for  (int i=0; i<NC; i++){
    b[i].re *= norm;
    b[i].im *= norm;
  }

#ifndef _RECONSTRUCT_TWELVE_
  //third vector 
  //orthogonal vector
  hmc_complex tmp;
  hmc_complex tmp2;
  tmp = complexmult(a[1], b[2]);
  tmp = complexconj(tmp);
  tmp2 = complexmult(a[2], b[1]);
  tmp2 = complexconj(tmp2);
  c[0] = complexsubtract(tmp, tmp2);
  tmp = complexmult(a[2], b[0]);
  tmp = complexconj(tmp);
  tmp2 = complexmult(a[0], b[2]);
  tmp2 = complexconj(tmp2);
  c[1] = complexsubtract(tmp, tmp2);
  tmp = complexmult(a[0], b[1]);
  tmp = complexconj(tmp);
  tmp2 = complexmult(a[1], b[0]);
  tmp2 = complexconj(tmp2);
  c[2] = complexsubtract(tmp, tmp2);
   //Set new values to matrix
  out.e20 = c[0];
  out.e21 = c[1];
  out.e22 = c[2];
#endif

   //Set new values to matrix
  out.e00 = a[0];
  out.e10 = b[0];
  out.e01 = a[1];
  out.e11 = b[1];
  out.e02 = a[2];
  out.e12 = b[2];
  
  return out;
}

Matrixsu2 reduction (const Matrix3x3 src, const int rand)
{
	Matrixsu2 out;
	if(rand == 1) {
		out.e00 = src.e00;
		out.e01 = src.e01;
		out.e10 = src.e10;
		out.e11 = src.e11;
	} else if (rand==2) {
		out.e00 = src.e11;
		out.e01 = src.e12;
		out.e10 = src.e21;
		out.e11 = src.e22;
	} else if (rand==3) {
		out.e00 = src.e00;
		out.e01 = src.e02;
		out.e10 = src.e20;
		out.e11 = src.e22;
	}
	return out;
}

//CP: I leave both version in here...
Matrixsu3 extend (const int random, Matrixsu2 src)
{
        Matrixsu3 out;
#ifdef _USE_GPU_

#ifdef _RECONSTRUCT_TWELVE_
switch(random){
        case 1:
                out.e00 = src.e00;
                out.e01 = src.e01;
                out.e02 = hmc_complex_zero;
                out.e10 = src.e10;
                out.e11 = src.e11;
                out.e12 = hmc_complex_zero;
                return out;
                break;
        case 2:
                out.e00 = hmc_complex_one;
                out.e01 = hmc_complex_zero;
                out.e02 = hmc_complex_zero;
                out.e10 = hmc_complex_zero;
                out.e11 = src.e00;
                out.e12 = src.e01;
                return out;
                break;
        case 3:
                out.e00 = src.e00;
                out.e01 = hmc_complex_zero;
                out.e02 = src.e01;
                out.e10= hmc_complex_zero;
                out.e11 = hmc_complex_one;
                out.e12 = hmc_complex_zero;
                return out;
                break;
        }
#else

switch(random){
        case 1:
                out.e00 = src.e00;
                out.e01 = src.e01;
                out.e02 = hmc_complex_zero;
                out.e10 = src.e10;
                out.e11 = src.e11;
                out.e12 = hmc_complex_zero;
                out.e20 = hmc_complex_zero;
                out.e21 = hmc_complex_zero;
                out.e22 = hmc_complex_one;
                return out;
                break;
        case 2:
                out.e00 = hmc_complex_one;
                out.e01 = hmc_complex_zero;
                out.e02 = hmc_complex_zero;
                out.e10 = hmc_complex_zero;
                out.e11 = src.e00;
                out.e12 = src.e01;
                out.e20 = hmc_complex_zero;
                out.e21 = src.e10;
                out.e22 = src.e11;
                return out;
                break;
        case 3:
                out.e00 = src.e00;
                out.e01 = hmc_complex_zero;
                out.e02 = src.e01;
                out.e10= hmc_complex_zero;
                out.e11 = hmc_complex_one;
                out.e12 = hmc_complex_zero;
                out.e20 = src.e10;
                out.e21 = hmc_complex_zero;
                out.e22 = src.e11;
                return out;
                break;
        }
#endif
	return;
}

#else
  
#ifdef _RECONSTRUCT_TWELVE_
	if (random == 1) {
		out.e00 = src.e00;
		out.e01 = src.e01;
		out.e02 = hmc_complex_zero;
		out.e10 = src.e10;
		out.e11 = src.e11;
		out.e12 = hmc_complex_zero;
	} else if (random == 2) {
		out.e00 = hmc_complex_one;
		out.e01 = hmc_complex_zero;
		out.e02 = hmc_complex_zero;
		out.e10 = hmc_complex_zero;
		out.e11 = src.e00;
		out.e12 = src.e01;
	} else if (random == 3) {
		out.e00 = src.e00;
		out.e01 = hmc_complex_zero;
		out.e02 = src.e01;
		out.e10= hmc_complex_zero;
		out.e11 = hmc_complex_one;
		out.e12 = hmc_complex_zero;
	}
#else
	if (random == 1) {
		out.e00 = src.e00;
		out.e01 = src.e01;
		out.e02 = hmc_complex_zero;
		out.e10 = src.e10;
		out.e11 = src.e11;
		out.e12 = hmc_complex_zero;
		out.e20 = hmc_complex_zero;
		out.e21 = hmc_complex_zero;
		out.e22 = hmc_complex_one;
	} else if (random == 2) {
		out.e00 = hmc_complex_one;
		out.e01 = hmc_complex_zero;
		out.e02 = hmc_complex_zero;
		out.e10 = hmc_complex_zero;
		out.e11 = src.e00;
		out.e12 = src.e01;
		out.e20 = hmc_complex_zero;
		out.e21 = src.e10;
		out.e22 = src.e11;
	} else if (random == 3) {
		out.e00 = src.e00;
		out.e01 = hmc_complex_zero;
		out.e02 = src.e01;
		out.e10= hmc_complex_zero;
		out.e11 = hmc_complex_one;
		out.e12 = hmc_complex_zero;
		out.e20 = src.e10;
		out.e21 = hmc_complex_zero;
		out.e22 = src.e11;
	}
#endif

	return out;
}
#endif //_USE_GPU_

// void gaugefield_apply_bc(__private hmc_ocl_su3matrix * in, hmc_float theta)
// {
// 	hmc_float tmp1,tmp2;
// #ifdef _RECONSTRUCT_TWELVE_
// 	for(int a=0; a<NC; a++) {
// 		for(int b=0; b<NC-1; b++) {
// // 			for(int n=0; n<NC*(NC-1); n++) {
// 			tmp1 = ((in)[ocl_su3matrix_element(a,b)]).re;
// 			tmp2 = ((in)[ocl_su3matrix_element(a,b)]).im;
// 			((in)[ocl_su3matrix_element(a,b)]).re = cos(theta)*tmp1 - sin(theta)*tmp2;
// 			((in)[ocl_su3matrix_element(a,b)]).im = sin(theta)*tmp1 + cos(theta)*tmp2;
// 		}
// 	}
// #else
// 	for(int a=0; a<NC; a++) {
// 		for(int b=0; b<NC; b++) {
// 			tmp1 = ((in)[ocl_su3matrix_element(a,b)]).re;
// 			tmp2 = ((in)[ocl_su3matrix_element(a,b)]).im;
// 			((in)[ocl_su3matrix_element(a,b)]).re = cos(theta)*tmp1 - sin(theta)*tmp2;
// 			((in)[ocl_su3matrix_element(a,b)]).im = sin(theta)*tmp1 + cos(theta)*tmp2;
// 		}
// 	}
// #endif
// 	return;
// }

//CP: the following two chem-pot functions are one in the host-code
// void gaugefield_apply_chem_pot_real(__private hmc_ocl_su3matrix * u, __private hmc_ocl_su3matrix * udagger, hmc_float chem_pot_re)
// {
// 	hmc_float tmp1,tmp2;
// #ifdef _RECONSTRUCT_TWELVE_
// 	for(int a=0; a<NC; a++) {
// 		for(int b=0; b<NC-1; b++) {
// //   for(int n=0; n<NC*(NC-1); n++) {
// 			tmp1 = ((u)[ocl_su3matrix_element(a,b)]).re;
// 			tmp2 = ((u)[ocl_su3matrix_element(a,b)]).im;
// 			((u)[ocl_su3matrix_element(a,b)]).re = exp(chem_pot_re)*tmp1;
// 			((u)[ocl_su3matrix_element(a,b)]).im = exp(chem_pot_re)*tmp2;
// 			tmp1 = ((udagger)[ocl_su3matrix_element(a,b)]).re;
// 			tmp2 = ((udagger)[ocl_su3matrix_element(a,b)]).im;
// 			((udagger)[ocl_su3matrix_element(a,b)]).re = exp(-chem_pot_re)*tmp1;
// 			((udagger)[ocl_su3matrix_element(a,b)]).im = exp(-chem_pot_re)*tmp2;
// 		}
// 	}
// #else
// 	for(int a=0; a<NC; a++) {
// 		for(int b=0; b<NC; b++) {
// 			tmp1 = ((u)[ocl_su3matrix_element(a,b)]).re;
// 			tmp2 = ((u)[ocl_su3matrix_element(a,b)]).im;
// 			((u)[ocl_su3matrix_element(a,b)]).re = exp(chem_pot_re)*tmp1;
// 			((u)[ocl_su3matrix_element(a,b)]).im = exp(chem_pot_re)*tmp2;
// 			tmp1 = ((udagger)[ocl_su3matrix_element(a,b)]).re;
// 			tmp2 = ((udagger)[ocl_su3matrix_element(a,b)]).im;
// 			((udagger)[ocl_su3matrix_element(a,b)]).re = exp(-chem_pot_re)*tmp1;
// 			((udagger)[ocl_su3matrix_element(a,b)]).im = exp(-chem_pot_re)*tmp2;
// 		}
// 	}
// #endif
// 	return;
// }

// void gaugefield_apply_chem_pot_imag(__private hmc_ocl_su3matrix * u, __private hmc_ocl_su3matrix * udagger, hmc_float chem_pot_im)
// {
// 	hmc_float tmp1,tmp2;
// #ifdef _RECONSTRUCT_TWELVE_
// 	for(int a=0; a<NC; a++) {
// 		for(int b=0; b<NC-1; b++) {
// //   for(int n=0; n<NC*(NC-1); n++) {
// 			tmp1 = ((u)[ocl_su3matrix_element(a,b)]).re;
// 			tmp2 = ((u)[ocl_su3matrix_element(a,b)]).im;
// 			((u)[ocl_su3matrix_element(a,b)]).re = ( cos(chem_pot_im)*tmp1 - sin(chem_pot_im)*tmp2 );
// 			((u)[ocl_su3matrix_element(a,b)]).im = ( sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
// 			tmp1 = ((udagger)[ocl_su3matrix_element(a,b)]).re;
// 			tmp2 = ((udagger)[ocl_su3matrix_element(a,b)]).im;
// 			((udagger)[ocl_su3matrix_element(a,b)]).re = ( cos(chem_pot_im)*tmp1 + sin(chem_pot_im)*tmp2 );
// 			((udagger)[ocl_su3matrix_element(a,b)]).im = ( -sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
// 		}
// 	}
// #else
// 	for(int a=0; a<NC; a++) {
// 		for(int b=0; b<NC; b++) {
// 			tmp1 = ((u)[ocl_su3matrix_element(a,b)]).re;
// 			tmp2 = ((u)[ocl_su3matrix_element(a,b)]).im;
// 			((u)[ocl_su3matrix_element(a,b)]).re = ( cos(chem_pot_im)*tmp1 - sin(chem_pot_im)*tmp2 );
// 			((u)[ocl_su3matrix_element(a,b)]).im = ( sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
// 			tmp1 = ((udagger)[ocl_su3matrix_element(a,b)]).re;
// 			tmp2 = ((udagger)[ocl_su3matrix_element(a,b)]).im;
// 			((udagger)[ocl_su3matrix_element(a,b)]).re = ( cos(chem_pot_im)*tmp1 + sin(chem_pot_im)*tmp2 );
// 			((udagger)[ocl_su3matrix_element(a,b)]).im = ( -sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
// 		}
// 	}
// #endif
// 	return;
// }

//calculate polyakov-loop matrix at spatial site n in time-direction
Matrixsu3 local_polyakov(__global ocl_s_gaugefield * field, const int n)
{
	Matrixsu3 out;
	out = unit_matrixsu3();
	for(int t=0; t<NTIME; t++) {
		Matrixsu3 tmp;
		tmp = get_matrixsu3(field,n,t,0);
		out = multiply_matrixsu3 (out, tmp);
	}
	return out;
}

//CP: final version, using one third of the original scratch-registers...
//calculate plaquette-matrix at site n,t in direction mu and nu
Matrixsu3 local_plaquette(__global ocl_s_gaugefield * field, const int n, const int t, const int mu, const int nu )
{
        Matrixsu3 out;
        int4 pos;
        if(mu==0) {
                pos.x = (t+1)%NTIME;
                pos.y = n;
        } else {
               pos.x = t;
               pos.y = get_neighbor(n, mu);
        }
        if(nu==0) {
                pos.z = (t+1)%NTIME;
                pos.w = n;
        } else {
               pos.z = t;
               pos.w = get_neighbor(n, nu);
        }
        out = multiply_matrixsu3 (get_matrixsu3(field,n,t,mu), get_matrixsu3(field,pos.y,pos.x,nu)      );
        out = multiply_matrixsu3_dagger(out, get_matrixsu3(field,pos.w,pos.z,mu) );
        out = multiply_matrixsu3_dagger(out, get_matrixsu3(field,n,t,nu) );

        return out;
}


//todo
Matrix3x3 local_Q_plaquette(__global ocl_s_gaugefield * field, const int n, const int t, const int mu, const int nu ){
  Matrix3x3 out;
  Matrixsu3 tmp;
  int newpos;
  
  //1st plaq
  Matrixsu3 plaq1;
  //u_mu(x)
  plaq1 = get_matrixsu3(field,n,t,mu);
  //u_nu(x+mu)
  if(mu==0) {
    int newt = (t+1)%NTIME;
    tmp = get_matrixsu3(field,n,newt,nu);
  }
  else
    tmp = get_matrixsu3(field,get_neighbor(n,mu),t,nu);
  plaq1 = multiply_matrixsu3 (plaq1, tmp);
  //adjoint(u_mu(x+nu))
  if(nu==0) {
    int newt = (t+1)%NTIME;
    tmp = get_matrixsu3(field,n,newt,mu);
  }
  else
    tmp = get_matrixsu3(field,get_neighbor(n,nu),t,mu);
  tmp = adjoint_matrixsu3(tmp);
  plaq1 = multiply_matrixsu3 (plaq1, tmp);
  //adjoint(u_nu(x))
  tmp = get_matrixsu3(field,n,t,nu);
  tmp = adjoint_matrixsu3(tmp);
  plaq1 = multiply_matrixsu3 (plaq1, tmp);

  //2nd plaq
  Matrixsu3 plaq2;
  //U_nu(x)
  plaq2 = get_matrixsu3(field,n,t,nu);
  //adj (u_mu(x-mu+nu))
  newpos = get_lower_neighbor(n, mu);
  if (nu==0){
    int newt =  (t+1)%NTIME;
    tmp = get_matrixsu3(field,newpos,newt,mu);
  }
  else if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    tmp = get_matrixsu3(field,get_neighbor(n,nu),newt,mu);
  }
  else
    tmp = get_matrixsu3(field,get_neighbor(newpos, nu),t,mu);
  
  tmp = adjoint_matrixsu3(tmp);
  plaq2 = multiply_matrixsu3 (plaq2, tmp);
  //adj (u_nu(x-mu))
  if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    tmp = get_matrixsu3(field,n,newt, nu);
  }
  else
    tmp = get_matrixsu3(field, get_lower_neighbor(n, mu),t, nu);
  tmp = adjoint_matrixsu3(tmp);
  plaq2 = multiply_matrixsu3 (plaq2, tmp);
  //u_mu(x-mu)
  if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    tmp = get_matrixsu3(field,n,newt, mu);
  }
  else
    tmp = get_matrixsu3(field, get_lower_neighbor(n, mu),t, mu);
  plaq2 = multiply_matrixsu3 (plaq2, tmp);

  //3rd plaq
  Matrixsu3 plaq3;
  //adj (u_mu(x-mu))
  if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    tmp = get_matrixsu3(field,n,newt, mu);
  }
  else
    tmp = get_matrixsu3(field, get_lower_neighbor(n, mu),t, mu);
  tmp = adjoint_matrixsu3(tmp);
  plaq3 = copy_matrixsu3(tmp); 
  //adj (u_nu(x-mu-nu))
  newpos = get_lower_neighbor(n, mu);
  if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    tmp = get_matrixsu3(field, newpos,newt,nu);
  }
  else if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
  tmp = get_matrixsu3(field,get_lower_neighbor(n,nu),newt,nu);
  }
  else
    tmp = get_matrixsu3(field,get_lower_neighbor(newpos, nu),t,nu);
  tmp = adjoint_matrixsu3(tmp);
  plaq3 = multiply_matrixsu3 (plaq3, tmp);
  //u_mu(x-mu-nu)
  newpos = get_lower_neighbor(n, mu);
  if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    tmp = get_matrixsu3(field,newpos,newt,mu);
  }
  else if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    tmp = get_matrixsu3(field,get_lower_neighbor(n,nu),newt,mu);
  }
  else
    tmp = get_matrixsu3 (field,get_lower_neighbor(newpos, nu),t,mu);
  plaq3 = multiply_matrixsu3(plaq3, tmp);
  //u_nu(x-nu)
  if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    tmp = get_matrixsu3(field, n,newt, nu);
  }
  else
    //this causes for nu=1 a speicherzugriffsfehler
    tmp = get_matrixsu3(field,get_lower_neighbor(n, nu),t,nu);
  plaq3 = multiply_matrixsu3 (plaq3, tmp);

  //4th plaq
  Matrixsu3 plaq4;
  //adj(u_nu(x-nu))
  tmp = adjoint_matrixsu3(tmp);
  plaq4 = copy_matrixsu3(tmp); 
  //u_mu(x-nu)
   if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    tmp = get_matrixsu3(field, n,newt, mu);
  }
  else
    tmp = get_matrixsu3(field, get_lower_neighbor(n, nu),t,mu);
  plaq4 = multiply_matrixsu3 (plaq4, tmp);
  //u_nu(x+mu-nu)
  newpos = get_lower_neighbor(n, nu);
  if (mu==0){
    int newt =  (t+1)%NTIME;
    tmp = get_matrixsu3(field,newpos,newt,nu);
  }
  else if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    tmp = get_matrixsu3(field,get_neighbor(n,mu),newt,nu);
  }
  else
    tmp = get_matrixsu3(field,get_neighbor(newpos, mu),t,nu);
  plaq4 = multiply_matrixsu3(plaq4, tmp);
  //adj (u_mu(x))
  tmp = get_matrixsu3(field,n,t,mu);
  tmp = adjoint_matrixsu3(tmp);
  plaq4 = multiply_matrixsu3 (plaq4, tmp);

  //Sum up
  Matrix3x3 tmp3x3;
  out = matrix_su3to3x3 (plaq1);
  tmp3x3 = matrix_su3to3x3 (plaq2);
  out = add_matrix3x3 (out, tmp3x3);
  tmp3x3 = matrix_su3to3x3 (plaq3);
  out = add_matrix3x3 (out, tmp3x3);
  tmp3x3 = matrix_su3to3x3 (plaq4);
  out = add_matrix3x3 (out, tmp3x3);
  
  return out;
}

Matrix3x3 calc_staple(__global ocl_s_gaugefield* field, const int pos, const int t, const int mu_in)
{
	Matrixsu3 prod;
	Matrixsu3 prod2;
	Matrixsu3 tmp;
	Matrix3x3 staple;
	int nu, newpos, newt;

 	staple = zero_matrix3x3();
	
	//iterate through the three directions other than mu
 	for(int i = 1; i<NDIM; i++) {

		prod = zero_matrixsu3();
		
		nu = (mu_in + i)%NDIM;
		//first staple
		//u_nu(x+mu)
		if(mu_in==0) {
			newt = (t+1)%NTIME;
			tmp = get_matrixsu3(field,pos,newt,nu);

		} else {
			tmp = get_matrixsu3(field,get_neighbor(pos,mu_in),t,nu);
		}
		prod = copy_matrixsu3(tmp);
		
		//adjoint(u_mu(x+nu))
		if(nu==0) {
			newt = (t+1)%NTIME;
			tmp = get_matrixsu3(field,pos,newt,mu_in);
		} else {
			tmp = get_matrixsu3(field,get_neighbor(pos,nu),t,mu_in);
		}
		tmp = adjoint_matrixsu3(tmp);
		
		prod = multiply_matrixsu3(prod, tmp);
		
		//adjoint(u_nu(x))
		tmp = get_matrixsu3(field,pos,t,nu);
		tmp = adjoint_matrixsu3(tmp);
		prod = multiply_matrixsu3 (prod, tmp);
		
		//second staple
		//adjoint (u_nu(x+mu-nu))
		//newpos is "pos-nu" (spatial)
		newpos = get_lower_neighbor(pos, nu);
		if(mu_in==0) {
			newt = (t+1)%NTIME;
			tmp = get_matrixsu3(field,newpos,newt,nu);
		} else if (nu == 0) {
			newt = (t-1+NTIME)%NTIME;
			tmp = get_matrixsu3(field,get_neighbor(pos,mu_in),newt,nu);
		} else {
			tmp = get_matrixsu3(field,get_neighbor(newpos,mu_in),t,nu);
		}
		prod2 = adjoint_matrixsu3(tmp);
		//adjoint(u_mu(x-nu))
		if(mu_in==0) {
			tmp = get_matrixsu3(field,newpos,t,mu_in);
		} else if (nu == 0) {
			newt = (t-1+NTIME)%NTIME;
			tmp = get_matrixsu3(field,pos,newt,mu_in);
		} else {
			tmp = get_matrixsu3(field,newpos,t,mu_in);
		}
		tmp = adjoint_matrixsu3(tmp);
		prod2 = multiply_matrixsu3(prod2, tmp);
		//adjoint(u_nu(x-nu))
		if(mu_in==0) {
			tmp = get_matrixsu3(field,newpos,t,nu);
		} else if (nu == 0) {
			newt = (t-1+NTIME)%NTIME;
			tmp = get_matrixsu3(field,pos,newt,nu);
		} else {
			tmp = get_matrixsu3(field,newpos,t,nu);
		}
		prod2 = multiply_matrixsu3(prod2, tmp);

		
		Matrix3x3 dummy;
		dummy = matrix_su3to3x3 (prod);
		staple = add_matrix3x3 (staple, dummy );
		dummy = matrix_su3to3x3 (prod2);
		staple = add_matrix3x3 (staple, dummy );
	}
	
	return staple;
}
