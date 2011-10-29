/** @file
 * Device code for operations on the fermion matrix
 */

//operations_gaugefield.cl

Matrixsu3 project_su3(const Matrixsu3 U)
{

	Matrixsu3 out;

	//Extract initial vectors
	hmc_complex a[NC];
	hmc_complex b[NC];
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

	//New SU3-Matrix
	//first vector
	//norm
	hmc_float norm = 0.;
	for (int i = 0; i < NC; i++) {
		hmc_complex tmp = complexconj(a[i]);
		tmp = complexmult (a[i], tmp);
		norm += tmp.re;
	}
	norm = 1. / sqrt(norm);
	//rescale
	for (int i = 0; i < NC; i++) {
		a[i].re *= norm;
		a[i].im *= norm;
	}

	//second vector
	//orthogonal vector
	hmc_complex factor;
	factor.re = 0.0;
	factor.im = 0.0;
	for (int i = 0; i < NC; i++) {
		hmc_complex tmp;
		tmp = complexconj (b[i]);
		tmp = complexmult (a[i], tmp);
		factor = complexadd (factor, tmp);
	}
	for (int i = 0; i < NC; i++) {
		hmc_complex tmp;
		tmp = complexmult(factor, a[i]);
		b[i] = complexsubtract(b[i], tmp);
	}

//norm
	norm = 0.;
	for (int i = 0; i < NC; i++) {
		hmc_complex tmp;
		tmp = complexconj(b[i]);
		tmp = complexmult (b[i], tmp);
		norm +=  tmp.re;
	}
	norm = 1. / sqrt(norm);
	//rescale
	for  (int i = 0; i < NC; i++) {
		b[i].re *= norm;
		b[i].im *= norm;
	}

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
	} else if (rand == 2) {
		out.e00 = src.e11;
		out.e01 = src.e12;
		out.e10 = src.e21;
		out.e11 = src.e22;
	} else if (rand == 3) {
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

	switch(random) {
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
		case 3:
			out.e00 = src.e00;
			out.e01 = hmc_complex_zero;
			out.e02 = src.e01;
			out.e10 = hmc_complex_zero;
			out.e11 = hmc_complex_one;
			out.e12 = hmc_complex_zero;
			out.e20 = src.e10;
			out.e21 = hmc_complex_zero;
			out.e22 = src.e11;
			return out;
	}
}

//calculate polyakov-loop matrix at spatial site n in time-direction
Matrixsu3 local_polyakov(__global ocl_s_gaugefield * field, const int n)
{
	Matrixsu3 out;
	out = unit_matrixsu3();
	for(int t = 0; t < NTIME; t++) {
		Matrixsu3 tmp;
		tmp = get_matrixsu3(field, n, t, 0);
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
	if(mu == 0) {
		pos.x = (t + 1) % NTIME;
		pos.y = n;
	} else {
		pos.x = t;
		pos.y = get_neighbor(n, mu);
	}
	if(nu == 0) {
		pos.z = (t + 1) % NTIME;
		pos.w = n;
	} else {
		pos.z = t;
		pos.w = get_neighbor(n, nu);
	}
	out = multiply_matrixsu3 (get_matrixsu3(field, n, t, mu), get_matrixsu3(field, pos.y, pos.x, nu)      );
	out = multiply_matrixsu3_dagger(out, get_matrixsu3(field, pos.w, pos.z, mu) );
	out = multiply_matrixsu3_dagger(out, get_matrixsu3(field, n, t, nu) );

	return out;
}


//todo
Matrix3x3 local_Q_plaquette(__global ocl_s_gaugefield * field, const int n, const int t, const int mu, const int nu )
{
	Matrix3x3 out;
	Matrixsu3 tmp;
	int newpos;

	//1st plaq
	Matrixsu3 plaq1;
	//u_mu(x)
	plaq1 = get_matrixsu3(field, n, t, mu);
	//u_nu(x+mu)
	if(mu == 0) {
		int newt = (t + 1) % NTIME;
		tmp = get_matrixsu3(field, n, newt, nu);
	} else
		tmp = get_matrixsu3(field, get_neighbor(n, mu), t, nu);
	plaq1 = multiply_matrixsu3 (plaq1, tmp);
	//adjoint(u_mu(x+nu))
	if(nu == 0) {
		int newt = (t + 1) % NTIME;
		tmp = get_matrixsu3(field, n, newt, mu);
	} else
		tmp = get_matrixsu3(field, get_neighbor(n, nu), t, mu);
	tmp = adjoint_matrixsu3(tmp);
	plaq1 = multiply_matrixsu3 (plaq1, tmp);
	//adjoint(u_nu(x))
	tmp = get_matrixsu3(field, n, t, nu);
	tmp = adjoint_matrixsu3(tmp);
	plaq1 = multiply_matrixsu3 (plaq1, tmp);

	//2nd plaq
	Matrixsu3 plaq2;
	//U_nu(x)
	plaq2 = get_matrixsu3(field, n, t, nu);
	//adj (u_mu(x-mu+nu))
	newpos = get_lower_neighbor(n, mu);
	if (nu == 0) {
		int newt =  (t + 1) % NTIME;
		tmp = get_matrixsu3(field, newpos, newt, mu);
	} else if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, get_neighbor(n, nu), newt, mu);
	} else
		tmp = get_matrixsu3(field, get_neighbor(newpos, nu), t, mu);

	tmp = adjoint_matrixsu3(tmp);
	plaq2 = multiply_matrixsu3 (plaq2, tmp);
	//adj (u_nu(x-mu))
	if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, n, newt, nu);
	} else
		tmp = get_matrixsu3(field, get_lower_neighbor(n, mu), t, nu);
	tmp = adjoint_matrixsu3(tmp);
	plaq2 = multiply_matrixsu3 (plaq2, tmp);
	//u_mu(x-mu)
	if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, n, newt, mu);
	} else
		tmp = get_matrixsu3(field, get_lower_neighbor(n, mu), t, mu);
	plaq2 = multiply_matrixsu3 (plaq2, tmp);

	//3rd plaq
	Matrixsu3 plaq3;
	//adj (u_mu(x-mu))
	if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, n, newt, mu);
	} else
		tmp = get_matrixsu3(field, get_lower_neighbor(n, mu), t, mu);
	tmp = adjoint_matrixsu3(tmp);
	plaq3 = copy_matrixsu3(tmp);
	//adj (u_nu(x-mu-nu))
	newpos = get_lower_neighbor(n, mu);
	if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, newpos, newt, nu);
	} else if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, get_lower_neighbor(n, nu), newt, nu);
	} else
		tmp = get_matrixsu3(field, get_lower_neighbor(newpos, nu), t, nu);
	tmp = adjoint_matrixsu3(tmp);
	plaq3 = multiply_matrixsu3 (plaq3, tmp);
	//u_mu(x-mu-nu)
	newpos = get_lower_neighbor(n, mu);
	if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, newpos, newt, mu);
	} else if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, get_lower_neighbor(n, nu), newt, mu);
	} else
		tmp = get_matrixsu3 (field, get_lower_neighbor(newpos, nu), t, mu);
	plaq3 = multiply_matrixsu3(plaq3, tmp);
	//u_nu(x-nu)
	if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, n, newt, nu);
	} else
		//this causes for nu=1 a speicherzugriffsfehler
		tmp = get_matrixsu3(field, get_lower_neighbor(n, nu), t, nu);
	plaq3 = multiply_matrixsu3 (plaq3, tmp);

	//4th plaq
	Matrixsu3 plaq4;
	//adj(u_nu(x-nu))
	tmp = adjoint_matrixsu3(tmp);
	plaq4 = copy_matrixsu3(tmp);
	//u_mu(x-nu)
	if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, n, newt, mu);
	} else
		tmp = get_matrixsu3(field, get_lower_neighbor(n, nu), t, mu);
	plaq4 = multiply_matrixsu3 (plaq4, tmp);
	//u_nu(x+mu-nu)
	newpos = get_lower_neighbor(n, nu);
	if (mu == 0) {
		int newt =  (t + 1) % NTIME;
		tmp = get_matrixsu3(field, newpos, newt, nu);
	} else if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, get_neighbor(n, mu), newt, nu);
	} else
		tmp = get_matrixsu3(field, get_neighbor(newpos, mu), t, nu);
	plaq4 = multiply_matrixsu3(plaq4, tmp);
	//adj (u_mu(x))
	tmp = get_matrixsu3(field, n, t, mu);
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

//this calculates the staple in nu direction given a direction mu of the link
//     under consideration:
//     s = U_nu(x + mu) * Udagger_mu(x + nu) * Udagger_nu(x) + Udagger_nu(x+mu - nu) * Udagger_mu(x-nu) * U_nu(x - nu)
Matrix3x3 local_staple(__global ocl_s_gaugefield * field, const int n, const int t, const int mu, const int nu )
{
	Matrix3x3 out;
	Matrixsu3 tmp;
	int4 pos;

	//first staple
	//calculate the coordinates for the matrices. this is the same as with	the plaquette
	if(mu == 0) {
		pos.x = (t + 1) % NTIME;
		pos.y = n;
	} else {
	        pos.x = t;
		pos.y = get_neighbor(n, mu);
	}
	if(nu == 0) {
		pos.z = (t + 1) % NTIME;
		pos.w = n;
	} else {
		pos.z = t;
		pos.w = get_neighbor(n, nu);
	}
	tmp = multiply_matrixsu3_dagger(get_matrixsu3(field, pos.y, pos.x, nu), get_matrixsu3(field, pos.w, pos.z, mu) );
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, n, t, nu) );
	
	out = matrix_su3to3x3(tmp);

	//second staple
	if(mu == 0) {
		 pos.x = (t + 1) % NTIME;
		 pos.y = get_lower_neighbor(n, nu);
		 pos.z = t;
		 pos.w = get_lower_neighbor(n, nu);
	} else if (nu == 0) {
	  	 pos.x = (t - 1 + NTIME) % NTIME;
		 pos.y = get_neighbor(n, mu);
		 pos.z = (t-1+NTIME)% NTIME;
		 pos.w = n;
	} else {
		 pos.x = t;
	       	 pos.y = get_neighbor(get_lower_neighbor(n, nu), mu);
		 pos.z = t;
		 pos.w = get_lower_neighbor(n, nu);
	}

	//@TODO for this one can write a new function that calculates Udagger*Vdagger
	tmp = adjoint_matrixsu3(get_matrixsu3(field, pos.y, pos.x, nu));
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, pos.w, pos.z, mu) );
	tmp = multiply_matrixsu3(tmp, get_matrixsu3(field, pos.w, pos.z, nu) );

	out = add_matrix3x3 (out, matrix_su3to3x3(tmp) );

	return out;
}


Matrix3x3 calc_staple(__global ocl_s_gaugefield* field, const int pos, const int t, const int mu_in)
{
	Matrix3x3 staple;
	int nu;

	staple = zero_matrix3x3();
	//iterate through the three directions other than mu
	for(int i = 1; i < NDIM; i++) {
		nu = (mu_in + i) % NDIM;
		staple = add_matrix3x3(staple,  local_staple(field, pos, t, mu_in, nu ));
	}

	return staple;
}

/*
Matrix3x3 calc_staple_sigma (__global ocl_s_gaugefield* field, const int pos, const int t, const int mu_in)
{
	Matrixsu3 prod;
	Matrixsu3 prod2;
	Matrixsu3 tmp;
	Matrix3x3 staple;
	int nu;

	staple = zero_matrix3x3();

	//iterate through the three directions other than mu, dont include the temporal direction since it is the spatial staple
	for(int i = 1; i < NDIM; i++) {

	  prod = zero_matrixsu3();

	  nu = (mu_in + i) % NDIM;

	  if(nu != 0)
	  {	    

		//first staple
		//u_nu(x+mu)
		tmp = get_matrixsu3(field, get_neighbor(pos, mu_in), t, nu);
		prod = copy_matrixsu3(tmp);

		//adjoint(u_mu(x+nu))
		tmp = get_matrixsu3(field, get_neighbor(pos, nu), t, mu_in);
		tmp = adjoint_matrixsu3(tmp);

		prod = multiply_matrixsu3(prod, tmp);

		//adjoint(u_nu(x))
		tmp = get_matrixsu3(field, pos, t, nu);
		tmp = adjoint_matrixsu3(tmp);
		prod = multiply_matrixsu3 (prod, tmp);

		//second staple
		//adjoint (u_nu(x+mu-nu))
		//newpos is "pos-nu" (spatial)
		int newpos = get_lower_neighbor(pos, nu);
		tmp = get_matrixsu3(field, get_neighbor(newpos, mu_in), t, nu);
		
		prod2 = adjoint_matrixsu3(tmp);
		//adjoint(u_mu(x-nu))
		tmp = get_matrixsu3(field, newpos, t, mu_in);
		
		tmp = adjoint_matrixsu3(tmp);
		prod2 = multiply_matrixsu3(prod2, tmp);
		//adjoint(u_nu(x-nu))
		tmp = get_matrixsu3(field, newpos, t, nu);
		prod2 = multiply_matrixsu3(prod2, tmp);

		Matrix3x3 dummy;
		dummy = matrix_su3to3x3 (prod);
		staple = add_matrix3x3 (staple, dummy );
		dummy = matrix_su3to3x3 (prod2);
		staple = add_matrix3x3 (staple, dummy );
	  }
	}
	return staple;
}
*/
/*

Matrix3x3 calc_staple_tau (__global ocl_s_gaugefield* field, const int pos, const int t, const int mu_in)
{
	Matrixsu3 prod;
	Matrixsu3 prod2;
	Matrixsu3 tmp;
	Matrix3x3 staple;
	int nu;

	staple = zero_matrix3x3();

	prod = zero_matrixsu3();

	nu = 0;

		//first staple
		//u_nu(x+mu)
		tmp = get_matrixsu3(field, get_neighbor(pos, mu_in), t, nu);
		prod = copy_matrixsu3(tmp);

		//adjoint(u_mu(x+nu))
		int newt = (t + 1) % NTIME;
		tmp = get_matrixsu3(field, pos, newt, mu_in);
		tmp = adjoint_matrixsu3(tmp);

		prod = multiply_matrixsu3(prod, tmp);

		//adjoint(u_nu(x))
		tmp = get_matrixsu3(field, pos, t, nu);
		tmp = adjoint_matrixsu3(tmp);
		prod = multiply_matrixsu3 (prod, tmp);

		//second staple
		//adjoint (u_nu(x+mu-nu))
		//newpos is "pos-nu" (spatial)
		int newpos = get_lower_neighbor(pos, nu);
		newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, get_neighbor(pos, mu_in), newt, nu);
		
		prod2 = adjoint_matrixsu3(tmp);
		//adjoint(u_mu(x-nu))
		newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, pos, newt, mu_in);
		
		tmp = adjoint_matrixsu3(tmp);
		prod2 = multiply_matrixsu3(prod2, tmp);
		//adjoint(u_nu(x-nu))
		newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, pos, newt, nu);
		
		prod2 = multiply_matrixsu3(prod2, tmp);

		Matrix3x3 dummy;
		dummy = matrix_su3to3x3 (prod);
		staple = add_matrix3x3 (staple, dummy );
		dummy = matrix_su3to3x3 (prod2);
		staple = add_matrix3x3 (staple, dummy );
	return staple;
}
*/


