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
	a[1] = U.e10;
	a[2] = U.e20;
	b[0] = U.e01;
	b[1] = U.e11;
	b[2] = U.e21;
	c[0] = U.e02;
	c[1] = U.e12;
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
		tmp = complexconj (a[i]);
		tmp = complexmult (b[i], tmp);
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
	tmp = complexconj(a[1]);
	tmp2 = complexconj (b[2]);
	c[0] = complexmult (tmp, tmp2);
	tmp = complexconj(a[2]);
	tmp2 = complexconj (b[1]);
	tmp = complexmult (tmp, tmp2);
	c[0] = complexsubtract (c[0],tmp);

	tmp = complexconj(a[2]);
	tmp2 = complexconj (b[0]);
	c[1] = complexmult (tmp, tmp2);
	tmp = complexconj(a[0]);
	tmp2 = complexconj (b[2]);
	tmp =  complexmult (tmp, tmp2);
	c[1] = complexsubtract (c[1],tmp);

	tmp = complexconj(a[0]);
	tmp2 = complexconj (b[1]);
	c[2] = complexmult (tmp, tmp2);
	tmp = complexconj(a[1]);
	tmp2 = complexconj (b[0]);
	tmp =  complexmult (tmp, tmp2);
	c[2] = complexsubtract (c[2],tmp);

	//Set new values to matrix
	out.e02 = c[0];
	out.e12 = c[1];
	out.e22 = c[2];
 
	//Set new values to matrix
	out.e01 = b[0];
	out.e11 = b[1];
	out.e21 = b[2];
	out.e00 = a[0];
	out.e10 = a[1];
	out.e20 = a[2];

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

//calculate plaquette-matrix at site n,t in direction mu and nu
Matrixsu3 local_plaquette(__global ocl_s_gaugefield * field, const int n, const int t, const int mu, const int nu )
{
	Matrixsu3 out;
	int4 pos;
	if(mu == 0) {
	        pos.x = get_neighbor_temporal(t);
		pos.y = n;
	} else {
		pos.x = t;
		pos.y = get_neighbor(n, mu);
	}
	if(nu == 0) {
                pos.z = get_neighbor_temporal(t);
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

Matrix3x3 local_Q_plaquette(__global ocl_s_gaugefield * field, const int n, const int t, const int mu, const int nu )
{
	//the Q-plaquette is a sum over four normal plaquettes
	Matrix3x3 qplaq = zero_matrix3x3();
	//first plaquette is at pos = (n,t)
	Matrixsu3 tmp = local_plaquette(field, n, t, mu, nu);
	qplaq = add_matrix3x3(qplaq, matrix_su3to3x3(tmp));
	//second plaquette is at pos - mu
	int4 coord;
	if(mu == 0){
		coord.x = n;
		coord.y = get_lower_neighbor_temporal(t);
	}
	else{
		coord.x = get_lower_neighbor(n, mu);
		coord.y = t;
	}
	tmp = local_plaquette(field, coord.x, coord.y, mu, nu);
	qplaq = add_matrix3x3(qplaq, matrix_su3to3x3(tmp));
	//third plaquette is at pos-mu-nu
	if(nu == 0) {
		coord.z = coord.x;
		coord.w = get_lower_neighbor_temporal(coord.y);
	}
	else{
		coord.z = get_lower_neighbor(coord.x, nu);
		coord.w = coord.y;
	}
	tmp = local_plaquette(field, coord.z, coord.w, mu, nu);
	qplaq = add_matrix3x3(qplaq, matrix_su3to3x3(tmp));
	//fourth plaquette is at pos-nu
	if(nu == 0) {
		coord.x = n;
		coord.y = get_lower_neighbor_temporal(t);
	}
	else{
		coord.x = get_lower_neighbor(n, nu);
		coord.y = t;
	}
	tmp = local_plaquette(field, coord.x, coord.y, mu, nu);
	qplaq = add_matrix3x3(qplaq, matrix_su3to3x3(tmp));
	return qplaq;
}

//this calculates the staple in nu direction given a direction mu of the link
//     under consideration:
//     s = U_nu(x + mu) * Udagger_mu(x + nu) * Udagger_nu(x) + Udagger_nu(x+mu - nu) * Udagger_mu(x-nu) * U_nu(x - nu)
Matrix3x3 local_staple(__global ocl_s_gaugefield * field, const int n, const int t, const int mu, const int nu )
{
	int4 pos;

	//first staple
	//calculate the coordinates for the matrices. this is the same as with	the plaquette
	if(mu == 0) {
  	        pos.x = get_neighbor_temporal(t);
		pos.y = n;
	} else {
	        pos.x = t;
		pos.y = get_neighbor(n, mu);
	}
	if(nu == 0) {
	        pos.z = get_neighbor_temporal(t);
		pos.w = n;
	} else {
		pos.z = t;
		pos.w = get_neighbor(n, nu);
	}
	Matrixsu3 tmp = multiply_matrixsu3_dagger(get_matrixsu3(field, pos.y, pos.x, nu), get_matrixsu3(field, pos.w, pos.z, mu) );
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, n, t, nu) );
	
	Matrix3x3 out = matrix_su3to3x3(tmp);

	//second staple
	if(mu == 0) {
	         pos.x = get_neighbor_temporal(t);
		 pos.y = get_lower_neighbor(n, nu);
		 pos.z = t;
		 pos.w = get_lower_neighbor(n, nu);
	} else if (nu == 0) {
	         pos.x = get_lower_neighbor_temporal(t);
		 pos.y = get_neighbor(n, mu);
		 pos.z = get_lower_neighbor_temporal(t);
		 pos.w = n;
	} else {
		 pos.x = t;
	       	 pos.y = get_neighbor(get_lower_neighbor(n, nu), mu);
		 pos.z = t;
		 pos.w = get_lower_neighbor(n, nu);
	}
	tmp = multiply_matrixsu3_dagger_dagger(get_matrixsu3(field, pos.y, pos.x, nu), get_matrixsu3(field, pos.w,pos.z, mu) );
	tmp = multiply_matrixsu3(tmp, get_matrixsu3(field, pos.w, pos.z, nu) );

	out = add_matrix3x3 (out, matrix_su3to3x3(tmp) );

	return out;
}

Matrix3x3 calc_staple(__global ocl_s_gaugefield* field, const int pos, const int t, const int mu_in)
{
	Matrix3x3 staple = zero_matrix3x3();
	//iterate through the three directions other than mu
	for(int i = 0; i < NDIM-1; i++) {
		int nu = (mu_in + i + 1) % NDIM;
		staple = add_matrix3x3(staple,  local_staple(field, pos, t, mu_in, nu ));
	}
	return staple;
}

//this is the staple only in the spatial directions only
Matrix3x3 calc_staple_sigma (__global ocl_s_gaugefield* field, const int pos, const int t, const int mu_in)
{
	Matrix3x3 staple = zero_matrix3x3();
	//iterate through the three directions other than mu
	for(int i = 0; i < NDIM-1; i++) {
		int nu = (mu_in + i + 1) % NDIM;
		if(nu!=0)
			staple = add_matrix3x3(staple,  local_staple(field, pos, t, mu_in, nu ));
	}
	return staple;
}

//this is the staple only in temporal direction only
Matrix3x3 calc_staple_tau (__global ocl_s_gaugefield* field, const int pos, const int t, const int mu_in)
{
	int nu = 0;
	return local_staple(field, pos, t, mu_in, nu );
}

// TODO document
Matrixsu3 getSU3SOA(__global const Matrixsu3StorageType * const restrict in, const uint idx)
{
#ifdef _USE_SOA_
	return (Matrixsu3) {
		in[0 * GAUGEFIELD_STRIDE + idx],
		in[1 * GAUGEFIELD_STRIDE + idx],
		in[2 * GAUGEFIELD_STRIDE + idx],
		in[3 * GAUGEFIELD_STRIDE + idx],
		in[4 * GAUGEFIELD_STRIDE + idx],
		in[5 * GAUGEFIELD_STRIDE + idx],
		in[6 * GAUGEFIELD_STRIDE + idx],
		in[7 * GAUGEFIELD_STRIDE + idx],
		in[8 * GAUGEFIELD_STRIDE + idx]
	};
#else
	return in[idx];
#endif
}

// TODO document
void putSU3SOA(__global Matrixsu3StorageType * const restrict out, const uint idx, const Matrixsu3 val)
{
#ifdef _USE_SOA_
	out[0 * GAUGEFIELD_STRIDE + idx] = val.e00;
	out[1 * GAUGEFIELD_STRIDE + idx] = val.e01;
	out[2 * GAUGEFIELD_STRIDE + idx] = val.e02;
	out[3 * GAUGEFIELD_STRIDE + idx] = val.e10;
	out[4 * GAUGEFIELD_STRIDE + idx] = val.e11;
	out[5 * GAUGEFIELD_STRIDE + idx] = val.e12;
	out[6 * GAUGEFIELD_STRIDE + idx] = val.e20;
	out[7 * GAUGEFIELD_STRIDE + idx] = val.e21;
	out[8 * GAUGEFIELD_STRIDE + idx] = val.e22;
#else
	out[idx] = val;
#endif
}

