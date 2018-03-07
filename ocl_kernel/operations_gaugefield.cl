/*
 * Copyright (c) 2011,2012 Christopher Pinke
 * Copyright (c) 2011-2013 Matthias Bach
 * Copyright (c) 2011 Christian Schäfer
 * Copyright (c) 2018 Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

/** @file
 * Device code for operations on the fermion matrix
 */

//operations_gaugefield.cl

// TODO document
inline Matrixsu3 getSU3(__global const Matrixsu3StorageType * const restrict in, const uint idx)
{
#ifdef _USE_SOA_
#ifdef _USE_REC12_
	//CP: this should be the most straightforward implementation
	return (Matrixsu3) {
		in[0 * GAUGEFIELD_STRIDE + idx],
		   in[1 * GAUGEFIELD_STRIDE + idx],
		   in[2 * GAUGEFIELD_STRIDE + idx],
		   in[3 * GAUGEFIELD_STRIDE + idx],
		   in[4 * GAUGEFIELD_STRIDE + idx],
		   in[5 * GAUGEFIELD_STRIDE + idx],
		   complexsubtract( complexmult(in[1 * GAUGEFIELD_STRIDE + idx], in[5 * GAUGEFIELD_STRIDE + idx]),  complexmult(in[2 * GAUGEFIELD_STRIDE + idx], in[4 * GAUGEFIELD_STRIDE + idx]) ),
		   complexsubtract( complexmult(in[2 * GAUGEFIELD_STRIDE + idx], in[3 * GAUGEFIELD_STRIDE + idx]),  complexmult(in[0 * GAUGEFIELD_STRIDE + idx], in[5 * GAUGEFIELD_STRIDE + idx]) ),
		   complexsubtract( complexmult(in[0 * GAUGEFIELD_STRIDE + idx], in[4 * GAUGEFIELD_STRIDE + idx]),  complexmult(in[1 * GAUGEFIELD_STRIDE + idx], in[3 * GAUGEFIELD_STRIDE + idx]) )
	};
#else // _USE_REC12_
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
#endif // _USE_REC12_
#else  // _USE_SOA_
	//printf("%i\n", idx);
	return in[idx];
#endif
}



// TODO document
inline void putSU3(__global Matrixsu3StorageType * const restrict out, const uint idx, const Matrixsu3 val)
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

inline Matrixsu3 get_matrixsu3(__global const Matrixsu3StorageType * const restrict field, const int spacepos, const int timepos, const int mu)
{
	uint idx = get_link_pos(mu, spacepos, timepos);
	return getSU3(field, idx);
}

inline void put_matrixsu3(__global Matrixsu3StorageType  * const restrict field, const Matrixsu3 in, const int spacepos, const int timepos, const int mu)
{
	uint idx = get_link_pos(mu, spacepos, timepos);
	putSU3(field, idx, in);
}

inline Matrix3x3 get3x3(__global const Matrix3x3StorageType * const restrict in, const uint idx)
{
#ifdef _USE_SOA_
	return (Matrix3x3) {
		in[0 * GAUGEFIELD_3X3_STRIDE + idx],
		   in[1 * GAUGEFIELD_3X3_STRIDE + idx],
		   in[2 * GAUGEFIELD_3X3_STRIDE + idx],
		   in[3 * GAUGEFIELD_3X3_STRIDE + idx],
		   in[4 * GAUGEFIELD_3X3_STRIDE + idx],
		   in[5 * GAUGEFIELD_3X3_STRIDE + idx],
		   in[6 * GAUGEFIELD_3X3_STRIDE + idx],
		   in[7 * GAUGEFIELD_3X3_STRIDE + idx],
		   in[8 * GAUGEFIELD_3X3_STRIDE + idx]
	};
#else  // _USE_SOA_
	//printf("%i\n", idx);
	return in[idx];
#endif
}

inline void put3x3(__global Matrix3x3StorageType * const restrict out, const uint idx, const Matrix3x3 val)
{
#ifdef _USE_SOA_
	out[0 * GAUGEFIELD_3X3_STRIDE + idx] = val.e00;
	out[1 * GAUGEFIELD_3X3_STRIDE + idx] = val.e01;
	out[2 * GAUGEFIELD_3X3_STRIDE + idx] = val.e02;
	out[3 * GAUGEFIELD_3X3_STRIDE + idx] = val.e10;
	out[4 * GAUGEFIELD_3X3_STRIDE + idx] = val.e11;
	out[5 * GAUGEFIELD_3X3_STRIDE + idx] = val.e12;
	out[6 * GAUGEFIELD_3X3_STRIDE + idx] = val.e20;
	out[7 * GAUGEFIELD_3X3_STRIDE + idx] = val.e21;
	out[8 * GAUGEFIELD_3X3_STRIDE + idx] = val.e22;
#else
	out[idx] = val;
#endif
}

inline Matrix3x3 get_matrix3x3(__global const Matrix3x3StorageType * const restrict field, const int spacepos, const int timepos, const int mu)
{
	uint idx = get_link_pos(mu, spacepos, timepos);
	return get3x3(field, idx);
}

inline void put_matrix3x3(__global Matrix3x3StorageType  * const restrict field, const Matrix3x3 in, const int spacepos, const int timepos, const int mu)
{
	uint idx = get_link_pos(mu, spacepos, timepos);
	put3x3(field, idx, in);
}
inline Matrixsu3 project_su3(const Matrixsu3 U)
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
	c[0] = complexsubtract (c[0], tmp);

	tmp = complexconj(a[2]);
	tmp2 = complexconj (b[0]);
	c[1] = complexmult (tmp, tmp2);
	tmp = complexconj(a[0]);
	tmp2 = complexconj (b[2]);
	tmp =  complexmult (tmp, tmp2);
	c[1] = complexsubtract (c[1], tmp);

	tmp = complexconj(a[0]);
	tmp2 = complexconj (b[1]);
	c[2] = complexmult (tmp, tmp2);
	tmp = complexconj(a[1]);
	tmp2 = complexconj (b[0]);
	tmp =  complexmult (tmp, tmp2);
	c[2] = complexsubtract (c[2], tmp);

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

inline Matrixsu2 reduction (const Matrix3x3 src, const int rand)
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

inline Matrixsu3 extend (const int random, Matrixsu2 src)
{
// Yes, it is poor madness having two different variants for CPU and GPU at this place.
// However, on Catalyst 12.3 the default value of the switch statement messes up the GPU
// results and the if-else change trigges a SEGFAULT in a completely different place of
// the heatbath kernel when running on the CPU.
#ifdef _USEGPU_
	if(random == 1) {
		Matrixsu3 out;
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
	} else if(random == 2) {
		Matrixsu3 out;
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
	} else if(random == 3) {
		Matrixsu3 out;
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
	} else {
		return (Matrixsu3) { {nan((uint) 0), nan((uint) 0)}, {nan((uint) 0), nan((uint) 0)}, {nan((uint) 0), nan((uint) 0)},
			{nan((uint) 0), nan((uint) 0)}, {nan((uint) 0), nan((uint) 0)}, {nan((uint) 0), nan((uint) 0)},
			{nan((uint) 0), nan((uint) 0)}, {nan((uint) 0), nan((uint) 0)}, {nan((uint) 0), nan((uint) 0)}
		};
	}
#else // _USEGPU_
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

	return (Matrixsu3) { {nan((uint) 0), nan((uint) 0)}, {nan((uint) 0), nan((uint) 0)}, {nan((uint) 0), nan((uint) 0)},
		{nan((uint) 0), nan((uint) 0)}, {nan((uint) 0), nan((uint) 0)}, {nan((uint) 0), nan((uint) 0)},
		{nan((uint) 0), nan((uint) 0)}, {nan((uint) 0), nan((uint) 0)}, {nan((uint) 0), nan((uint) 0)}
	};
#endif // _USEGPU_
}

//calculate polyakov-loop matrix at spatial site n in time-direction
inline Matrixsu3 local_polyakov(__global const Matrixsu3StorageType * const restrict field, const int n)
{
	Matrixsu3 out;
	out = unit_matrixsu3();
	for(int t = 0; t < NTIME_LOCAL; t++) {
		Matrixsu3 tmp;
		tmp = get_matrixsu3(field, n, t, 0);
		out = multiply_matrixsu3 (out, tmp);
	}
	return out;
}

//calculate plaquette-matrix at site n,t in direction mu and nu
inline Matrixsu3 local_plaquette(__global const Matrixsu3StorageType * const restrict field, const int n, const int t, const int mu, const int nu )
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

//calculate rectangle-matrix at site i = (n,t) in direction mu and nu
//	The rectangle is then: U_mu(i) * U_nu(i+mu) * U_nu(i + mu + nu) * U_mu(i + nu + nu)dagger * U_nu(i + nu)dagger * U_nu(i)dagger
inline Matrixsu3 local_rectangles(__global const Matrixsu3StorageType * const restrict field, const int n, const int t, const int mu, const int nu )
{
	Matrixsu3 out;
	int4 pos;
	int4 pos2;
	//calc positions of neighbors
	//(x,y) = i + mu
	if(mu == 0) {
		pos.x = get_neighbor_temporal(t);
		pos.y = n;
	} else {
		pos.x = t;
		pos.y = get_neighbor(n, mu);
	}
	//(w,z) = i + nu
	if(nu == 0) {
		pos.z = get_neighbor_temporal(t);
		pos.w = n;
	} else {
		pos.z = t;
		pos.w = get_neighbor(n, nu);
	}
	//(x2, y2) = i + mu + nu = (x,y) + nu
	if(nu == 0) {
		pos2.x = get_neighbor_temporal(pos.x);
		pos2.y = pos.y;
	} else {
		pos2.x = pos.x;
		pos2.y = get_neighbor(pos.y, nu);
	}
	//(w2, z2) = i + nu + nu = (w,z) + nu
	if(nu == 0) {
		pos2.z = get_neighbor_temporal(pos.z);
		pos2.w = pos.w;
	} else {
		pos2.z = pos.z;
		pos2.w = get_neighbor(pos.w, nu);
	}
	out = multiply_matrixsu3 (get_matrixsu3(field, n, t, mu), get_matrixsu3(field, pos.y, pos.x, nu)      );
	out = multiply_matrixsu3 (out, get_matrixsu3(field, pos2.y, pos2.x, nu)      );
	out = multiply_matrixsu3_dagger(out, get_matrixsu3(field, pos2.w, pos2.z, mu) );
	out = multiply_matrixsu3_dagger(out, get_matrixsu3(field, pos.w, pos.z, nu) );
	out = multiply_matrixsu3_dagger(out, get_matrixsu3(field, n, t, nu) );

	return out;
}


inline Matrix3x3 local_Q_plaquette(__global const Matrixsu3StorageType * const restrict field, const int n, const int t, const int mu, const int nu )
{
	//the Q-plaquette is a sum over four normal plaquettes
	Matrix3x3 qplaq = zero_matrix3x3();
	//first plaquette is at pos = (n,t)
	Matrixsu3 tmp = local_plaquette(field, n, t, mu, nu);
	qplaq = add_matrix3x3(qplaq, matrix_su3to3x3(tmp));
	//second plaquette is at pos - mu
	int4 coord;
	if(mu == 0) {
		coord.x = n;
		coord.y = get_lower_neighbor_temporal(t);
	} else {
		coord.x = get_lower_neighbor(n, mu);
		coord.y = t;
	}
	tmp = local_plaquette(field, coord.x, coord.y, mu, nu);
	qplaq = add_matrix3x3(qplaq, matrix_su3to3x3(tmp));
	//third plaquette is at pos-mu-nu
	if(nu == 0) {
		coord.z = coord.x;
		coord.w = get_lower_neighbor_temporal(coord.y);
	} else {
		coord.z = get_lower_neighbor(coord.x, nu);
		coord.w = coord.y;
	}
	tmp = local_plaquette(field, coord.z, coord.w, mu, nu);
	qplaq = add_matrix3x3(qplaq, matrix_su3to3x3(tmp));
	//fourth plaquette is at pos-nu
	if(nu == 0) {
		coord.x = n;
		coord.y = get_lower_neighbor_temporal(t);
	} else {
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
inline void local_staple(Matrix3x3 * const restrict aggregate, __global const Matrixsu3StorageType * const restrict field, const int n, const int t, const int mu, const int nu )
{
	int4 pos;

	//first staple
	//calculate the coordinates for the matrices. this is the same as with  the plaquette
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
	Matrixsu3 tmp = get_matrixsu3(field, pos.y, pos.x, nu);
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, pos.w, pos.z, mu));
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, n, t, nu) );

	aggregate_matrix3x3(aggregate, matrix_su3to3x3(tmp));

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
	tmp = get_matrixsu3(field, pos.y, pos.x, nu);
	tmp = multiply_matrixsu3_dagger_dagger(tmp, get_matrixsu3(field, pos.w, pos.z, mu));
	tmp = multiply_matrixsu3(tmp, get_matrixsu3(field, pos.w, pos.z, nu) );

	aggregate_matrix3x3(aggregate, matrix_su3to3x3(tmp));
}

/**     This calculates the rectangles staple in nu direction given a direction mu of the link
 *  under consideration. There are six ingredients to the staple
 *  (I write them in a way that they form a path, altough this is not really necessary since they are under a trace):
 *  1.  U_nu(x + mu) * U_nu(x + mu + nu) * Udagger_mu(x + nu + nu) * Udagger_nu(x + nu) * Udagger_nu(x)
 *  2.  U_nu(x - nu) * U_mu(x + mu) * Udagger_nu(x - nu + mu + mu) Udagger_mu(x - nu + mu) Udagger_mu(x - nu)
 *  ->  U_mu(x + mu) * Udagger_nu(x - nu + mu + mu) Udagger_mu(x - nu + mu) Udagger_mu(x - nu) * U_nu(x - nu)
 *  3.  U_nu(x - mu - nu) * U_mu(x - mu) * Udagger_nu(x - nu + mu) * Udagger_mu(x - nu) * Udagger_mu(x - mu - nu)
 *  ->  Udagger_nu(x - nu + mu) * Udagger_mu(x - nu) * Udagger_mu(x - mu - nu) * U_nu(x - mu - nu) * U_mu(x - mu)
 *      NOTE: The last three have to be "daggered"
 *  4.  U_mu(x - nu - nu) * U_nu(x - nu - nu + mu) * U_nu(x - nu + mu) * Udagger_nu(x - nu) * Udagger_nu(x - nu - nu)
 *  ->  Udagger_nu(x - nu) * Udagger_nu(x - nu - nu) * U_mu(x - nu - nu) * U_nu(x - nu - nu + mu) * U_nu(x - nu + mu)
 *  ^+  Udagger_nu(x - nu + mu) * Udagger_nu(x - nu - nu + mu) * Udagger_mu(x - nu - nu) * Ud_nu(x - nu - nu) * U_nu(x - nu)
 *  5.  U_nu(x - mu) * U_mu(x - mu + nu) * U_mu(x + nu) * Udagger_nu(x + mu) * Udagger_mu(x - mu)
 *  ->  Udagger_mu(x - mu) * U_nu(x - mu) * U_mu(x - mu + nu) * U_mu(x + nu) * Udagger_nu(x + mu)
 *  ^+  U_nu(x + mu) * Udagger_mu(x + nu) * Udagger_mu(x - mu + nu) * Udagger_nu(x - mu) * U_mu(x - mu)
 *  6.  U_nu(x) * U_mu(x + nu) * U_mu(x + nu + mu) * Udagger_nu(x + mu + mu) * Udagger_mu(x + mu)
 *  ^+  U_mu(x + mu)  U_nu(x + mu + mu) Udagger_mu(x + nu + mu) Udagger_mu(x + nu) Udagger_nu(x)
 */
inline Matrix3x3 local_rectangles_staple_1(__global const Matrixsu3StorageType * const restrict field, const int n, const int t, const int mu, const int nu )
{
	int4 pos;
	int4 pos2;
	Matrixsu3 tmp;

	//first ingredient
	//1.  U_nu(x + mu) * U_nu(x + mu + nu) * Udagger_mu(x + nu + nu) * Udagger_nu(x + nu) * Udagger_nu(x)
	//calculate the coordinates for the matrices. this is the same as with  the rectangles
	//(x, y) = i + mu (call the site idx "x" "i" for now)
	if(mu == 0) {
		pos.x = get_neighbor_temporal(t);
		pos.y = n;
	} else {
		pos.x = t;
		pos.y = get_neighbor(n, mu);
	}
	//(w,z) = x + nu
	if(nu == 0) {
		pos.z = get_neighbor_temporal(t);
		pos.w = n;
	} else {
		pos.z = t;
		pos.w = get_neighbor(n, nu);
	}
	//(x2, y2) = i + mu + nu = (x,y) + nu
	if(nu == 0) {
		pos2.x = get_neighbor_temporal(pos.x);
		pos2.y = pos.y;
	} else {
		pos2.x = pos.x;
		pos2.y = get_neighbor(pos.y, nu);
	}
	//(w2, z2) = i + nu + nu = (w,z) + nu
	if(nu == 0) {
		pos2.z = get_neighbor_temporal(pos.z);
		pos2.w = pos.w;
	} else {
		pos2.z = pos.z;
		pos2.w = get_neighbor(pos.w, nu);
	}

	tmp = multiply_matrixsu3 (get_matrixsu3(field, pos.y, pos.x, nu), get_matrixsu3(field, pos2.y, pos2.x, nu)      );
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, pos2.w, pos2.z, mu) );
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, pos.w, pos.z, nu) );
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, n, t, nu) );

	return matrix_su3to3x3(tmp);
}

inline Matrix3x3 local_rectangles_staple_2(__global const Matrixsu3StorageType * const restrict field, const int n, const int t, const int mu, const int nu )
{
	int4 pos;
	int4 pos2;
	Matrixsu3 tmp;

	//second ingredient
	//2.  U_mu(x + mu) * Udagger_nu(x - nu + mu + mu) Udagger_mu(x - nu + mu) Udagger_mu(x - nu) * U_nu(x - nu)
	//calculate the coordinates for the matrices. this is the same as with  the rectangles
	//(x, y) = i + mu (call the site idx "x" "i" for now)
	if(mu == 0) {
		pos.x = get_neighbor_temporal(t);
		pos.y = n;
	} else {
		pos.x = t;
		pos.y = get_neighbor_spatial(n, mu);
	}
	//(w,z) = i - nu
	if(nu == 0) {
		pos.z = get_lower_neighbor_temporal(t);
		pos.w = n;
	} else {
		pos.z = t;
		pos.w = get_lower_neighbor_spatial(n, nu);
	}
	//(x2, y2) = i + mu - nu = (x,y) - nu
	if(nu == 0) {
		pos2.x = get_lower_neighbor_temporal(pos.x);
		pos2.y = pos.y;
	} else {
		pos2.x = pos.x;
		pos2.y = get_lower_neighbor_spatial(pos.y, nu);
	}
	//(w2, z2) = i + mu - nu + mu = (x,y) - nu + mu = (x2,y2) + mu
	if(mu == 0) {
		pos2.z = get_neighbor_temporal(pos2.x);
		pos2.w = pos2.y;
	} else {
		pos2.z = pos2.x;
		pos2.w = get_neighbor(pos2.y, mu);
	}

	tmp = multiply_matrixsu3_dagger(get_matrixsu3(field, pos.y, pos.x, mu), get_matrixsu3(field, pos2.w, pos2.z, nu)      );
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, pos2.y, pos2.x, mu) );
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, pos.w, pos.z, mu) );
	tmp = multiply_matrixsu3(tmp, get_matrixsu3(field, pos.w, pos.z, nu) );

	return matrix_su3to3x3(tmp);
}

inline Matrix3x3 local_rectangles_staple_3(__global const Matrixsu3StorageType * const restrict field, const int n, const int t, const int mu, const int nu )
{
	int4 pos;
	int4 pos2;
	Matrixsu3 tmp;

	//3. Udagger_nu(x - nu + mu) * Udagger_mu(x - nu) * Udagger_mu(x - mu - nu) * U_nu(x - mu - nu) * U_mu(x - mu)
	//(x, y) = i - mu (call the site idx "x" "i" for now)
	if(mu == 0) {
		pos.x = get_lower_neighbor_temporal(t);
		pos.y = n;
	} else {
		pos.x = t;
		pos.y = get_lower_neighbor_spatial(n, mu);
	}
	//(w,z) = i - nu
	if(nu == 0) {
		pos.z = get_lower_neighbor_temporal(t);
		pos.w = n;
	} else {
		pos.z = t;
		pos.w = get_lower_neighbor_spatial(n, nu);
	}
	//(x2, y2) = i - mu - nu = (x,y) - nu
	if(nu == 0) {
		pos2.x = get_lower_neighbor_temporal(pos.x);
		pos2.y = pos.y;
	} else {
		pos2.x = pos.x;
		pos2.y = get_lower_neighbor_spatial(pos.y, nu);
	}
	//(w2, z2) = i - nu + mu = (w,z) + mu
	if(mu == 0) {
		pos2.z = get_neighbor_temporal(pos.z);
		pos2.w = pos.w;
	} else {
		pos2.z = pos.z;
		pos2.w = get_neighbor(pos.w, mu);
	}

	tmp = multiply_matrixsu3_dagger_dagger(get_matrixsu3(field, pos2.w, pos2.z, nu), get_matrixsu3(field, pos.w, pos.z, mu)      );
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, pos2.y, pos2.x, mu) );
	tmp = multiply_matrixsu3(tmp, get_matrixsu3(field, pos2.y, pos2.x, nu) );
	tmp = multiply_matrixsu3(tmp, get_matrixsu3(field, pos.y, pos.x, mu) );

	return matrix_su3to3x3(tmp);
}

inline Matrix3x3 local_rectangles_staple_4(__global const Matrixsu3StorageType * const restrict field, const int n, const int t, const int mu, const int nu )
{
	int4 pos;
	int4 pos2;
	Matrixsu3 tmp;

	//4. Udagger_nu(x - nu + mu) * Udagger_nu(x - nu - nu + mu) * Udagger_mu(x - nu - nu) * Ud_nu(x - nu - nu) * U_nu(x - nu)
	//(x, y) = i - nu (call the site idx "x" "i" for now)
	if(nu == 0) {
		pos.x = get_lower_neighbor_temporal(t);
		pos.y = n;
	} else {
		pos.x = t;
		pos.y = get_lower_neighbor_spatial(n, nu);
	}
	//(w,z) = i - nu - nu = (x,y) - nu
	if(nu == 0) {
		pos.z = get_lower_neighbor_temporal(pos.x);
		pos.w = pos.y;
	} else {
		pos.z = pos.x;
		pos.w = get_lower_neighbor_spatial(pos.y, nu);
	}
	//(x2, y2) = i - nu + mu = (x,y) + mu
	if(mu == 0) {
		pos2.x = get_neighbor_temporal(pos.x);
		pos2.y = pos.y;
	} else {
		pos2.x = pos.x;
		pos2.y = get_neighbor_spatial(pos.y, mu);
	}
	//(w2, z2) = i - nu  - nu + mu = (w,z) + mu
	if(mu == 0) {
		pos2.z = get_neighbor_temporal(pos.z);
		pos2.w = pos.w;
	} else {
		pos2.z = pos.z;
		pos2.w = get_neighbor(pos.w, mu);
	}

	tmp = multiply_matrixsu3_dagger_dagger(get_matrixsu3(field, pos2.y, pos2.x, nu), get_matrixsu3(field, pos2.w, pos2.z, nu)      );
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, pos.w, pos.z, mu) );
	tmp = multiply_matrixsu3(tmp, get_matrixsu3(field, pos.w, pos.z, nu) );
	tmp = multiply_matrixsu3(tmp, get_matrixsu3(field, pos.y, pos.x, nu) );

	return matrix_su3to3x3(tmp);
}

inline Matrix3x3 local_rectangles_staple_5(__global const Matrixsu3StorageType * const restrict field, const int n, const int t, const int mu, const int nu )
{
	int4 pos;
	int4 pos2;
	Matrixsu3 tmp;

	//5.    U_nu(x + mu) * Udagger_mu(x + nu) * Udagger_mu(x - mu + nu) * Udagger_nu(x - mu) * U_mu(x - mu)
	//(x, y) = i - mu (call the site idx "x" "i" for now)
	if(mu == 0) {
		pos.x = get_lower_neighbor_temporal(t);
		pos.y = n;
	} else {
		pos.x = t;
		pos.y = get_lower_neighbor_spatial(n, mu);
	}
	//(w,z) = i + nu
	if(nu == 0) {
		pos.z = get_neighbor_temporal(t);
		pos.w = n;
	} else {
		pos.z = t;
		pos.w = get_neighbor_spatial(n, nu);
	}
	//(x2, y2) = i - mu + nu = (x,y) + nu
	if(nu == 0) {
		pos2.x = get_neighbor_temporal(pos.x);
		pos2.y = pos.y;
	} else {
		pos2.x = pos.x;
		pos2.y = get_neighbor_spatial(pos.y, nu);
	}
	//(w2, z2) = i + mu
	if(mu == 0) {
		pos2.z = get_neighbor_temporal(t);
		pos2.w = n;
	} else {
		pos2.z = t;
		pos2.w = get_neighbor(n, mu);
	}

	tmp = multiply_matrixsu3_dagger( get_matrixsu3(field, pos2.w, pos2.z, nu), get_matrixsu3(field, pos.w, pos.z, mu) );
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, pos2.y, pos2.x, mu) );
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, pos.y, pos.x, nu) );
	tmp = multiply_matrixsu3(tmp, get_matrixsu3(field, pos.y, pos.x, mu) );

	return matrix_su3to3x3(tmp);
}

inline Matrix3x3 local_rectangles_staple_6(__global const Matrixsu3StorageType * const restrict field, const int n, const int t, const int mu, const int nu )
{
	int4 pos;
	int4 pos2;
	Matrixsu3 tmp;

	//6.  U_mu(x + mu)  U_nu(x + mu + mu) Udagger_mu(x + nu + mu) Udagger_mu(x + nu) Udagger_nu(x)
	//(x, y) = i + mu (call the site idx "x" "i" for now)
	if(mu == 0) {
		pos.x = get_neighbor_temporal(t);
		pos.y = n;
	} else {
		pos.x = t;
		pos.y = get_neighbor(n, mu);
	}
	//(w,z) = x + nu
	if(nu == 0) {
		pos.z = get_neighbor_temporal(t);
		pos.w = n;
	} else {
		pos.z = t;
		pos.w = get_neighbor(n, nu);
	}
	//(x2, y2) = i + mu + nu = (x,y) + nu
	if(nu == 0) {
		pos2.x = get_neighbor_temporal(pos.x);
		pos2.y = pos.y;
	} else {
		pos2.x = pos.x;
		pos2.y = get_neighbor(pos.y, nu);
	}
	//(w2, z2) = i + mu + mu = (x,y) + mu
	if(mu == 0) {
		pos2.z = get_neighbor_temporal(pos.x);
		pos2.w = pos.y;
	} else {
		pos2.z = pos.x;
		pos2.w = get_neighbor(pos.y, mu);
	}

	tmp = multiply_matrixsu3 (get_matrixsu3(field, pos.y, pos.x, mu), get_matrixsu3(field, pos2.w, pos2.z, nu)      );
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, pos2.y, pos2.x, mu) );
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, pos.w, pos.z, mu) );
	tmp = multiply_matrixsu3_dagger(tmp, get_matrixsu3(field, n, t, nu) );

	return matrix_su3to3x3(tmp);
}

inline Matrix3x3 calc_staple(__global const Matrixsu3StorageType * const restrict field, const int pos, const int t, const int mu_in)
{
	Matrix3x3 staple = zero_matrix3x3();
	//iterate through the three directions other than mu
#ifdef _USEGPU_
#pragma unroll 3 // unroll required for proper register reuse when using newer Catalysts on Cypress
#endif
	for(int i = 0; i < NDIM - 1; i++) {
		int nu = (mu_in + i + 1) % NDIM;
		local_staple(&staple, field, pos, t, mu_in, nu );
	}
	return staple;
}

//this is the staple only in the spatial directions only
inline Matrix3x3 calc_staple_sigma (__global const Matrixsu3StorageType * const restrict field, const int pos, const int t, const int mu_in)
{
	Matrix3x3 staple = zero_matrix3x3();
	//iterate through the three directions other than mu
#ifdef _USEGPU_
#pragma unroll 3 // unroll required for proper register reuse when using newer Catalysts on Cypress
#endif
	for(int i = 0; i < NDIM - 1; i++) {
		int nu = (mu_in + i + 1) % NDIM;
		if(nu != 0) {
			local_staple(&staple, field, pos, t, mu_in, nu );
		}
	}
	return staple;
}

//this is the staple only in temporal direction only
inline Matrix3x3 calc_staple_tau (__global const Matrixsu3StorageType * const restrict field, const int pos, const int t, const int mu_in)
{
	int nu = 0;
	Matrix3x3 staple = zero_matrix3x3();
	local_staple(&staple, field, pos, t, mu_in, nu );
	return staple;
}

//this is the rectangles staple
inline Matrix3x3 calc_rectangles_staple(__global const Matrixsu3StorageType * const restrict field, const int pos, const int t, const int mu_in)
{
	Matrix3x3 staple = zero_matrix3x3();
	//iterate through the three directions other than mu
#ifdef _USEGPU_
#pragma unroll 3 // unroll required for proper register reuse when using newer Catalysts on Cypress
#endif
	for(int i = 0; i < NDIM - 1 ; i++) {
		int nu = (mu_in + i + 1) % NDIM;
		staple = add_matrix3x3(staple, local_rectangles_staple_1(field, pos, t, mu_in, nu ));
		staple = add_matrix3x3(staple, local_rectangles_staple_2(field, pos, t, mu_in, nu ));
		staple = add_matrix3x3(staple, local_rectangles_staple_3(field, pos, t, mu_in, nu ));
		staple = add_matrix3x3(staple, local_rectangles_staple_4(field, pos, t, mu_in, nu ));
		staple = add_matrix3x3(staple, local_rectangles_staple_5(field, pos, t, mu_in, nu ));
		staple = add_matrix3x3(staple, local_rectangles_staple_6(field, pos, t, mu_in, nu ));
	}

	return staple;
}
