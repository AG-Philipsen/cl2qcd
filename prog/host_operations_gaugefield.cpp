#include "host_operations_gaugefield.h"

//taken from the opencl-file

Matrixsu3 multiply_matrixsu3(const Matrixsu3 p, const Matrixsu3 q)
{
  Matrixsu3 out;

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

  return out;
}

Matrixsu3 multiply_matrixsu3_dagger(const Matrixsu3 p, const Matrixsu3 q)
{
  Matrixsu3 out;

  out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e01.re + p.e02.re * q.e02.re
    + p.e00.im * q.e00.im + p.e01.im * q.e01.im + p.e02.im * q.e02.im;
  out.e00.im = -p.e00.re * q.e00.im - p.e01.re * q.e01.im - p.e02.re * q.e02.im
    + p.e00.im * q.e00.re + p.e01.im * q.e01.re + p.e02.im * q.e02.re;

  out.e01.re = p.e00.re * q.e10.re + p.e01.re * q.e11.re + p.e02.re * q.e12.re
    + p.e00.im * q.e10.im + p.e01.im * q.e11.im + p.e02.im * q.e12.im;
  out.e01.im = -p.e00.re * q.e10.im - p.e01.re * q.e11.im - p.e02.re * q.e12.im
    + p.e00.im * q.e10.re + p.e01.im * q.e11.re + p.e02.im * q.e12.re;

  out.e02.re = p.e00.re * q.e20.re + p.e01.re * q.e21.re + p.e02.re * q.e22.re
    + p.e00.im * q.e20.im + p.e01.im * q.e21.im + p.e02.im * q.e22.im;
  out.e02.im = -p.e00.re * q.e20.im - p.e01.re * q.e21.im - p.e02.re * q.e22.im
    + p.e00.im * q.e20.re + p.e01.im * q.e21.re + p.e02.im * q.e22.re;

  out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e01.re + p.e12.re * q.e02.re
    + p.e10.im * q.e00.im + p.e11.im * q.e01.im + p.e12.im * q.e02.im;
  out.e10.im = -p.e10.re * q.e00.im - p.e11.re * q.e01.im - p.e12.re * q.e02.im
    + p.e10.im * q.e00.re + p.e11.im * q.e01.re + p.e12.im * q.e02.re;

  out.e11.re = p.e10.re * q.e10.re + p.e11.re * q.e11.re + p.e12.re * q.e12.re
    + p.e10.im * q.e10.im + p.e11.im * q.e11.im + p.e12.im * q.e12.im;
  out.e11.im = -p.e10.re * q.e10.im - p.e11.re * q.e11.im - p.e12.re * q.e12.im
    + p.e10.im * q.e10.re + p.e11.im * q.e11.re + p.e12.im * q.e12.re;

  out.e12.re = p.e10.re * q.e20.re + p.e11.re * q.e21.re + p.e12.re * q.e22.re
    + p.e10.im * q.e20.im + p.e11.im * q.e21.im + p.e12.im * q.e22.im;
  out.e12.im = -p.e10.re * q.e20.im - p.e11.re * q.e21.im - p.e12.re * q.e22.im
    + p.e10.im * q.e20.re + p.e11.im * q.e21.re + p.e12.im * q.e22.re;

  out.e20.re = p.e20.re * q.e00.re + p.e21.re * q.e01.re + p.e22.re * q.e02.re
    + p.e20.im * q.e00.im + p.e21.im * q.e01.im + p.e22.im * q.e02.im;
  out.e20.im = -p.e20.re * q.e00.im - p.e21.re * q.e01.im - p.e22.re * q.e02.im
    + p.e20.im * q.e00.re + p.e21.im * q.e01.re + p.e22.im * q.e02.re;

  out.e21.re = p.e20.re * q.e10.re + p.e21.re * q.e11.re + p.e22.re * q.e12.re
    + p.e20.im * q.e10.im + p.e21.im * q.e11.im + p.e22.im * q.e12.im;
  out.e21.im = -p.e20.re * q.e10.im - p.e21.re * q.e11.im - p.e22.re * q.e12.im
    + p.e20.im * q.e10.re + p.e21.im * q.e11.re + p.e22.im * q.e12.re;

  out.e22.re = p.e20.re * q.e20.re + p.e21.re * q.e21.re + p.e22.re * q.e22.re
    + p.e20.im * q.e20.im + p.e21.im * q.e21.im + p.e22.im * q.e22.im;
  out.e22.im = -p.e20.re * q.e20.im - p.e21.re * q.e21.im - p.e22.re * q.e22.im
    + p.e20.im * q.e20.re + p.e21.im * q.e21.re + p.e22.im * q.e22.re;

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

  out.e20.re = 0.;
  out.e20.im = 0.;
  out.e21.re = 0.;
  out.e21.im = 0.;
  out.e22.re = 1.;
  out.e22.im = 0.;

  return out;
}



Matrixsu3 get_matrixsu3(Matrixsu3 * in, int spacepos, int timepos, int mu, const inputparameters * const parameters){
	Matrixsu3 tmp;
	size_t link_pos = get_global_link_pos(mu, spacepos, timepos, parameters);
	tmp = in[link_pos];
	return tmp;
}

void put_matrixsu3(Matrixsu3 * field, Matrixsu3 in, int spacepos, int timepos, int mu, const inputparameters * const parameters){
  size_t link_pos = get_global_link_pos(mu, spacepos, timepos, parameters);
	field[link_pos] = in;
}

void get_su3matrix(hmc_su3matrix * out, hmc_complex * in, int spacepos, int timepos, int mu, const inputparameters * const parameters)
{
	for(int a = 0; a < NC; a++) {
		for(int b = 0; b < NC; b++) {
			(*out)[a][b] = in[get_hmc_gaugefield_index(a, b, spacepos, timepos, mu, parameters)];
		}
	}
	return;
}

void put_su3matrix(hmc_complex * field, hmc_su3matrix * in, int spacepos, int timepos, int mu, const inputparameters * const parameters)
{
	for(int a = 0; a < NC; a++) {
		for(int b = 0; b < NC; b++) {
			size_t index = get_hmc_gaugefield_index(a, b, spacepos, timepos, mu, parameters);
			field[index] = (*in)[a][b];
		}
	}
	return;
}

Matrixsu3 local_polyakov(Matrixsu3 * field, int n, const inputparameters * const parameters)
{
        Matrixsu3 res, prod;
	res = unit_matrixsu3();
	for(int t = 0; t < parameters->get_nt(); t++) {
		prod = get_matrixsu3(field, n, t, 0, parameters);
		res = multiply_matrixsu3(res, prod);
	}
	return res;
}

//CP: I introduced explicit calculations of the neighbors because this does not rely on any geometric conventions!
Matrixsu3 local_plaquette(Matrixsu3 * field, int coord_in[NDIM], int mu, int nu, const inputparameters * const parameters)
{
        Matrixsu3 res, tmp;
	//coordinates of neighbors
	int coord[NDIM];
	coord[0] = coord_in[0];
	coord[1] = coord_in[1];
	coord[2] = coord_in[2];
	coord[3] = coord_in[3];
	//spatial index
	int n;
	
	const size_t NTIME = parameters->get_nt();
	const size_t NSPACE = parameters->get_ns();
	//u_mu(x)
	n = get_nspace(coord_in, parameters);
	res = get_matrixsu3(field, n, coord[0], mu, parameters);
	//u_nu(x+mu)
	if(mu == 0) {
		coord[mu] = (coord_in[mu] + 1) % NTIME;
		n = get_nspace(coord, parameters);
		tmp = get_matrixsu3(field, n, coord[mu], nu, parameters);
		coord[mu] = coord_in[mu];
	} else {
		coord[mu] = (coord_in[mu] + 1) % NSPACE;
		int newn = get_nspace(coord, parameters);
		tmp = get_matrixsu3(field, newn, coord[0], nu, parameters);
		coord[mu] = coord_in[mu];
	}
	//accumulate_su3matrix_prod(&prod, &tmp);
	res = multiply_matrixsu3(res, tmp);
	//adjoint(u_mu(x+nu))
	if(nu == 0) {
		coord[nu] = (coord_in[nu] + 1) % NTIME;
		n = get_nspace(coord, parameters);
	        tmp = get_matrixsu3(field, n, coord[nu], mu, parameters);
		coord[nu] = coord_in[nu];
	} else {
		coord[nu] = (coord_in[nu] + 1) % NSPACE;
		int newn = get_nspace(coord, parameters);
	        tmp = get_matrixsu3(field, newn, coord[0], mu, parameters);
		coord[nu] = coord_in[nu];
	}
	res = multiply_matrixsu3_dagger(res, tmp);

	//adjoint(u_nu(x))
	n = get_nspace(coord_in, parameters);
        tmp = get_matrixsu3(field, n, coord[0], nu, parameters);
	res = multiply_matrixsu3_dagger(res, tmp);

	return res;
}
