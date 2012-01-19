#include "host_operations_gaugefield.h"

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

//this is a temporal thing...
void get_su3matrix_tmp(hmc_su3matrix * out, Matrixsu3 * in, int spacepos, int timepos, int mu, const inputparameters * const parameters)
{
	Matrixsu3 tmp = get_matrixsu3(in, spacepos, timepos, mu, parameters);
	(*out)[0][0] = tmp.e00;
	(*out)[0][1] = tmp.e01;
	(*out)[0][2] = tmp.e02;
	(*out)[1][0] = tmp.e10;
	(*out)[1][1] = tmp.e11;
	(*out)[1][2] = tmp.e12;
	(*out)[2][0] = tmp.e20;
	(*out)[2][1] = tmp.e21;
	(*out)[2][2] = tmp.e22;
	return;
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

//temporal:
Matrixsu3 convert_hmc_matrixsu3_to_Matrixsu3(hmc_su3matrix in){
	Matrixsu3 res;
	res.e00 = in[0][0];
	res.e01 = in[0][1];
	res.e02 = in[0][2];
	res.e10 = in[1][0];
	res.e11 = in[1][1];
	res.e12 = in[1][2];
	res.e20 = in[2][0];
	res.e21 = in[2][1];
	res.e22 = in[2][2];
	return res;
}

Matrixsu3 local_polyakov(Matrixsu3 * field, int n, const inputparameters * const parameters)
{
	Matrixsu3 res;
	
	//CP: using the old methods...
	hmc_su3matrix prod;
	unit_su3matrix(&prod);
	for(int t = 0; t < parameters->get_nt(); t++) {
		hmc_su3matrix tmp;
		get_su3matrix_tmp(&tmp, field, n, t, 0, parameters);
		accumulate_su3matrix_prod(&prod, &tmp);
	}
	
	res = convert_hmc_matrixsu3_to_Matrixsu3(prod);
	return res;
}

//CP: I introduced explicit calculations of the neighbors because this does not rely on any geometric conventions!
Matrixsu3 local_plaquette(Matrixsu3 * field, int coord_in[NDIM], int mu, int nu, const inputparameters * const parameters)
{
	Matrixsu3 res;
	//coordinates of neighbors
	int coord[NDIM];
	coord[0] = coord_in[0];
	coord[1] = coord_in[1];
	coord[2] = coord_in[2];
	coord[3] = coord_in[3];
	//spatial index
	int n;
	
	//using the old methods
	hmc_su3matrix prod, tmp;
	const size_t NTIME = parameters->get_nt();
	const size_t NSPACE = parameters->get_ns();
	//u_mu(x)
	n = get_nspace(coord_in, parameters);
	get_su3matrix_tmp(&prod, field, n, coord[0], mu, parameters);
	//u_nu(x+mu)
	if(mu == 0) {
		coord[mu] = (coord_in[mu] + 1) % NTIME;
		n = get_nspace(coord, parameters);
		get_su3matrix_tmp(&tmp, field, n, coord[mu], nu, parameters);
		coord[mu] = coord_in[mu];
	} else {
		coord[mu] = (coord_in[mu] + 1) % NSPACE;
		int newn = get_nspace(coord, parameters);
		get_su3matrix_tmp(&tmp, field, newn, coord[0], nu, parameters);
		coord[mu] = coord_in[mu];
	}
	accumulate_su3matrix_prod(&prod, &tmp);
	//adjoint(u_mu(x+nu))
	if(nu == 0) {
		coord[nu] = (coord_in[nu] + 1) % NTIME;
		n = get_nspace(coord, parameters);
		get_su3matrix_tmp(&tmp, field, n, coord[nu], mu, parameters);
		coord[nu] = coord_in[nu];
	} else {
		coord[nu] = (coord_in[nu] + 1) % NSPACE;
		int newn = get_nspace(coord, parameters);
		get_su3matrix_tmp(&tmp, field, newn, coord[0], mu, parameters);
		coord[nu] = coord_in[nu];
	}
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(&prod, &tmp);
	//adjoint(u_nu(x))
	n = get_nspace(coord_in, parameters);
	get_su3matrix_tmp(&tmp, field, n, coord[0], nu, parameters);
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(&prod, &tmp);

	res = convert_hmc_matrixsu3_to_Matrixsu3(prod);	
	return res;
}
