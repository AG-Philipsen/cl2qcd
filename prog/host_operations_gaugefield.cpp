#include "host_operations_gaugefield.h"

void set_gaugefield_cold(hmc_complex * field, const inputparameters * const parameters)
{
	for(int t = 0; t < parameters->get_nt(); t++) {
		for(int n = 0; n < parameters->get_volspace(); n++) {
			for(int mu = 0; mu < NDIM; mu++) {
				hmc_su3matrix tmp;
				unit_su3matrix(&tmp);
				put_su3matrix(field, &tmp, n, t, mu, parameters);
			}
		}
	}
	return;
}

void set_gaugefield_hot(hmc_complex * field, const inputparameters * const parameters)
{
	for(int t = 0; t < parameters->get_nt(); t++) {
		for(int n = 0; n < parameters->get_volspace(); n++) {
			for(int mu = 0; mu < NDIM; mu++) {
				hmc_su3matrix tmp;
				random_su3matrix(&tmp);
				put_su3matrix(field, &tmp, n, t, mu, parameters);
			}
		}
	}
	return;
}

Matrixsu3 get_matrixsu3(Matrixsu3 * in, int spacepos, int timepos, int mu, const inputparameters * const parameters){
	Matrixsu3 tmp;
	size_t link_pos = get_global_link_pos(spacepos, mu, timepos, parameters);
	tmp = in[link_pos];
	return tmp;
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

Matrixsu3 local_plaquette(Matrixsu3 * field, int n, int t, int mu, int nu, const inputparameters * const parameters)
{
	Matrixsu3 res;
	
	//using the old methods
	hmc_su3matrix prod, tmp;
	const size_t NTIME = parameters->get_nt();
	//u_mu(x)
	get_su3matrix_tmp(&prod, field, n, t, mu, parameters);
	//u_nu(x+mu)
	if(mu == 0) {
		int newt = (t + 1) % NTIME; //(haha)
		get_su3matrix_tmp(&tmp, field, n, newt, nu, parameters);
	} else {
		get_su3matrix_tmp(&tmp, field, get_neighbor(n, mu, parameters), t, nu, parameters);
	}
	accumulate_su3matrix_prod(&prod, &tmp);
	//adjoint(u_mu(x+nu))
	if(nu == 0) {
		int newt = (t + 1) % NTIME;
		get_su3matrix_tmp(&tmp, field, n, newt, mu, parameters);
	} else {
		get_su3matrix_tmp(&tmp, field, get_neighbor(n, nu, parameters), t, mu, parameters);
	}
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(&prod, &tmp);
	//adjoint(u_nu(x))
	get_su3matrix_tmp(&tmp, field, n, t, nu, parameters);
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(&prod, &tmp);

	res = convert_hmc_matrixsu3_to_Matrixsu3(prod);	
	return res;
}

size_t get_hmc_gaugefield_index(size_t m, size_t n, size_t spacepos, size_t timepos, size_t mu, const inputparameters * const parameters)
{
	const size_t VOLSPACE = parameters->get_volspace();
	const size_t NTIME = parameters->get_nt();
	size_t result = (mu * VOLSPACE + spacepos ) * NTIME + timepos;
	result += (m * NC + n) * NDIM * VOLSPACE * NTIME;
	return result;
}
