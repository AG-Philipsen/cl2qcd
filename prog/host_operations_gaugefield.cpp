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

void local_polyakov(hmc_complex * field, hmc_su3matrix * prod, int n, const inputparameters * const parameters)
{
	unit_su3matrix(prod);
	for(int t = 0; t < parameters->get_nt(); t++) {
		hmc_su3matrix tmp;
		get_su3matrix(&tmp, field, n, t, 0, parameters);
		accumulate_su3matrix_prod(prod, &tmp);
	}
	return;
}

void local_plaquette(hmc_complex * field, hmc_su3matrix * prod, int n, int t, int mu, int nu, const inputparameters * const parameters)
{
	const size_t NTIME = parameters->get_nt();
	hmc_su3matrix tmp;
	//u_mu(x)
	get_su3matrix(prod, field, n, t, mu, parameters);
	//u_nu(x+mu)
	if(mu == 0) {
		int newt = (t + 1) % NTIME; //(haha)
		get_su3matrix(&tmp, field, n, newt, nu, parameters);
	} else {
		get_su3matrix(&tmp, field, get_neighbor(n, mu, parameters), t, nu, parameters);
	}
	accumulate_su3matrix_prod(prod, &tmp);
	//adjoint(u_mu(x+nu))
	if(nu == 0) {
		int newt = (t + 1) % NTIME;
		get_su3matrix(&tmp, field, n, newt, mu, parameters);
	} else {
		get_su3matrix(&tmp, field, get_neighbor(n, nu, parameters), t, mu, parameters);
	}
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(prod, &tmp);
	//adjoint(u_nu(x))
	get_su3matrix(&tmp, field, n, t, nu, parameters);
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(prod, &tmp);

	return;
}

size_t get_hmc_gaugefield_index(size_t m, size_t n, size_t spacepos, size_t timepos, size_t mu, const inputparameters * const parameters)
{
	const size_t VOLSPACE = parameters->get_volspace();
	const size_t NTIME = parameters->get_nt();
	size_t result = (mu * VOLSPACE + spacepos ) * NTIME + timepos;
	result += (m * NC + n) * NDIM * VOLSPACE * NTIME;
	return result;
}
