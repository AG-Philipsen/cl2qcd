#include "host_operations_gaugefield.h"

void copy_to_ocl_format(ocl_s_gaugefield* host_gaugefield, Matrixsu3* gaugefield, const inputparameters * const params)
{
	for(int spacepos = 0; spacepos < NSPACE * NSPACE * NSPACE; spacepos++) {
		for(int t = 0; t < NTIME; t++) {
			for(int mu = 0; mu < NDIM; mu++) {
				const size_t index = get_global_link_pos(mu, spacepos, t, params);
				host_gaugefield[index] = gaugefield[index];
			}
		}
	}
	return;
}

void copy_from_ocl_format(Matrixsu3* gaugefield, ocl_s_gaugefield* host_gaugefield, const inputparameters * const params)
{
	for(int spacepos = 0; spacepos < NSPACE * NSPACE * NSPACE; spacepos++) {
		for(int t = 0; t < NTIME; t++) {
			for(int mu = 0; mu < NDIM; mu++) {
				const size_t index = get_global_link_pos(mu, spacepos, t, params);
				gaugefield[index] = host_gaugefield[index];
			}
		}
	}
	return;
}

void set_gaugefield_cold(hmc_complex * field)
{
	for(int t = 0; t < NTIME; t++) {
		for(int n = 0; n < VOLSPACE; n++) {
			for(int mu = 0; mu < NDIM; mu++) {
				hmc_su3matrix tmp;
				unit_su3matrix(&tmp);
				put_su3matrix(field, &tmp, n, t, mu);
			}
		}
	}
	return;
}

void set_gaugefield_hot(hmc_complex * field)
{
	for(int t = 0; t < NTIME; t++) {
		for(int n = 0; n < VOLSPACE; n++) {
			for(int mu = 0; mu < NDIM; mu++) {
				hmc_su3matrix tmp;
				random_su3matrix(&tmp);
				put_su3matrix(field, &tmp, n, t, mu);
			}
		}
	}
	return;
}

void copy_gaugefield_from_ildg_format(hmc_complex * gaugefield, hmc_float * gaugefield_tmp, int check, const inputparameters * const params)
{
	//little check if arrays are big enough
	if (params->get_vol4d() *NDIM*NC*NC * 2 != check) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values!!\nCheck global settings!!";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}

	int cter = 0;
	//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
	for (int t = 0; t < NTIME; t++) {
		//if the alg is known to be correct, the next three for-loops could also be changed to one
		for (int i = 0; i < NSPACE; i++) {
			for (int j = 0; j < NSPACE; j++) {
				for (int k = 0; k < NSPACE; k++) {
					for (int l = 0; l < NDIM; l++) {
						int spacepos = k + j * NSPACE + i * NSPACE * NSPACE;
						int globalpos = l + spacepos * NDIM + t * VOLSPACE * NDIM;
						hmc_su3matrix tmp;
#ifdef _RECONSTRUCT_TWELVE_
						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
								int ncindex = m + (NC - 1) * n;
								//ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
								//which is stored in one single array here
								//skip NC*NC*2 cmplx numbers
								int pos = 2 * n + 2 * m * NC + globalpos * NC * NC * 2;
								if(m < NC - 1) {
									tmp[ncindex].re = gaugefield_tmp[pos];
									tmp[ncindex].im = gaugefield_tmp[pos + 1];
								}
								cter++;
							}
						}
#else
						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
								//ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
								//which is stored in one single array here
								//skip NC*NC*2 cmplx numbers
								int pos = 2 * n + 2 * m * NC + globalpos * NC * NC * 2;
								tmp[m][n].re = gaugefield_tmp[pos];
								tmp[m][n].im = gaugefield_tmp[pos + 1];
								cter++;
							}
						}
#endif
						put_su3matrix(gaugefield, &tmp, spacepos, t, (l + 1) % NDIM);
					}
				}
			}
		}
	}

	if(cter * 2 != check) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values! there were " << cter * 2 << " vals set and not " << check << ".";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}

	return;
}

void copy_gaugefield_to_ildg_format(hmc_float * dest, hmc_complex * source)
{

	int cter = 0;
	//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
	for (int t = 0; t < NTIME; t++) {
		//if the alg is known to be correct, the next three for-loops could also be changed to one
		for (int i = 0; i < NSPACE; i++) {
			for (int j = 0; j < NSPACE; j++) {
				for (int k = 0; k < NSPACE; k++) {
					for (int l = 0; l < NDIM; l++) {
						int spacepos = k + j * NSPACE + i * NSPACE * NSPACE;
						int globalpos = l + spacepos * NDIM + t * VOLSPACE * NDIM;
						hmc_su3matrix tmp;
						get_su3matrix(&tmp, source, spacepos, t, (l + 1) % NDIM);
#ifdef _RECONSTRUCT_TWELVE_
						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
								int ncindex = m + (NC - 1) * n;
								//ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
								//which is stored in one single array here
								//skip NC*NC*2 cmplx numbers

								int pos = 2 * n + 2 * m * NC + globalpos * NC * NC * 2;
								if(m < NC - 1) {


									dest[pos]     = tmp[ncindex].re;
									dest[pos + 1] = tmp[ncindex].im;
								}
								if (m == NC - 1) {
									hmc_complex tmpc = reconstruct_su3 (&tmp, n);
									dest[pos]     = tmpc.re;
									dest[pos + 1] = tmpc.im;
								}

								cter++;
							}
						}
#else
						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
								//ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
								//which is stored in one single array here
								//skip NC*NC*2 cmplx numbers
								int pos = 2 * n + 2 * m * NC + globalpos * NC * NC * 2;
								dest[pos]     = tmp[m][n].re;
								dest[pos + 1] = tmp[m][n].im;
								cter++;
							}
						}
#endif
					}
				}
			}
		}
	}

	return;
}

hmc_complex global_trace_su3(hmc_complex * field, int mu)
{
	hmc_complex sum;
	sum.re = 0;
	sum.im = 0;
	for(int t = 0; t < NTIME; t++) {
		for(int n = 0; n < VOLSPACE; n++) {
			hmc_su3matrix tmp;
			get_su3matrix(&tmp, field, n, t, mu);
			sum.re += trace_su3matrix(&tmp).re;
			sum.im += trace_su3matrix(&tmp).im;
		}
	}
	return sum;
}




void get_su3matrix(hmc_su3matrix * out, hmc_complex * in, int spacepos, int timepos, int mu)
{
#ifdef _RECONSTRUCT_TWELVE_
	for(int n = 0; n < NC*(NC - 1); n++) (*out)[n] = in[get_hmc_gaugefield_index(n, spacepos, timepos, mu)];
#else
	for(int a = 0; a < NC; a++) {
		for(int b = 0; b < NC; b++) {
			(*out)[a][b] = in[get_hmc_gaugefield_index(a, b, spacepos, timepos, mu)];
		}
	}
#endif
	return;
}

void put_su3matrix(hmc_complex * field, hmc_su3matrix * in, int spacepos, int timepos, int mu)
{
#ifdef _RECONSTRUCT_TWELVE_
	for(int n = 0; n < NC*(NC - 1); n++) field[get_hmc_gaugefield_index(n, spacepos, timepos, mu)] = (*in)[n];
#else
	for(int a = 0; a < NC; a++) {
		for(int b = 0; b < NC; b++) {
			size_t index = get_hmc_gaugefield_index(a, b, spacepos, timepos, mu);
			field[index] = (*in)[a][b];
		}
	}
#endif
	return;
}

void local_polyakov(hmc_complex * field, hmc_su3matrix * prod, int n)
{
	unit_su3matrix(prod);
	for(int t = 0; t < NTIME; t++) {
		hmc_su3matrix tmp;
		get_su3matrix(&tmp, field, n, t, 0);
		accumulate_su3matrix_prod(prod, &tmp);
	}
	return;
}

void local_plaquette(hmc_complex * field, hmc_su3matrix * prod, int n, int t, int mu, int nu, const inputparameters * const params )
{
	hmc_su3matrix tmp;
	//u_mu(x)
	get_su3matrix(prod, field, n, t, mu);
	//u_nu(x+mu)
	if(mu == 0) {
		int newt = (t + 1) % NTIME; //(haha)
		get_su3matrix(&tmp, field, n, newt, nu);
	} else {
		get_su3matrix(&tmp, field, get_neighbor(n, mu, params), t, nu);
	}
	accumulate_su3matrix_prod(prod, &tmp);
	//adjoint(u_mu(x+nu))
	if(nu == 0) {
		int newt = (t + 1) % NTIME;
		get_su3matrix(&tmp, field, n, newt, mu);
	} else {
		get_su3matrix(&tmp, field, get_neighbor(n, nu, params), t, mu);
	}
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(prod, &tmp);
	//adjoint(u_nu(x))
	get_su3matrix(&tmp, field, n, t, nu);
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(prod, &tmp);

	return;
}

/** @todo memcpy ... */
void copy_gaugefield(hmc_complex * source, hmc_complex * dest, const inputparameters * const params)
{
	// copies source to destination within cpu memory, layer for gaugefield array
	return complexcopy(source, dest, params->get_gaugefieldsize()); // SL: not tested
}

void local_Q_plaquette(hmc_3x3matrix * out, hmc_complex * field, int n, int t, int mu, int nu, const inputparameters * const params )
{
	hmc_su3matrix tmp;
	int newpos;

	//1st plaq
	hmc_su3matrix plaq1;
	//u_mu(x)
	get_su3matrix(&plaq1, field, n, t, mu);
	//u_nu(x+mu)
	if(mu == 0) {
		int newt = (t + 1) % NTIME;
		get_su3matrix(&tmp, field, n, newt, nu);
	} else
		get_su3matrix(&tmp, field, get_neighbor(n, mu, params), t, nu);
	accumulate_su3matrix_prod(&plaq1, &tmp);
	//adjoint(u_mu(x+nu))
	if(nu == 0) {
		int newt = (t + 1) % NTIME;
		get_su3matrix(&tmp, field, n, newt, mu);
	} else
		get_su3matrix(&tmp, field, get_neighbor(n, nu, params), t, mu);
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(&plaq1, &tmp);
	//adjoint(u_nu(x))
	get_su3matrix(&tmp, field, n, t, nu);
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(&plaq1, &tmp);

	//2nd plaq
	hmc_su3matrix plaq2;
	//U_nu(x)
	get_su3matrix(&plaq2, field, n, t, nu);
	//adj (u_mu(x-mu+nu))
	newpos = get_lower_neighbor(n, mu, params);
	if (nu == 0) {
		int newt =  (t + 1) % NTIME;
		get_su3matrix(&tmp, field, newpos, newt, mu);
	} else if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		get_su3matrix(&tmp, field, get_neighbor(n, nu, params), newt, mu);
	} else
		get_su3matrix(&tmp, field, get_neighbor(newpos, nu, params), t, mu);
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(&plaq2, &tmp);
	//adj (u_nu(x-mu))
	if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		get_su3matrix(&tmp, field, n, newt, nu);
	} else
		get_su3matrix(&tmp, field, get_lower_neighbor(n, mu, params), t, nu);
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(&plaq2, &tmp);
	//u_mu(x-mu)
	if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		get_su3matrix(&tmp, field, n, newt, mu);
	} else
		get_su3matrix(&tmp, field, get_lower_neighbor(n, mu, params), t, mu);
	accumulate_su3matrix_prod(&plaq2, &tmp);

	//3rd plaq
	hmc_su3matrix plaq3;
	//adj (u_mu(x-mu))
	if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		get_su3matrix(&tmp, field, n, newt, mu);
	} else
		get_su3matrix(&tmp, field, get_lower_neighbor(n, mu, params), t, mu);
	adjoin_su3matrix(&tmp);
	copy_su3matrix(&plaq3, &tmp);
	//adj (u_nu(x-mu-nu))
	newpos = get_lower_neighbor(n, mu, params);
	if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		get_su3matrix(&tmp, field, newpos, newt, nu);
	} else if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		get_su3matrix(&tmp, field, get_lower_neighbor(n, nu, params), newt, nu);
	} else
		get_su3matrix(&tmp, field, get_lower_neighbor(newpos, nu, params), t, nu);
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(&plaq3, &tmp);
	//u_mu(x-mu-nu)
	newpos = get_lower_neighbor(n, mu, params);
	if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		get_su3matrix(&tmp, field, newpos, newt, mu);
	} else if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		get_su3matrix(&tmp, field, get_lower_neighbor(n, nu, params), newt, mu);
	} else
		get_su3matrix(&tmp, field, get_lower_neighbor(newpos, nu, params), t, mu);
	accumulate_su3matrix_prod(&plaq3, &tmp);
	//u_nu(x-nu)
	if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		get_su3matrix(&tmp, field, n, newt, nu);
	} else
		get_su3matrix(&tmp, field, get_lower_neighbor(n, nu, params), t, nu);
	accumulate_su3matrix_prod(&plaq3, &tmp);

	//4th plaq
	hmc_su3matrix plaq4;
	//adj(u_nu(x-nu))
	adjoin_su3matrix(&tmp);
	copy_su3matrix(&plaq4, &tmp);
	//u_mu(x-nu)
	if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		get_su3matrix(&tmp, field, n, newt, mu);
	} else
		get_su3matrix(&tmp, field, get_lower_neighbor(n, nu, params), t, mu);
	accumulate_su3matrix_prod(&plaq4, &tmp);
	//u_nu(x+mu-nu)
	newpos = get_lower_neighbor(n, nu, params);
	if (mu == 0) {
		int newt =  (t + 1) % NTIME;
		get_su3matrix(&tmp, field, newpos, newt, nu);
	} else if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		get_su3matrix(&tmp, field, get_neighbor(n, mu, params), newt, nu);
	} else
		get_su3matrix(&tmp, field, get_neighbor(newpos, mu, params), t, nu);
	accumulate_su3matrix_prod(&plaq4, &tmp);
	//adj (u_mu(x))
	get_su3matrix(&tmp, field, n, t, mu);
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(&plaq4, &tmp);

	//Sum up
	su3matrix_to_3x3matrix (out, &plaq1);
	accumulate_su3matrix_3x3_add(out, &plaq2);
	accumulate_su3matrix_3x3_add(out, &plaq3);
	accumulate_su3matrix_3x3_add(out, &plaq4);
}

size_t get_hmc_gaugefield_index(size_t spacepos, size_t timepos, size_t mu)
{
	size_t result = (mu * VOLSPACE + spacepos ) * NTIME + timepos;

#ifdef _RECONSTRUCT_TWELVE_
	result *= NC * (NC - 1);
#else
	result *= NC * NC;
#endif
	return result;
}
#ifdef _RECONSTRUCT_TWELVE_
size_t get_hmc_gaugefield_index(size_t ncolor, size_t spacepos, size_t timepos, size_t mu)
{
	size_t result = (mu * VOLSPACE + spacepos ) * NTIME + timepos;
	result += ncolor * NDIM * VOLSPACE * NTIME;
	return result;
}
#else
size_t get_hmc_gaugefield_index(size_t m, size_t n, size_t spacepos, size_t timepos, size_t mu)
{
	size_t result = (mu * VOLSPACE + spacepos ) * NTIME + timepos;
	result += (m * NC + n) * NDIM * VOLSPACE * NTIME;
	return result;
}
#endif
