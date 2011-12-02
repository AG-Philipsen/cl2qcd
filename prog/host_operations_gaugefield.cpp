#include "host_operations_gaugefield.h"

void copy_to_ocl_format(ocl_s_gaugefield* host_gaugefield, Matrixsu3* gaugefield, const inputparameters * const parameters)
{
	const size_t NSPACE = parameters->get_ns();
	const size_t NTIME = parameters->get_nt();
	for(size_t spacepos = 0; spacepos < NSPACE * NSPACE * NSPACE; spacepos++) {
		for(size_t t = 0; t < NTIME; t++) {
			for(int mu = 0; mu < NDIM; mu++) {
				const size_t index = get_global_link_pos(mu, spacepos, t, parameters);
				host_gaugefield[index] = gaugefield[index];
			}
		}
	}
	return;
}

void copy_from_ocl_format(Matrixsu3* gaugefield, ocl_s_gaugefield* host_gaugefield, const inputparameters * const parameters)
{
	const size_t NSPACE = parameters->get_ns();
	const size_t NTIME = parameters->get_nt();
	for(size_t spacepos = 0; spacepos < NSPACE * NSPACE * NSPACE; spacepos++) {
		for(size_t t = 0; t < NTIME; t++) {
			for(int mu = 0; mu < NDIM; mu++) {
				const size_t index = get_global_link_pos(mu, spacepos, t, parameters);
				gaugefield[index] = host_gaugefield[index];
			}
		}
	}
	return;
}

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

void copy_gaugefield_from_ildg_format(hmc_complex * gaugefield, hmc_float * gaugefield_tmp, int check, const inputparameters * const parameters)
{
	//little check if arrays are big enough
	if (parameters->get_vol4d() *NDIM*NC*NC * 2 != check) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values!!\nCheck global settings!!";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}

	const size_t VOLSPACE = parameters->get_volspace();
	const size_t NSPACE = parameters->get_ns();

	int cter = 0;
	//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
	for (int t = 0; t < parameters->get_nt(); t++) {
		//if the alg is known to be correct, the next three for-loops could also be changed to one
		for (size_t i = 0; i < NSPACE; i++) {
			for (size_t j = 0; j < NSPACE; j++) {
				for (size_t k = 0; k < NSPACE; k++) {
					for (int l = 0; l < NDIM; l++) {
						int spacepos = k + j * NSPACE + i * NSPACE * NSPACE;
						int globalpos = l + spacepos * NDIM + t * VOLSPACE * NDIM;
						hmc_su3matrix tmp;
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
						put_su3matrix(gaugefield, &tmp, spacepos, t, (l + 1) % NDIM, parameters);
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

void copy_gaugefield_to_ildg_format(hmc_float * dest, hmc_complex * source, const inputparameters * const parameters)
{

	int cter = 0;
	const size_t NSPACE = parameters->get_ns();

	//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
	for (int t = 0; t < parameters->get_nt(); t++) {
		//if the alg is known to be correct, the next three for-loops could also be changed to one
		for (size_t i = 0; i < NSPACE; i++) {
			for (size_t j = 0; j < NSPACE; j++) {
				for (size_t k = 0; k < NSPACE; k++) {
					for (int l = 0; l < NDIM; l++) {
						int spacepos = k + j * NSPACE + i * NSPACE * NSPACE;
						int globalpos = l + spacepos * NDIM + t * parameters->get_volspace() * NDIM;
						hmc_su3matrix tmp;
						get_su3matrix(&tmp, source, spacepos, t, (l + 1) % NDIM, parameters);
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
					}
				}
			}
		}
	}

	return;
}

hmc_complex global_trace_su3(hmc_complex * field, int mu, const inputparameters * const parameters)
{
	hmc_complex sum;
	sum.re = 0;
	sum.im = 0;
	for(int t = 0; t < parameters->get_nt(); t++) {
		for(int n = 0; n < parameters->get_volspace(); n++) {
			hmc_su3matrix tmp;
			get_su3matrix(&tmp, field, n, t, mu, parameters);
			sum.re += trace_su3matrix(&tmp).re;
			sum.im += trace_su3matrix(&tmp).im;
		}
	}
	return sum;
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

/** @todo memcpy ... */
void copy_gaugefield(hmc_complex * source, hmc_complex * dest, const inputparameters * const parameters)
{
	// copies source to destination within cpu memory, layer for gaugefield array
	return complexcopy(source, dest, parameters->get_gaugefieldsize()); // SL: not tested
}

size_t get_hmc_gaugefield_index(size_t m, size_t n, size_t spacepos, size_t timepos, size_t mu, const inputparameters * const parameters)
{
	const size_t VOLSPACE = parameters->get_volspace();
	const size_t NTIME = parameters->get_nt();
	size_t result = (mu * VOLSPACE + spacepos ) * NTIME + timepos;
	result += (m * NC + n) * NDIM * VOLSPACE * NTIME;
	return result;
}
