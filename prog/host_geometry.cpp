#include "host_geometry.h"

/** @todo Remove these undefs, only to here to hide global defs until those are removed */
#undef NSPACE
#undef VOLSPACE
#undef NTIME

int get_neighbor(int nspace, int dir, const inputparameters * const params)
{
	int coord[NDIM];
	coord[0] = 0;
	for(int j = 1; j < NDIM; j++) coord[j] = get_spacecoord(nspace, j, params);
	coord[dir] = (coord[dir] + 1) % params->get_ns();
	return get_nspace(coord, params);
}
int get_lower_neighbor(int nspace, int dir, const inputparameters * const params)
{
	const size_t NSPACE = params->get_ns();
	int coord[NDIM];
	coord[0] = 0;
	for(int j = 1; j < NDIM; j++) coord[j] = get_spacecoord(nspace, j, params);
	coord[dir] = (coord[dir] - 1 + NSPACE) % NSPACE;
	return get_nspace(coord, params);
}
int get_nspace(int* coord, const inputparameters * const params)
{
	int n = 0;
	for(int j = 1; j < NDIM; j++) n += pow(params->get_ns(), j - 1) * coord[j];
	return n;
}
int get_spacecoord(int nspace, int dir, const inputparameters * const params)
{
	const size_t NSPACE = params->get_ns();
	int nred = NDIM - 1;
	int res = int(nspace / pow(NSPACE, nred - 1));
	if(dir == nred) return res;
	int acc = res;
	for(int j = 1; j < nred; j++) {
		res = int(nspace / pow(NSPACE, nred - 1 - j)) - NSPACE * acc;
		if(dir == nred - j) {
			return res;
		}
		acc = NSPACE * acc + res;
	}
	return -99;
}

/** @todo CP: the general structure here should be overviewed. It would be good to have the spinors structured like
 * (cv0, cv1, cv2, cv3) with cv being a colorvector and not like (dv0, dv1, dv2) with dv being a Diracvector!!
 * At least in the OpenCL-Code this should be done! */
int spinor_color(int spinor_element, const inputparameters * const)
{
	return (int)(spinor_element / NSPIN);
}

int spinor_spin(int spinor_element, int color, const inputparameters * const)
{
	return spinor_element - NSPIN * color;
}

int spinor_element(int alpha, int color)
{
	return alpha + NSPIN * color;
}

int eoprec_spinor_field_element(int alpha, int color, int nspace, int t, const inputparameters * const params)
{
	return alpha + NSPIN * color + NSPIN * NC * get_n_eoprec(t, nspace, params);
}

int eoprec_spinor_field_element(int alpha, int color, int n_eoprec, const inputparameters * const)
{
	return alpha + NSPIN * color + NSPIN * NC * n_eoprec;
}

//!!CP: changed the args to fit all other functions!!
int get_n_eoprec(int spacepos, int timepos, const inputparameters * const params)
{
	return (int)((get_global_pos(spacepos, timepos, params)) / 2);
}

int spinor_field_element(int alpha, int color, int nspace, int t, const inputparameters * const params)
{
	return alpha + NSPIN * color + NSPIN * NC * (get_global_pos(nspace, t, params));
}

//functions that have explicite spatial and temporal positions in them
//make it:
//site = pos + VOLSPACE*t =  x + y*NSPACE + z*NSPACE*NSPACE + VOLSPACE*t
//idx = mu + NDIM*site
int get_global_pos(int spacepos, int t, const inputparameters * const params)
{
	return spacepos + params->get_volspace() * t;
}

int get_global_link_pos(int mu, int spacepos, int t, const inputparameters * const params)
{
	return mu + NDIM * get_global_pos(spacepos, t, params);
}

//dispensable
/*
int ocl_gaugefield_element(int c, int a, int b, int mu, int spacepos, int t){
#ifdef _RECONSTRUCT_TWELVE_
  //old:
//  return c + 2*a + 2*(NC-1)*b+2*NC*(NC-1)*mu+2*NC*(NC-1)*NDIM*spacepos+2*NC*(NC-1)*NDIM*VOLSPACE*t;
  return c + 2*( a + (NC-1)*b + NC*(NC-1)* ( get_global_link_pos(mu, spacepos, t) ) );
  #else
  //old:
//  return c + 2*a + 2*NC*b+2*NC*NC*mu+2*NC*NC*NDIM*spacepos+2*NC*NC*NDIM*VOLSPACE*t;
   return c + 2* (a + NC*b+ NC*NC*( get_global_link_pos(mu, spacepos, t) ) );
#endif
}
*/

//it is assumed that idx iterates only over half the number of sites
void get_even_site(int idx, int * out_space, int * out_t, const inputparameters * const params)
{
	const size_t NSPACE = params->get_ns();
	const size_t VOLSPACE = params->get_volspace();
	int x, y, z, t;
	x = idx;
	t = (int)(idx / (VOLSPACE / 2));
	x -= t * VOLSPACE / 2;
	z = (int)(x / (NSPACE * NSPACE / 2));
	x -= z * NSPACE * NSPACE / 2;
	y = (int)(x / NSPACE);
	x -= y * NSPACE;
	(*out_space) =  (int)((z + t) % 2) * (1 + 2 * x - (int) (2 * x / NSPACE)) + (int)((t + z + 1) % 2) * (2 * x + (int) (2 * x / NSPACE)) + 2 * NSPACE * y + NSPACE * NSPACE * z;
	(*out_t) = t;
}

//it is assumed that idx iterates only over half the number of sites
void get_odd_site(int idx, int * out_space, int * out_t, const inputparameters * const params)
{
	const size_t NSPACE = params->get_ns();
	const size_t VOLSPACE = params->get_volspace();
	int x, y, z, t;
	x = idx;
	t = (int)(idx / (VOLSPACE / 2));
	x -= t * VOLSPACE / 2;
	z = (int)(x / (NSPACE * NSPACE / 2));
	x -= z * NSPACE * NSPACE / 2;
	y = (int)(x / NSPACE);
	x -= y * NSPACE;
	(*out_space) =  (int)((z + t + 1) % 2) * (1 + 2 * x - (int) (2 * x / NSPACE)) + (int)((t + z) % 2) * (2 * x + (int) (2 * x / NSPACE)) + 2 * NSPACE * y + NSPACE * NSPACE * z;
	(*out_t) = t;
}
