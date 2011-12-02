#include "host_geometry.h"

int get_neighbor(int nspace, int dir, const inputparameters * const params)
{
	int coord[NDIM];
	coord[0] = 0;
	for(int j = 1; j < NDIM; j++) coord[j] = get_spacecoord(nspace, j, params);
	coord[dir] = (coord[dir] + 1) % params->get_ns();
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
