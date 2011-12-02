#include "host_geometry.h"

int get_nspace(int* coord, const inputparameters * const params)
{
	int n = 0;
	for(int j = 1; j < NDIM; j++) n += pow(params->get_ns(), j - 1) * coord[j];
	return n;
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
