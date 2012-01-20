#include "host_geometry.h"

int get_nspace(int* coord, const inputparameters * const params)
{
	int n = 0;
	n = params->get_ns() * params->get_ns() * coord[3] + params->get_ns() * coord[2] + coord[1];
	//for(int j = 1; j < NDIM; j++) n += pow(params->get_ns(), j - 1) * coord[j];
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

size_t get_su3_idx_ildg_format(size_t n, size_t m, size_t x, size_t y, size_t z, size_t t, size_t mu, const inputparameters * const parameters){
        //ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
        size_t link_idx_ildg = get_link_idx_ildg_format(x,y,z,t,mu, parameters);        
	//the 2 is introduced because we are dealing with complex numbers here..
	return 2 * n + 2 * m * NC + link_idx_ildg * NC * NC * 2;
}

size_t get_link_idx_ildg_format(size_t x, size_t y, size_t z, size_t t, size_t mu, const inputparameters * const parameters){
	const size_t VOLSPACE = parameters->get_volspace();
	const size_t NSPACE = parameters->get_ns();
	size_t spacepos_ildg = z + y * NSPACE + x * NSPACE * NSPACE;
	size_t link_idx_ildg = mu + spacepos_ildg * NDIM + t * VOLSPACE * NDIM;
        return link_idx_ildg;
}
