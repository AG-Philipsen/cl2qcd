#ifndef _GEOMETRYH_
#define _GEOMETRYH_

#include "globaldefs.h"
#include "cstdlib"
#include "cmath"
#include "hmcerrs.h"

//coord[0] = t
//coord[1] = x
//coord[2] = y
//coord[3] = z

//switch between (x,y,z) <-> nspace=0,...,VOLSPACE-1

hmc_error init_geometry();
hmc_error finalize_geometry();

int get_nspace(int* coord);
int get_spacecoord(int nspace, int dir);

int get_neighbor(int nspace, int dir);
int get_lower_neighbor(int nspace, int dir);

int get_n_eoprec(int timepos, int spacepos);

int get_nspace_from_eoprecindex(int n, int which);
int get_ntime_from_eoprecindex(int n, int which);

extern int* nspace_from_even_index;
extern int* ntime_from_even_index;
extern int* nspace_from_odd_index;
extern int* ntime_from_odd_index;

#endif
