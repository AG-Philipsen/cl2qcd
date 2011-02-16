#ifndef _GEOMETRYH_
#define _GEOMETRYH_

#include "globaldefs.h"
#include "hmcerrs.h"
#include "cstdlib"
#include "cmath"

//coord[0] = t
//coord[1] = x
//coord[2] = y
//coord[3] = z

int get_n_eoprec(int timepos, int spacepos);

//switch between (x,y,z) <-> nspace=0,...,VOLSPACE-1
int get_nspace(int* coord);
int get_spacecoord(int nspace, int dir);

int get_neighbor(int nspace, int dir);
int get_lower_neighbor(int nspace, int dir);

extern int* nspace_from_even_index;
extern int* ntime_from_even_index;
extern int* nspace_from_odd_index;
extern int* ntime_from_odd_index;
int get_ntime_from_eoprecindex(int n, int which);
int get_nspace_from_eoprecindex(int n, int which);

//Checkerboard:
void get_even_site(int idx, int * out_space, int * out_t);
void get_odd_site(int idx, int * out_space, int * out_t);

int get_global_pos(int spacepos, int t);

int ocl_gaugefield_element(int c, int a, int b, int mu, int spacepos, int t);

//Spinor functions
int spinor_element(int alpha, int color);
int spinor_field_element(int alpha, int color, int nspace, int t);
int spinor_color(int spinor_element);
int spinor_spin(int spinor_element,int color);
int eoprec_spinor_field_element(int alpha, int color, int nspace, int t);
int eoprec_spinor_field_element(int alpha, int color, int n_eoprec);

#endif
