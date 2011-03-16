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

int get_n_eoprec(int spacepos,int timepos);

//switch between (x,y,z) <-> nspace=0,...,VOLSPACE-1
int get_nspace(int* coord);
int get_spacecoord(int nspace, int dir);
//get spatial neighbors
int get_neighbor(int nspace, int dir);
int get_lower_neighbor(int nspace, int dir);
//Checkerboard: get real coordinates from EVEN/ODD-index
void get_even_site(int idx, int * out_space, int * out_t);
void get_odd_site(int idx, int * out_space, int * out_t);

int get_global_pos(int spacepos, int t);
//get gaugefield element from long array
int ocl_gaugefield_element(int c, int a, int b, int mu, int spacepos, int t);

//Spinor functions
int spinor_element(int alpha, int color);
int spinor_field_element(int alpha, int color, int nspace, int t);
int spinor_color(int spinor_element);
int spinor_spin(int spinor_element,int color);
int eoprec_spinor_field_element(int alpha, int color, int nspace, int t);
int eoprec_spinor_field_element(int alpha, int color, int n_eoprec);

#endif
