#include "host_geometry.h"

int get_neighbor(int nspace, int dir) {
  int coord[NDIM];
  coord[0]=0;
  for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(nspace,j);
  coord[dir] = (coord[dir] + 1)%NSPACE;
  return get_nspace(coord);
}
int get_lower_neighbor(int nspace, int dir) {
  int coord[NDIM];
  coord[0]=0;
  for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(nspace,j);
  coord[dir] = (coord[dir] - 1 + NSPACE)%NSPACE;
  return get_nspace(coord);
}
int get_nspace(int* coord){
  int n=0;
  for(int j=1; j<NDIM; j++) n += pow(NSPACE,j-1)*coord[j];
  return n;
}
int get_spacecoord(int nspace, int dir){
  int nred=NDIM-1;
  int res = int(nspace/pow(NSPACE,nred-1));
  if(dir==nred) return res;
  int acc = res;
  for(int j=1; j<nred; j++) {
      res = int(nspace/pow(NSPACE,nred-1-j)) - NSPACE*acc;
      if(dir==nred-j) {
	return res;
      }
      acc = NSPACE*acc + res;
  }
  return -99;
}

int spinor_color(int spinor_element){
  return (int)(spinor_element/NSPIN);
}

int spinor_spin(int spinor_element,int color){
  return spinor_element - NSPIN*color;
}

int spinor_element(int alpha, int color) {
  return alpha + NSPIN*color;
}

int eoprec_spinor_field_element(int alpha, int color, int nspace, int t) {
  return alpha + NSPIN*color + NSPIN*NC*get_n_eoprec(t,nspace);
}

int eoprec_spinor_field_element(int alpha, int color, int n_eoprec) {
  return alpha + NSPIN*color + NSPIN*NC*n_eoprec;
}

//!!CP: changed the args to fit all other functions!!
int get_n_eoprec(int spacepos, int timepos){
  return (int)((get_global_pos(spacepos, timepos))/2);
}

int spinor_field_element(int alpha, int color, int nspace, int t) {
  return alpha + NSPIN*color + NSPIN*NC*(get_global_pos(nspace, t));
}

//functions that have explicite spatial and temporal positions in them
//make it:
//site = pos + VOLSPACE*t =  x + y*NSPACE + z*NSPACE*NSPACE + VOLSPACE*t
//idx = mu + NDIM*site
int get_global_pos(int spacepos, int t){
  return spacepos + VOLSPACE * t;
}

int get_global_link_pos(int mu, int spacepos, int t){
  return mu + NDIM*get_global_pos(spacepos, t);
}

int ocl_gaugefield_element(int c, int a, int b, int mu, int spacepos, int t){
#ifdef _RECONSTRUCT_TWELVE_
	//old: 
//	return c + 2*a + 2*(NC-1)*b+2*NC*(NC-1)*mu+2*NC*(NC-1)*NDIM*spacepos+2*NC*(NC-1)*NDIM*VOLSPACE*t;
	return c + 2*( a + (NC-1)*b + NC*(NC-1)* ( get_global_link_pos(mu, spacepos, t) ) );
	#else
	//old:
//	return c + 2*a + 2*NC*b+2*NC*NC*mu+2*NC*NC*NDIM*spacepos+2*NC*NC*NDIM*VOLSPACE*t;
	 return c + 2* (a + NC*b+ NC*NC*( get_global_link_pos(mu, spacepos, t) ) );
#endif
}

//it is assumed that idx iterates only over half the number of sites
void get_even_site(int idx, int * out_space, int * out_t){
  int x,y,z,t;
  x = idx;
  t = (int)(idx/(VOLSPACE/2));
  x -= t*VOLSPACE/2;
  z = (int)(x/(NSPACE*NSPACE/2));
  x -= z*NSPACE*NSPACE/2;
  y = (int)(x/NSPACE);
  x -= y*NSPACE;
  (*out_space) =  (int)((z+t)%2)*(1 + 2*x - (int) (2*x/NSPACE)) + (int)((t+z+1)%2)*(2*x + (int) (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
  (*out_t) = t;
}

//it is assumed that idx iterates only over half the number of sites
void get_odd_site(int idx, int * out_space, int * out_t){
  int x,y,z,t;
  x = idx;
  t = (int)(idx/(VOLSPACE/2));
  x -= t*VOLSPACE/2;
  z = (int)(x/(NSPACE*NSPACE/2));
  x -= z*NSPACE*NSPACE/2;
  y = (int)(x/NSPACE);
  x -= y*NSPACE;
  (*out_space) =  (int)((z+t+1)%2)*(1 + 2*x - (int) (2*x/NSPACE)) + (int)((t+z)%2)*(2*x + (int) (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
  (*out_t) = t;
}
