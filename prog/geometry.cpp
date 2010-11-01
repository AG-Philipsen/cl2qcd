#include "geometry.h"

int get_neighbor(int nspace, int dir) {
  int coord[NDIM];
  coord[0]=0;
  for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(nspace,j);
  coord[dir] = (coord[dir] + 1)%NSPACE;
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



