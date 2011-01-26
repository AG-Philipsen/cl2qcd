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


/*
int get_nspace(int* coord){
  int n = coord[1] +  NSPACE*coord[2] + NSPACE*NSPACE*coord[3];
  return n;
}

int get_spacecoord(int nspace, int dir){
  int res = int(nspace/(NSPACE*NSPACE));
  if(dir==3) return res;
  int acc = res;
  res = int(nspace/(NSPACE)) - NSPACE*acc;
  if(dir==2) return res;
  acc = NSPACE*acc + res;
  res = nspace - NSPACE*acc;
  return res;
}

void get_allspacecoord(int nspace, int coord[NDIM]){
  int res = int(nspace/(NSPACE*NSPACE));
  coord[3] = res;
  int acc = res;
  res = int(nspace/(NSPACE)) - NSPACE*acc;
  coord[2] = res;
  acc = NSPACE*acc + res;
  res = nspace - NSPACE*acc;
  coord[1] = res;
}

int get_neighbor(int nspace,int dir) {
  int coord[NDIM];
//   coord[0]=0;
//   for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(nspace,j);
  get_allspacecoord(nspace, coord);
  coord[dir] = (coord[dir] + 1)%NSPACE;
  return get_nspace(coord);
}

int get_lower_neighbor(int nspace, int dir) {
  int coord[NDIM];
//   coord[0]=0;
//   for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(nspace,j);
  get_allspacecoord(nspace, coord);
  coord[dir] = (coord[dir] - 1 + NSPACE)%NSPACE;
  return get_nspace(coord);
}
*/
