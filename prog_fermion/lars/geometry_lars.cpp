#include "geometry.h"

hmc_error init_geometry(){
  nspace_from_even_index = new int[VOL4D/2];
  ntime_from_even_index = new int[VOL4D/2];
  nspace_from_odd_index = new int[VOL4D/2];
  ntime_from_odd_index = new int[VOL4D/2];
  for(int x=0; x<NSPACE; x++) {
    for(int y=0; y<NSPACE; y++) {
      for(int z=0; z<NSPACE; z++) {
	for(int t=0; t<NTIME; t++) {
	  int coord[NDIM] = {t,x,y,z};
	  int nspace = get_nspace(coord);
	  int n = (int)( (t + NTIME*nspace)/2 );
	  int check = t+x+y+z;
	  if(check%2==0) {
	    nspace_from_even_index[n] = nspace;
	    ntime_from_even_index[n]  = t;
	  } else {
	    nspace_from_odd_index[n] = nspace;
	    ntime_from_odd_index[n]  = t;
	  }
	}
      }
    }
  }
  return HMC_SUCCESS;
}

hmc_error finalize_geometry(){
  delete [] nspace_from_even_index;
  delete [] ntime_from_even_index;
  delete [] nspace_from_odd_index;
  delete [] ntime_from_odd_index;
  return HMC_SUCCESS;
}


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

int get_n_eoprec(int timepos, int spacepos){
  return (int)((timepos+NTIME*spacepos)/2);
}

int get_nspace_from_eoprecindex(int n, int which){
  if(which==EVEN) return nspace_from_even_index[n];
  if(which==ODD) return nspace_from_odd_index[n];
  exit(HMC_EOERROR);
}

int get_ntime_from_eoprecindex(int n, int which){
  if(which==EVEN) return ntime_from_even_index[n];
  if(which==ODD) return ntime_from_odd_index[n];
  exit(HMC_EOERROR);
}
