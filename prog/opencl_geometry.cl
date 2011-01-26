
int get_nspace(const int* coord){
  int n = coord[1] +  NSPACE*coord[2] + NSPACE*NSPACE*coord[3];
  return n;
}

int get_spacecoord(const int nspace, const int dir){
  int res = convert_int(nspace/(NSPACE*NSPACE));
  if(dir==3) return res;
  int acc = res;
  res = convert_int(nspace/(NSPACE)) - NSPACE*acc;
  if(dir==2) return res;
  acc = NSPACE*acc + res;
  res = nspace - NSPACE*acc;
  return res;
}

void get_allspacecoord(const int nspace, int coord[NDIM]){
  int res = convert_int(nspace/(NSPACE*NSPACE));
  coord[3] = res;
  int acc = res;
  res = convert_int(nspace/(NSPACE)) - NSPACE*acc;
  coord[2] = res;
  acc = NSPACE*acc + res;
  res = nspace - NSPACE*acc;
  coord[1] = res;
}

int get_neighbor(const int nspace,const int dir) {
  int coord[NDIM];
//   coord[0]=0;
//   for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(nspace,j);
  get_allspacecoord(nspace, coord);
  coord[dir] = (coord[dir] + 1)%NSPACE;
  return get_nspace(coord);
}

int get_lower_neighbor(const int nspace, int const dir) {
  int coord[NDIM];
//   coord[0]=0;
//   for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(nspace,j);
  get_allspacecoord(nspace, coord);
  coord[dir] = (coord[dir] - 1 + NSPACE)%NSPACE;
  return get_nspace(coord);
}
