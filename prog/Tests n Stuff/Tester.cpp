#include "globaldefs.h"

#include <stdlib.h>
#include <stdio.h>


int main (int argc, char ** argv){
  int idx;
//   idx = atoi(argv[1]);
  printf("%i %i %i %i\n", NSPACE, NTIME, VOLSPACE, VOLSPACE/2);
  int x, y, z, t;
  int m, n, l, o;
//   x = idx; printf("%i \n", x);
//   t = idx - idx%(VOLSPACE/2);printf("%i \n", t);
//   x -= t*VOLSPACE/2;printf("%i \n", x);
//   z = x%(NSPACE*NSPACE/2);printf("%i \n", z);
//   x -= z*(NSPACE*NSPACE/2);printf("%i \n", x);
//   y = x%(NSPACE);printf("%i \n", y);
//   x -= y*NSPACE;printf("%i \n", x);
//   printf("%i %i %i %i \n", x, y, z, t);


    for (t = 0; t<NTIME; t++){
      for (z = 0; z< NSPACE; z++){ 
        for (y = 0; y<NSPACE/2; y++){
	  for (x = 0; x<NSPACE; x++){
          
	    idx = x + NSPACE*y + NSPACE*NSPACE/2*z + VOLSPACE/2*t;
	    m= idx;
	    o = int(idx/(VOLSPACE/2));
	    m -= o*VOLSPACE/2;
	    l = int(m/(NSPACE*NSPACE/2));
	    m -= l*(NSPACE*NSPACE/2);
	    n = int(m/(NSPACE));
	    m -= n*NSPACE;
// int pos = int((x+t)%2)*(1 + 2*z - int (2*z/NSPACE)) + int((x+1+t)%2)*(2*z + int (2*z/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*x;
int pos = int((z+t)%2)*(1 + 2*x - int (2*x/NSPACE)) + int((z+1+t)%2)*(2*x + int (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
  printf("%i\t%i %i %i %i \t", pos, m, n, l, o);
  printf("%i\t%i %i %i %i \n", idx, x, y, z, t);
	  }}}}
  return 0;
}