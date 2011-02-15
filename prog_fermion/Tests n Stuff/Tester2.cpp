#include "globaldefs.h"

#include <stdlib.h>
#include <stdio.h>


int main (int argc, char ** argv){
int x, y, z, t;
int pos, pos2;
int cter = 0;
int check = 0;

for (int i = 0; i<VOL4D; i++) check += i;

for(int idx = 0; idx<VOL4D/2; idx++){
  x = idx;
  t = int(idx/(VOLSPACE/2));
  x -= t*VOLSPACE/2;
  z = int(x/(NSPACE*NSPACE/2));
  x -= z*NSPACE*NSPACE/2;
  y = int(x/NSPACE);
  x -= y*NSPACE;
  
 
  
  //red
  pos  = int((z+t)%2)*(1 + 2*x - int (2*x/NSPACE)) + int((z+t+1)%2)*(2*x + int (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
  //black
  pos2 = int((z+t+1)%2)*(1 + 2*x - int (2*x/NSPACE)) + int((z+t)%2)*(2*x + int (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
  cter += pos + VOLSPACE*t;
  cter += pos2 + VOLSPACE*t;
}
	  printf("%i\t%i\n", cter, check);
  return 0;
}