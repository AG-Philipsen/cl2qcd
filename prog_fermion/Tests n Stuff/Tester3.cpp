#include "globaldefs.h"

#include <stdlib.h>
#include <stdio.h>

//test even-odd preconditioning
int inline ocl_gaugefield_element(int c, int a, int b, int mu, int spacepos, int t){
#ifdef _RECONSTRUCT_TWELVE_
  return c + 2*a + 2*(NC-1)*b+2*NC*(NC-1)*mu+2*NC*(NC-1)*NDIM*spacepos+2*NC*(NC-1)*NDIM*VOLSPACE*t;
#else
  return c + 2*a + 2*NC*b+2*NC*NC*mu+2*NC*(NC-1)*NDIM*spacepos+2*NC*(NC-1)*NDIM*VOLSPACE*t;
#endif
}

void inline get_even_site(int idx, int * out_space, int * out_t){
  //one starts with one idx and has to get the four indices back
  int x,y,z,t;
  x = idx;
  t = int(idx/(VOLSPACE/2));
  x -= t*VOLSPACE/2;
  z = int(x/(NSPACE*NSPACE/2));
  x -= z*NSPACE*NSPACE/2;
  y = int(x/NSPACE);
  x -= y*NSPACE;
  (*out_space) =  int((z+t)%2)*(1 + 2*x - int (2*x/NSPACE)) + int((t+z+1)%2)*(2*x + int (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
  (*out_t) = t;
}

void inline get_odd_site(int idx, int * out_space, int * out_t){
  //one starts with one idx and has to get the four indices back
  int x,y,z,t;
  x = idx;
  t = int(idx/(VOLSPACE/2));
  x -= t*VOLSPACE/2;
  z = int(x/(NSPACE*NSPACE/2));
  x -= z*NSPACE*NSPACE/2;
  y = int(x/NSPACE);
  x -= y*NSPACE;
  (*out_space) =  int((z+t+1)%2)*(1 + 2*x - int (2*x/NSPACE)) + int((t+z)%2)*(2*x + int (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
  (*out_t) = t;
}

int inline get_global_pos(int spacepos, int t){
  return spacepos + VOLSPACE * t;
}

int main (int argc, char ** argv){
int x, y, z, t;
int pos, pos2;
int globalpos, globalpos2;
int cter = 0;
int check = 0;

for (int i = 0; i<VOL4D; i++) check += i;
printf("x y z t\t\teven odd\t\t global_even global_odd\n");
for(int idx = 0; idx<VOL4D/2; idx++){
  //one starts with one idx and has to get the four indices back
  //new: t, x, y, z
  //idx = t + NTIME*pos = t + NTIME( x + y*NSPACE/2 + z*NSPACE*NSPACE/2 )
//   t = idx;
//   z = int(idx/(NTIME*NSPACE*NSPACE/2));
//   t -= z*NTIME*NSPACE*NSPACE/2;
//   x = int(t/(NTIME*NSPACE/2));
//   t -= x*NSPACE*NTIME/2;
//   y = int(t/NTIME);
//   t -= y*NTIME;
  
  //old: x,y,z,t
  //idx = pos + VOLSPACE/2*t =  x + y*NSPACE/2 + z*NSPACE*NSPACE/2 + VOLSPACE/2*t
//   x = idx;
//   t = int(idx/(VOLSPACE/2));
//   x -= t*VOLSPACE/2;
//   z = int(x/(NSPACE*NSPACE/2));
//   x -= z*NSPACE*NSPACE/2;
//   y = int(x/NSPACE);
//   x -= y*NSPACE;
  
  //old: 
  //red
//   pos  = int((z+t)%2)*(1 + 2*x - int (2*x/NSPACE)) + int((t+z+1)%2)*(2*x + int (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
  //black
//   pos2 = int((z+t+1)%2)*(1 + 2*x - int (2*x/NSPACE)) + int((z+t)%2)*(2*x + int (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
  
  //new
//   //red
//   int newt;
//   newt = int((z+y)%2)*(1 + 2*t - int (2*t/NTIME)) + int((y+z+1)%2)*(2*t + int (2*t/NTIME));
//   pos  =  2*x + NSPACE*y + NSPACE*NSPACE*z;
//   pos2 = pos;
//   //black
//   int newt2; 
//   newt2 = int((z+y+1)%2)*(1 + 2*t - int (2*t/NTIME)) + int((z+y)%2)*(2*t + int (2*t/NTIME)); 
  
  //new
//     globalpos = t + NTIME*pos;
//     globalpos2 = t + NTIME*pos2 ;
  
  //old
//   globalpos = pos + VOLSPACE * t;
//   globalpos2 = pos2 + VOLSPACE * t;
  
  //part with inline functions:
  
  get_even_site(idx, &pos, &t);
  globalpos = get_global_pos(pos, t);
  
  get_odd_site(idx, &pos2, &t);
  globalpos2 = get_global_pos(pos2,t);
  
  
  printf("%i %i %i %i\t\t%i %i\t\t %i %i\n",x, y, z, t,  pos, pos2, globalpos, globalpos2);

  cter += globalpos;
  cter += globalpos2;
}
	  printf("%i\t%i\n", cter, check);

  return 0;
}

//formula-explanation:
  //(for the parameter order x,y,z,t)
  //this formula basically iterates through the x-y-plane and picks out the even/odd sites in this plane
  //this is done in dependence of the other two directions, which just move the starting point
  //(it is assumed that one has even directions)
  //example: 4^4
  //   in the x-y plane it looks like this:
  //   ...
  //   4 5 6 7
  //   0 1 2 3 
  //   so one has to pick the sites (0,2,5,7) and (1,3,4,6), respectivly, covering two rows in the y-direction. 
  //   After that, the numbers are repeated with an additional offset (here 8*, or 2*NSPACE in general). 
  //   if one now goes to 3 dimensions (z=1), the same is repeated, but one has to pick the even sites as (1,3,4,6) 
  //   and the odd sites like (0,2,5,7) because one would otherwise get the direct neighbours of the prior picked sites,
  //   so its just interchanged.
  //   In 4 dimensions(z and t), this means that one has to pick the first series of sites at an even z + t and the 
  //   second series at an odd z+t for the even sites and vice versa for odd sites.
  //
  //the term (2*x + int (2*x/NSPACE)      iterates through the x-direction and selects the even sites (0,2,..2*NSPACE-1)
  //the term (1 + 2*x - int (2*x/NSPACE)) iterates through the x-direction and selects the odd sites (1,3,..2*NSPACE)
  //the term 2*NSPACE*y                   iterates throght the y direction, leaving out one row (0,2,4..)
  //  (this covers the whole y direction since it is only iterates to NSPACE/2)
  //the term NSPACE*NSPACE*z              makes the spatial position complete
  //the term int((z+t)%2)                 switches between which series of numbers has to be taken
  //the term int((z+t+1)%2)                 switches between which series of numbers has to be taken
  