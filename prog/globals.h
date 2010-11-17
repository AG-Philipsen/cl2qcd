#ifndef _GLOBALSH_
#define _GLOBALSH_

#define NC 3
#define NSPIN 4
#define NDIM 4

//it could be a good idea to define Nt and Ns at compile time
//usually you stick to one volume for quite a while anyways...
#define NSPACE 32
#define NTIME 12
int const VOLSPACE  = NSPACE*NSPACE*NSPACE;
int const VOL4D = VOLSPACE*NTIME;

#endif
