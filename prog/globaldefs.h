#ifndef _GLOBALSH_
#define _GLOBALSH_

#define NC 3
#define NSPIN 4
#define NDIM 4

#ifndef _INKENREL_
//this should not be read by kernel at runtime
//NTIME and NSPACE must be fixed at compile time

//it could be a good idea to define Nt and Ns at compile time
//usually you stick to one volume for quite a while anyways...
#define NSPACE 4
#define NTIME 4
#define VOLSPACE NSPACE*NSPACE*NSPACE
#define VOL4D VOLSPACE*NTIME

#define su2_entries 4

#define PI 	3.14159265358979

//startconditions:
#define START_FROM_SOURCE 2
#define COLD_START 0
#define HOT_START 1

#endif //_INKERNEL_

#endif
