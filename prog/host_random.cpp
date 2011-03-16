#include "host_random.h"

inline int random_123 ()
{
  return rnd.int64() % 3 +1;
}

void random_1_2_3 (int rand[3])
{
  rand[0] = random_123();
  do
    {rand[1] = random_123();}
  while (rand[1] == rand[0]);
  rand[2] = 6 - rand[1] - rand[0];
}

#define CLU_VEC( vec, idx ) (vec).s[idx]

//CP: NR3-PRNG converted for OpenCL by MB
inline cl_ulong nr3_int64( hmc_ocl_ran * state )
{
	CLU_VEC(*state,0) = CLU_VEC(*state,0) * 2862933555777941757L + 7046029254386353087L;
	CLU_VEC(*state,1) ^= CLU_VEC(*state,1) >> 17; CLU_VEC(*state,1) ^= CLU_VEC(*state,1) << 31; CLU_VEC(*state,1) ^= CLU_VEC(*state,1) >> 8;
	CLU_VEC(*state,2) = 4294957665U*(CLU_VEC(*state,2) & 0xffffffff) + (CLU_VEC(*state,2) >> 32);
	cl_ulong tmp = CLU_VEC(*state,0) ^ (CLU_VEC(*state,0) << 21); tmp ^= tmp >> 35; tmp ^= tmp << 4;
	return (tmp + CLU_VEC(*state,1)) ^ CLU_VEC(*state,2);
}

inline void nr3_init_state( hmc_ocl_ran * state, cl_ulong seed )
{
	CLU_VEC(*state,1) = 4101842887655102017L; 
	CLU_VEC(*state,2) = 1;
	// TODO abort if seed > y
	CLU_VEC(*state,0) = seed ^ CLU_VEC(*state,1); nr3_int64( state );
	CLU_VEC(*state,1) = CLU_VEC(*state,0); nr3_int64( state );
	CLU_VEC(*state,2) = CLU_VEC(*state,1); nr3_int64( state );
}

void init_random_seeds(Random random, hmc_ocl_ran * out, usetimer * timer){
  (*timer).reset();
  int dummy;
  hmc_ocl_ran initializer;
  nr3_init_state( &initializer, 50000 );
  //it is assumed that the array is always of size NUMTHREADS
  for(int i = 0; i<NUMTHREADS; i++){
    while( (dummy = nr3_int64(&initializer) ) >= 4101842887655102017L ) { };
		nr3_init_state( &out[i], dummy );
  }
   
  (*timer).add();
  return;
}

void SU2Update(hmc_float dst [su2_entries], const hmc_float alpha)
{
  hmc_float delta;
  hmc_float a0 ;
  hmc_float eta ;
  do
  {
    delta = -log(rnd.doub())/alpha*pow(cos(2. * PI * rnd.doub()), 2.) -log(rnd.doub())/alpha;
    a0 = 1.-delta;
    eta = rnd.doub();
  }while ( (1.-0.5*delta) < eta*eta); 
  hmc_float phi = 2.*PI*rnd.doub();
  hmc_float theta = asin(2.*rnd.doub() - 1.);
  dst[0] = a0;
  dst[1] = sqrt(1.-a0 * a0)*cos(theta) * cos(phi);
  dst[2] = sqrt(1.-a0 * a0)*cos(theta) * sin(phi);
  dst[3] = sqrt(1.-a0 * a0)*sin(theta);
}

