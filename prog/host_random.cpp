#include "host_random.h"

//this should be avoidable by using the extern Random rnd
// Random zufall(seed);

//Zufallszahl 1,2,3 vom Typ int 
inline int random_123 ()
{
  //return zufall.int64() % 3 +1;
  return rnd.int64() % 3 +1;
}

//Gibt drei Zufallszahlen 1,2,3
void random_1_2_3 (int rand[3])
{
  rand[0] = random_123();
  do
    {rand[1] = random_123();}
  while (rand[1] == rand[0]);
  rand[2] = 6 - rand[1] - rand[0];
}

//old
/*
void init_random_seeds(Random random, hmc_ocl_ran * out, const int NUM, usetimer * timer){
  (*timer).reset();
  int dummy;
  //The initial values for the first the x, y and z components of the state should be > 128
  for(int i = 0; i<NUM; i++){
    do{
      dummy = random.int64();
    } while(dummy < 128);
    (out[i]).x = dummy;
    do{
      dummy = random.int64();
    } while(dummy < 128);
    (out[i]).y = dummy;
    do{
      dummy = random.int64();
    } while(dummy < 128);
    (out[i]).z = dummy;    
    dummy = random.int64();
    (out[i]).w = dummy;    
  }
  (*timer).add();
  return;
}
*/

#define CLU_VEC( vec, idx ) (vec).s[idx]

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

//CP: Original by MB
/*
inline void nr3_init_state( nr3_state * state, cl_ulong seed )
{
	CLU_VEC(*state,1) = 4101842887655102017L; CLU_VEC(*state,2) = 1;
	// TODO abort if seed > y
	CLU_VEC(*state,0) = seed ^ CLU_VEC(*state,1); nr3_int64( state );
	CLU_VEC(*state,1) = CLU_VEC(*state,0); nr3_int64( state );
	CLU_VEC(*state,2) = CLU_VEC(*state,1); nr3_int64( state );
}
*/

void init_random_seeds(Random random, hmc_ocl_ran * out, const int NUM, usetimer * timer){
  (*timer).reset();
  int dummy;
  hmc_ocl_ran initializer;
  nr3_init_state( &initializer, 50000 );
  for(int i = 0; i<NUM; i++){
    while( (dummy = nr3_int64(&initializer) ) >= 4101842887655102017L ) { };
		nr3_init_state( &out[i], nr3_int64(&initializer) );
  }
   
  (*timer).add();
  return;
}

//CP: new random code by MB
/*

inline void nr3_init_states( size_t num_states, nr3_state states[], cl_ulong seed )
{
	// create generator to genereate random seeds
	nr3_state initializer;
	nr3_init_state( &initializer, seed );

	// initialize using seeds from a single generator
	for( size_t i = 0; i < num_states; ++i )
	{
		while( (seed = nr3_int64(&initializer) ) >= 4101842887655102017L ) { };
		nr3_init_state( &states[i], nr3_int64(&initializer) );
	}
}

*/
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

