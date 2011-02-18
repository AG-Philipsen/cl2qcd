//CP: old generator
/*
typedef uint4 hmc_ocl_ran;

inline unsigned _ocl_taus_step( unsigned *z, int S1, int S2, int S3, unsigned M)
{
	unsigned b = ( ( ( (*z) << S1 )^(*z) ) >> S2 );
	return *z = ( ( ( (*z) & M ) << S3 )^b );
}

inline unsigned _ocl_LCG_step( unsigned *z, unsigned A, unsigned C)
{
	return *z= ( A * (*z) + C );
}

inline float ocl_new_ran( __global hmc_ocl_ran* state )
{
	uint x,y,z,w;
	x = (*state).x;
	y = (*state).y;
	z = (*state).z;
	w = (*state).w;
	float res = 2.328306436538696e-10f*( _ocl_taus_step( &x, 13, 19, 12, 4294967294ul)^
	                                     _ocl_taus_step( &y, 2, 25, 4,   4294967288ul)^
	                                     _ocl_taus_step( &z, 3, 11, 17,  4294967280ul)^
	                                     _ocl_LCG_step(  &w, 1664525,    1013904223ul) );
	*state = (hmc_ocl_ran)(x,y,z,w);
	return res;
}
*/

//CP: new one

/**
 * RNG state for the NR3
 */
// typedef ulong4 nr3_state;
typedef ulong4 hmc_ocl_ran;


/**
 * Calculate the next random number as described in NR3
 */
// inline ulong nr3_int64( nr3_state * state ) {
inline ulong nr3_int64(__global hmc_ocl_ran * state ) {
	(*state).x = (*state).x * 2862933555777941757L + 7046029254386353087L;
	(*state).y ^= (*state).y >> 17; (*state).y ^= (*state).y << 31; (*state).y ^= (*state).y >> 8;
	(*state).z = 4294957665U*((*state).z & 0xffffffff) + ((*state).z >> 32);
	ulong tmp = (*state).x ^ ((*state).x << 21); tmp ^= tmp >> 35; tmp ^= tmp << 4;
	return (tmp + (*state).y) ^ (*state).z;
}

/**
 * Calculate the next random number and return as float
 * FIXME this conversion is probably broken
 */
// inline float nr3_float( nr3_state * state )
inline float ocl_new_ran(__global hmc_ocl_ran * state )
{
	return 5.42101086242752217E-20f * nr3_int64( state );
}

/**
 * Calculate the next random number and return as int32
 */
// inline uint nr3_int32( nr3_state * state )
inline uint nr3_int32(__global hmc_ocl_ran * state )
{
	return (uint) nr3_int64( state );
}




int random_int( int range, __global hmc_ocl_ran* state )
{
	//return convert_int( ocl_new_ran( state ) * range );
	return (nr3_int64( state ) % range);
}

//returns 1,2,3 in a random way
void random_1_2_3 (int rand[3], __global hmc_ocl_ran * rnd) { 

	// first value can be any value 1..3
	rand[ 0 ] = random_int( 3, rnd ) + 1;

	// second value must be 1..3 but not equal to the first value
	do
	{
	  rand[ 1 ] = random_int( 3, rnd ) + 1;
	} while( rand[ 0 ] == rand[ 1 ] );

	// third value must take the remaining value from 1..3
	rand[ 2 ] = 6 - rand[ 1 ] - rand[ 0 ];
}

hmc_float calc_delta (__global hmc_ocl_ran* rnd , hmc_float alpha){
  return  (-log(ocl_new_ran(rnd))/alpha*pow(cos(2. * PI * ocl_new_ran(rnd)), 2.) -log(ocl_new_ran(rnd))/alpha );
}

hmc_float calc_eta (__global hmc_ocl_ran* rnd ){
  return  (ocl_new_ran(rnd) );
}

// Construct new SU2 matrix using improved alg by Kennedy Pendleton
void SU2Update(__private hmc_float dst [su2_entries], const hmc_float alpha, __global hmc_ocl_ran * rnd)
{
  hmc_float delta;
  hmc_float a0 ;
  hmc_float eta ; 
  do
  {
    delta = -log(ocl_new_ran(rnd))/alpha*pow(cos(2. * PI * ocl_new_ran(rnd)), 2.) -log(ocl_new_ran(rnd))/alpha;
//     delta = calc_delta(rnd, alpha);
    a0 = 1.-delta;
    eta = ocl_new_ran(rnd);
//     eta = calc_eta(rnd);
  }while ( (1.-0.5*delta) < eta*eta); 
  hmc_float phi = 2.*PI*ocl_new_ran(rnd);
  hmc_float theta = asin(2.*ocl_new_ran(rnd) - 1.);
  dst[0] = a0;
  dst[1] = sqrt(1.-a0 * a0)*cos(theta) * cos(phi);
  dst[2] = sqrt(1.-a0 * a0)*cos(theta) * sin(phi);
  dst[3] = sqrt(1.-a0 * a0)*sin(theta);
}

