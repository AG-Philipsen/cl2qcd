/** @file
 * Device code for random number generation.
 */
//opencl_random.cl

/** Type for random number generator state */
typedef ulong4 hmc_ocl_ran;

/**
 * Draw a 64-bit random integer using the algorithm described in Numerical Recipes 3.
 *
 * @param[in,out] state Pointer to this threads random number generator state in global memory.
 * @return A pseudo-random integer
 */
inline ulong nr3_int64(__global hmc_ocl_ran * state ) {
	(*state).x = (*state).x * 2862933555777941757L + 7046029254386353087L;
	(*state).y ^= (*state).y >> 17; (*state).y ^= (*state).y << 31; (*state).y ^= (*state).y >> 8;
	(*state).z = 4294957665U*((*state).z & 0xffffffff) + ((*state).z >> 32);
	ulong tmp = (*state).x ^ ((*state).x << 21); tmp ^= tmp >> 35; tmp ^= tmp << 4;
	return (tmp + (*state).y) ^ (*state).z;
}
/**
 * Draw a 32-bit random float using the algorithm described in Numerical Recipes 3.
 *
 * @param[in,out] state Pointer to this threads random number generator state in global memory
 * @return A pseudo-random float
 */
inline float ocl_new_ran(__global hmc_ocl_ran * state ){
	return 5.42101086242752217E-20f * nr3_int64( state );
}
/**
 * Draw a 32-bit random integer using the algorithm described in Numerical Recipes 3.
 *
 * @param[in,out] state Pointer to this threads random number generator state in global memory
 * @return A pseudo-random integer
 */
inline uint nr3_int32(__global hmc_ocl_ran * state ){
	return (uint) nr3_int64( state );
}
/**
 * Draw a 32-bit random integer in the range [0,range) using the algorithm described in Numerical Recipes 3.
 *
 * @param[in] range Upper bound for the drawn number, nummber will be one less than this at maximum
 * @param[in,out] state Pointer to this threads random number generator state in global memory
 * @return A pseudo-random integer
 */
int random_int( int range, __global hmc_ocl_ran* state ){
	return (nr3_int64( state ) % range);
}
/**
 * Get 1,2,3 in random order
 *
 * @param[out] rand Storage location for the result
 * @param[in,out] rnd Pointer to this threads random number generator state in global memory
 */
void random_1_2_3 (int rand[3], __global hmc_ocl_ran * rnd) { 
  rand[ 0 ] = random_int( 3, rnd ) + 1;
  do{
    rand[ 1 ] = random_int( 3, rnd ) + 1;
  } while( rand[ 0 ] == rand[ 1 ] );
  rand[ 2 ] = 6 - rand[ 1 ] - rand[ 0 ];
}
/**
 * Construct new SU2 matrix using improved alg by Kennedy Pendleton
 *
 * @param[out] dst Storage location for the resulting SU2
 * @param[in,out] rnd Pointer to this threads random number generator state in global memory
 */
void SU2Update(__private hmc_float dst [su2_entries], const hmc_float alpha, __global hmc_ocl_ran * rnd){
  hmc_float delta;
  hmc_float a0 ;
  hmc_float eta ; 
  do {
    delta = -log(ocl_new_ran(rnd))/alpha*pow(cos(2. * PI * ocl_new_ran(rnd)), 2.) -log(ocl_new_ran(rnd))/alpha;
    a0 = 1.-delta;
    eta = ocl_new_ran(rnd);
  } while ( (1.-0.5*delta) < eta*eta); 
  hmc_float phi = 2.*PI*ocl_new_ran(rnd);
  hmc_float theta = asin(2.*ocl_new_ran(rnd) - 1.);
  dst[0] = a0;
  dst[1] = sqrt(1.-a0 * a0)*cos(theta) * cos(phi);
  dst[2] = sqrt(1.-a0 * a0)*cos(theta) * sin(phi);
  dst[3] = sqrt(1.-a0 * a0)*sin(theta);
}
