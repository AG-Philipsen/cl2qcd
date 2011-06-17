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
inline ulong nr3_int64(__global hmc_ocl_ran * state )
{
	(*state).x = (*state).x * 2862933555777941757L + 7046029254386353087L;
	(*state).y ^= (*state).y >> 17;
	(*state).y ^= (*state).y << 31;
	(*state).y ^= (*state).y >> 8;
	(*state).z = 4294957665U*((*state).z & 0xffffffff) + ((*state).z >> 32);
	ulong tmp = (*state).x ^ ((*state).x << 21);
	tmp ^= tmp >> 35;
	tmp ^= tmp << 4;
	return (tmp + (*state).y) ^ (*state).z;
}
/**
 * Draw a 32-bit random float using the algorithm described in Numerical Recipes 3.
 *
 * @param[in,out] state Pointer to this threads random number generator state in global memory
 * @return A pseudo-random float
 */
inline float ocl_new_ran(__global hmc_ocl_ran * state )
{
	return 5.42101086242752217E-20f * nr3_int64( state );
}
/**
 * Draw a 32-bit random integer using the algorithm described in Numerical Recipes 3.
 *
 * @param[in,out] state Pointer to this threads random number generator state in global memory
 * @return A pseudo-random integer
 */
inline uint nr3_int32(__global hmc_ocl_ran * state )
{
	return (uint) nr3_int64( state );
}
/**
 * Draw a 32-bit random integer in the range [0,range) using the algorithm described in Numerical Recipes 3.
 *
 * @param[in] range Upper bound for the drawn number, nummber will be one less than this at maximum
 * @param[in,out] state Pointer to this threads random number generator state in global memory
 * @return A pseudo-random integer
 */
int random_int( int range, __global hmc_ocl_ran* state )
{
	return (nr3_int64( state ) % range);
}
/**
 * Get 1,2,3 in random order
 *
 * @param[in,out] rnd Pointer to this threads random number generator state in global memory
 * @return A random permutation of 1,2,3 in x,y,z.
 */
void random_1_2_3(__private int seq[3], __global hmc_ocl_ran * const restrict state)
{
	/// @todo using uint3 as a return type instead of using an arg for return value
	/// would reduce register usage on cypress by two GPR

	// calculate everythin in (0..2) and add 1 in the end
	//
	// 1. is drawn from 0..2
	// 2. is drawn from 1..2 and processed as follows
	//  - 1. == 0
	//  -- 0 + 1 = 1
	//  -- 0 + 2 = 2
	//  - 1. == 1
	//  -- 1 + 1 = 2
	//  -- 1 + 2 = 3 -> 0
	//  - 1. == 2
	//  -- 2 + 1 = 3 -> 0
	//  -- 2 + 2 = 4 -> 0
	//    note how this keeps the probabilities correct by doing only a 50% roll and then a bijective mapping.
	// 3. as above the remaing number by difference to the fixed sum of 0+1+2=3.

// using the commented lines instead of the actually used ones
// will increase register usage by 2 GPRs on the Cypress, so it
// is essentially the same code
//	seq[0] = nr3_int64(state) % 3;
//	seq[1] = (seq[0] + (nr3_int64(state) % 2) + 1) % 3;
	seq[0] = random_int(3, state);
	seq[1] = (seq[0] + random_int(2, state) + 1) % 3;
	seq[2] = 3 - seq[1] - seq[0];

	++seq[0];
	++seq[1];
	++seq[2];
}