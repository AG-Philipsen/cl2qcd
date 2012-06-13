/** @file
 * Device code for random number generation.
 */
//opencl_random.cl

#ifdef USE_PRNG_NR3
/** Type for random number generator state */
typedef ulong4 prng_state;
typedef prng_state rngStateStorageType;
#elif defined(USE_PRNG_RANLUX)
typedef ranluxcl_state_t prng_state;
typedef float4 rngStateStorageType;
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX

#ifdef USE_PRNG_NR3
/**
 * Draw a 64-bit random integer using the algorithm described in Numerical Recipes 3.
 *
 * @param[in,out] state Pointer to this threads random number generator state in global memory.
 * @return A pseudo-random integer
 */
inline ulong nr3_int64(prng_state * state )
{
	(*state).x = (*state).x * 2862933555777941757L + 7046029254386353087L;
	(*state).y ^= (*state).y >> 17;
	(*state).y ^= (*state).y << 31;
	(*state).y ^= (*state).y >> 8;
	(*state).z = 4294957665U * ((*state).z & 0xffffffff) + ((*state).z >> 32);
	ulong tmp = (*state).x ^ ((*state).x << 21);
	tmp ^= tmp >> 35;
	tmp ^= tmp << 4;
	return (tmp + (*state).y) ^ (*state).z;
}
#ifdef _USEDOUBLEPREC_
/**
 * Draw a 64-bit random float using the algorithm described in Numerical Recipes 3.
 *
 * @param[in,out] state Pointer to this threads random number generator state in global memory
 * @return A pseudo-random float
 */
inline double nr3_double(prng_state * state)
{
	return 5.42101086242752217E-20 * nr3_int64( state );
}
#endif /* _USEDOUBLEPREC_ */
/**
 * Draw a 42-bit random float using the algorithm described in Numerical Recipes 3.
 *
 * @param[in,out] state Pointer to this threads random number generator state in global memory
 * @return A pseudo-random float
 */
inline float nr3_float(prng_state * state)
{
	return 5.328306437e-10 * nr3_int32( state );
}
/**
 * Draw a 32-bit random integer using the algorithm described in Numerical Recipes 3.
 *
 * @param[in,out] state Pointer to this threads random number generator state in global memory
 * @return A pseudo-random integer
 */
inline uint nr3_int32(prng_state * state )
{
	return (uint) nr3_int64( state );
}
#endif // USE_PRNG_NR3

/**
 * Read the random number generate state from global mamory
 */
void prng_loadState(prng_state * const restrict state, __global rngStateStorageType * const restrict states)
{
#ifdef USE_PRNG_NR3
	*state = states[get_global_id(0)];
#elif defined(USE_PRNG_RANLUX)
	ranluxcl_download_seed(state, states);
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}

/**
 * Read the random number generate state from global mamory
 */
void prng_storeState(__global rngStateStorageType * const restrict states, prng_state * const restrict state)
{
#ifdef USE_PRNG_NR3
	states[get_global_id(0)] = *state;
#elif defined(USE_PRNG_RANLUX)
	ranluxcl_upload_seed(state, states);
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}

/**
 * Draw a 32-bit random integer in the range [0,range).
 *
 * @param[in] range Upper bound for the drawn number, nummber will be one less than this at maximum
 * @param[in,out] state Pointer to this threads random number generator state in global memory
 * @return A pseudo-random integer
 */
uint prng_int32(uint range, prng_state * const restrict state)
{
#ifdef USE_PRNG_NR3
	return nr3_int64(state) % range;
#elif defined(USE_PRNG_RANLUX)
	float4 tmp = ranluxcl(state);
	return tmp.x * range;
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}

#ifdef _USEDOUBLEPREC_
/**
 * Draw a double precision floating point number.
 *
 * @param[in] range Upper bound for the drawn number, nummber will be one less than this at maximum
 * @param[in,out] state Pointer to this threads random number generator state in global memory
 * @return A pseudo-random integer
 */
double prng_double(prng_state * const restrict state)
{
#ifdef USE_PRNG_NR3
	return nr3_double(state);
#elif defined(USE_PRNG_RANLUX)
	float4 tmp = ranluxcl(state);
	return tmp.x;
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}
#endif

/**
 * Draw a single precision floating point number.
 *
 * @param[in] range Upper bound for the drawn number, nummber will be one less than this at maximum
 * @param[in,out] state Pointer to this threads random number generator state in global memory
 * @return A pseudo-random integer
 */
float prng_float(prng_state * const restrict state)
{
#ifdef USE_PRNG_NR3
	return nr3_float(state);
#elif defined(USE_PRNG_RANLUX)
	float4 tmp = ranluxcl(state);
	return tmp.x;
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}

#ifdef _USEDOUBLEPREC_
/**
 * Draw 4 double precision floating point numbers.
 *
 * Depending on the PRNG this might be more efficient than pulling multiple numbers via seperate calls.
 *
 * @param[in] range Upper bound for the drawn number, nummber will be one less than this at maximum
 * @param[in,out] state Pointer to this threads random number generator state in global memory
 * @return A pseudo-random integer
 */
double4 prng_double4(prng_state * const restrict state)
{
#ifdef USE_PRNG_NR3
	return (double4) (nr3_double(state), nr3_double(state), nr3_double(state), nr3_double(state));
#elif defined(USE_PRNG_RANLUX)
	float4 tmp = ranluxcl(state);
	return (double4) (tmp.x, tmp.y, tmp.z, tmp.w);
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}
#endif /* _USEDOUBLEPREC */

/**
 * Draw 4 single precision floating point numbers.
 *
 * Depending on the PRNG this might be more efficient than pulling multiple numbers via seperate calls.
 *
 * @param[in] range Upper bound for the drawn number, nummber will be one less than this at maximum
 * @param[in,out] state Pointer to this threads random number generator state in global memory
 * @return A pseudo-random integer
 */
float4 prng_float4(prng_state * const restrict state)
{
#ifdef USE_PRNG_NR3
	return (float4) (nr3_float(state), nr3_float(state), nr3_float(state), nr3_float(state));
#elif defined(USE_PRNG_RANLUX)
	return ranluxcl(state);
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}

/**
 * Get PRNG back into a SIMD-friendly state in case different threads requested different amounts
 * of random numbers.
 */
void prng_synchronize(prng_state * const restrict state)
{
#ifdef USE_PRNG_NR3
	// nothing to do for NR3
#elif defined(USE_PRNG_RANLUX)
	ranluxcl_synchronize(state);
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}

/**
 * Get 1,2,3 in random order
 *
 * @param[in,out] rnd Pointer to this threads random number generator state in global memory
 * @return A random permutation of 1,2,3 in x,y,z.
 */
int3 prng_123(prng_state * const restrict state)
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

#ifdef USE_PRNG_NR3
	int3 res;
	res.x = nr3_int64(state) % 3;
	res.y = (res.x + (nr3_int64(state) % 2) + 1) % 3;
	res.z = 3 - res.x - res.y;
	++res;
	return res;
#elif defined(USE_PRNG_RANLUX)
	float4 tmp = ranluxcl(state);
	int3 res;
	res.x = tmp.x * 3.f;
	res.y = (res.x + (int) (tmp.y * 2) + 1) % 3;
	res.z = 3 - res.x - res.y;
	++res;
	return res;
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}
/**
 * Get a normal distributed complex number
 */
hmc_complex inline gaussianNormalPair(prng_state * const restrict rnd)
{
#ifdef USE_PRNG_NR3
	// Box-Muller method, cartesian form, for extracting two independent normal standard real numbers
	hmc_complex tmp;
	hmc_float u1_tmp;
	//CP: if u1 == 1., p will be "inf"
	do {
#ifdef _USE_DOUBLE_PREC_
		u1_tmp = nr3_double(rnd);
#else
		u1_tmp = nr3_float(rnd);
#endif
		if(u1_tmp < 1.) break;
	} while (1 > 0);

	hmc_float u1 = 1.0 - u1_tmp;
	//  hmc_float u2 = 1.0 - nr3_double(rnd);
#ifdef _USE_DOUBLE_PREC_
	hmc_float u2 = nr3_double(rnd);
#else
	hmc_float u2 = nr3_float(rnd);
#endif
	//CP: this is the standard Box-Müller way:
	hmc_float p  = sqrt(- 2.* log(u1));
	//CP: without the 2, one gets sigma = 0.5 rightaway, this is done in tmlqcd
	//hmc_float p  = sqrt(-log(u1));

	//CP: in tmlqcd, sin and cos are interchanged!!
	tmp.re = p * cos(2 * PI * u2);
	tmp.im = p * sin(2 * PI * u2);
	return tmp;
#elif defined(USE_PRNG_RANLUX)
	// TODO update to current ranluxcl!
	//double4 tmp = ranluxcl64norm(state);
	//return (hmc_complex) {tmp.x, tmp.y};

	// LEGEACY CODE
	// Box-Muller method, cartesian form, for extracting two independent normal standard real numbers
	float4 rands = ranluxcl(rnd);

	//CP: if u1 == 1., p will be "inf"
	hmc_float u1 = 1.0 - rands.x;
	//  hmc_float u3 = 1.0 - nr3_double(rnd);
	hmc_float u2 = rands.y;
	//CP: this is the standard Box-Müller way:
	hmc_float p  = sqrt(- 2.* log(u1));
	//CP: without the 2, one gets sigma = 0.5 rightaway, this is done in tmlqcd
	//hmc_float p  = sqrt(-log(u1));

	//CP: in tmlqcd, sin and cos are interchanged!!
	hmc_complex tmp;
	tmp.re = p * cos(2 * PI * u2);
	tmp.im = p * sin(2 * PI * u2);
	return tmp;
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
}
