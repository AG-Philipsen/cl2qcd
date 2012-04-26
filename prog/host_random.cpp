#include "host_random.h"

#include <cstdio>

/** Seed for the singleton random number generator rnd */
const unsigned long long int seed = 500000;

#ifdef USE_PRNG_NR3
/** Utility 64-bit integer for random number generation */
typedef unsigned long long int Ullong;
/** Utility 32-bit integer for random number generation */
typedef unsigned int Uint;

/** The Random number generator described in Numerical Recipes 3 */
struct Random {

	/** Random number state */
	Ullong u, v, w;

	/**
	 * Initializes the random number generator.
	 *
	 * @param j Seed for the random number generator state
	 */
	Random(Ullong j) : v(4101842887655102017LL), w(1) {
		u = j ^ v;
		int64();
		v = u;
		int64();
		w = v;
		int64();
	}

	/**
	 * Generate a random 64-bit integer.
	 */
	inline Ullong int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17;
		v ^= v << 31;
		v ^= v >> 8;
		w = 4294957665U * (w & 0xffffffff) + (w >> 32);
		Ullong x = u ^ (u << 21);
		x ^= x >> 35;
		x ^= x << 4;
		return (x + v) ^ w;
	}
	/**
	 * Generate a random 64-bit floating point number.
	 */
	inline double doub() {
		return 5.42101086242752217E-20 * int64();
	}
	/**
	 * Generate a random 32-bit integer.
	 */
	inline Uint int32() {
		return (Uint)int64();
	}
};

/** The singleton single-threaded random number generator */
Random rnd(seed);

inline int random_123 ()
{
	return rnd.int64() % 3 + 1;
}

#endif // USE_PRNG_NR3

// FIXME GPU PRNG should also respect selection

#if defined(__APPLE__) && !defined(CL_VERSION_1_1)
#define CLU_VEC( vec, idx ) (vec)[idx]
#else
#define CLU_VEC( vec, idx ) (vec).s[idx]
#endif

//CP: NR3-PRNG converted for OpenCL by MB
inline cl_ulong nr3_int64(prng_state_dev * state )
{
	CLU_VEC(*state, 0) = CLU_VEC(*state, 0) * 2862933555777941757L + 7046029254386353087L;
	CLU_VEC(*state, 1) ^= CLU_VEC(*state, 1) >> 17;
	CLU_VEC(*state, 1) ^= CLU_VEC(*state, 1) << 31;
	CLU_VEC(*state, 1) ^= CLU_VEC(*state, 1) >> 8;
	CLU_VEC(*state, 2) = 4294957665U * (CLU_VEC(*state, 2) & 0xffffffff) + (CLU_VEC(*state, 2) >> 32);
	cl_ulong tmp = CLU_VEC(*state, 0) ^ (CLU_VEC(*state, 0) << 21);
	tmp ^= tmp >> 35;
	tmp ^= tmp << 4;
	return (tmp + CLU_VEC(*state, 1)) ^ CLU_VEC(*state, 2);
}

inline void nr3_init_state(prng_state_dev * const state, const cl_ulong seed )
{
	CLU_VEC(*state, 1) = 4101842887655102017L;
	CLU_VEC(*state, 2) = 1;
	// TODO abort if seed > y
	CLU_VEC(*state, 0) = seed ^ CLU_VEC(*state, 1);
	nr3_int64( state );
	CLU_VEC(*state, 1) = CLU_VEC(*state, 0);
	nr3_int64( state );
	CLU_VEC(*state, 2) = CLU_VEC(*state, 1);
	nr3_int64( state );

	//CP: if one wants to compare exact numbers to the tmlqcd code, one
	//    should use these to have the same seeds for the random numbers
	///@todo take this out in the end
	bool compare_to_tmlqcd = false;
	if(compare_to_tmlqcd) {
		CLU_VEC(*state, 0) = 1234567;
		CLU_VEC(*state, 1) = 8912345;
		CLU_VEC(*state, 2) = 6789123;
	}

}

void init_random_seeds(prng_state_dev * const hmc_rndarray, char const * const seedfile, int const num_rndstates)
{
	const cl_ulong MAX_SEED = 4101842887655102017L;

	FILE * file = fopen( seedfile, "rb" );
	if( !file) {
		logger.debug() << "No random seeds in work directory. Try random seeds from source directory." ;
		std::stringstream file_in_sourcedir;
		file_in_sourcedir << SOURCEDIR << '/' << seedfile ;
		file = fopen( file_in_sourcedir.str().c_str(), "rb" );
		if( ! file ) throw File_Exception(seedfile);
	}

	size_t bytes_read = 0;
	for(size_t i_state = 0; i_state < (size_t)num_rndstates; ++i_state) {
		cl_ulong seed;
		int f_err = 1;

		// read bytes until we find there is an error or we found a
		// a working one
		do {
			f_err = fread( &seed, sizeof( cl_ulong ), 1, file );
			bytes_read += sizeof( cl_ulong );
		} while( f_err == 1 && seed >= MAX_SEED );

		if( f_err != 1 ) {
			bytes_read -= sizeof( cl_ulong ); // the last read was unsuccessfull, but we incremented anyways -> correct that.
			std::stringstream errstr;
			errstr << "Ran out of bytes after initializing " << i_state << " states using " << bytes_read << " bytes.";
			throw Print_Error_Message(errstr.str());
		}

		// we successfully got bytes for this state -> initialize
		nr3_init_state( &hmc_rndarray[i_state], seed );
	}

	return;
}

void prng_init(uint32_t seed)
{
#ifdef USE_PRNG_NR3
	rnd = Random(seed);
#else // USE_PRNG_NR3
#error 'No implemented PRNG chosen'
#endif // USE_PRNG_NR3
}

double prng_double()
{
#ifdef USE_PRNG_NR3
	return rnd.doub();
#else // USE_PRNG_NR3
#error 'No implemented PRNG chosen'
#endif // USE_PRNG_NR3
}

void SU2Update(hmc_float dst [su2_entries], const hmc_float alpha)
{
	hmc_float delta;
	hmc_float a0 ;
	hmc_float eta ;
#ifdef USE_PRNG_NR3
	do {
		delta = -log(rnd.doub()) / alpha * pow(cos(2. * PI * rnd.doub()), 2.) - log(rnd.doub()) / alpha;
		a0 = 1. - delta;
		eta = rnd.doub();
	} while ( (1. - 0.5 * delta) < eta * eta);
	hmc_float phi = 2.*PI * rnd.doub();
	hmc_float theta = asin(2.*rnd.doub() - 1.);
#else // USE_PRNG_NR3
#error 'No implemented PRNG chosen'
#endif // USE_PRNG_NR3
	dst[0] = a0;
	dst[1] = sqrt(1. - a0 * a0) * cos(theta) * cos(phi);
	dst[2] = sqrt(1. - a0 * a0) * cos(theta) * sin(phi);
	dst[3] = sqrt(1. - a0 * a0) * sin(theta);
}

void gaussianComplexVector(hmc_complex * vector, int length, hmc_float sigma)
{

	// SL: this fills real and imaginary part of a vector of "length" complex numbers
	//     with components drawn with a Gaussian distribution and variance sigma
	for(int idx = 0; idx < length; idx++) {
		gaussianNormalPair(&vector[idx].re, &vector[idx].im);
		vector[idx].re *= sigma;
		vector[idx].im *= sigma;
	}
	return;
	// SL: not yet tested
}

void gaussianNormalPair(hmc_float * z1, hmc_float * z2)
{
#ifdef USE_PRNG_NR3
	// Box-Muller method, cartesian form, for extracting two independent normal standard real numbers
	hmc_float u1 = 1.0 - rnd.doub();
	hmc_float u2 = 1.0 - rnd.doub();
#else // USE_PRNG_NR3
#error 'No implemented PRNG chosen'
#endif // USE_PRNG_NR3
	hmc_float p  = sqrt(-2 * log(u1));
	*z1 = p * cos(2 * PI * u2);
	*z2 = p * sin(2 * PI * u2);
	return;
	// SL: not yet tested
}

void random_1_2_3 (int rand[3])
{
#ifdef USE_PRNG_NR3
	rand[0] = random_123();
	do {
		rand[1] = random_123();
	} while (rand[1] == rand[0]);
	rand[2] = 6 - rand[1] - rand[0];
#else // USE_PRNG_NR3
#error 'No implemented PRNG chosen'
#endif // USE_PRNG_NR3
}
