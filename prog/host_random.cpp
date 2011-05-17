#include "host_random.h"

#include <cstdio>

inline int random_123 ()
{
	return rnd.int64() % 3 + 1;
}

void random_1_2_3 (int rand[3])
{
	rand[0] = random_123();
	do {
		rand[1] = random_123();
	} while (rand[1] == rand[0]);
	rand[2] = 6 - rand[1] - rand[0];
}

#ifdef __APPLE__
#define CLU_VEC( vec, idx ) (vec)[idx]
#else
#define CLU_VEC( vec, idx ) (vec).s[idx]
#endif

//CP: NR3-PRNG converted for OpenCL by MB
inline cl_ulong nr3_int64( hmc_ocl_ran * state )
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

inline void nr3_init_state( hmc_ocl_ran * state, cl_ulong seed )
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
}

int init_random_seeds(hmc_ocl_ran * hmc_rndarray, char * seedfile, usetimer * timer)
{
	const cl_ulong MAX_SEED = 4101842887655102017L;
	timer->reset();

	FILE* file = fopen( seedfile, "rb" );

	if( ! file ) {
		std::cerr << "Unable to open file " << seedfile << std::endl;
		return HMC_FILEERROR;
	}

	size_t bytes_read = 0;
	for(size_t i_state = 0; i_state < NUMRNDSTATES; ++i_state) {
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
			std::cerr << "Ran out of bytes after initializing " << i_state << " states using " << bytes_read << " bytes." << std::endl;
			return HMC_INVALIDVALUE;
		}

		// we successfully got bytes for this state -> initialize
		nr3_init_state( &hmc_rndarray[i_state], seed );
	}

	timer->add();
	return HMC_SUCCESS;
}

void SU2Update(hmc_float dst [su2_entries], const hmc_float alpha)
{
	hmc_float delta;
	hmc_float a0 ;
	hmc_float eta ;
	do {
		delta = -log(rnd.doub()) / alpha * pow(cos(2. * PI * rnd.doub()), 2.) - log(rnd.doub()) / alpha;
		a0 = 1. - delta;
		eta = rnd.doub();
	} while ( (1. - 0.5 * delta) < eta * eta);
	hmc_float phi = 2.*PI * rnd.doub();
	hmc_float theta = asin(2.*rnd.doub() - 1.);
	dst[0] = a0;
	dst[1] = sqrt(1. - a0 * a0) * cos(theta) * cos(phi);
	dst[2] = sqrt(1. - a0 * a0) * cos(theta) * sin(phi);
	dst[3] = sqrt(1. - a0 * a0) * sin(theta);
}

