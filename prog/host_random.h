#ifndef _RANDOMH_
#define _RANDOMH_

#include <cstdlib>
#include <cstdio>

#include "hmcerrs.h"
#include "globaldefs.h"
#include "types.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_geometry.h"
#include "host_use_timer.h"


typedef unsigned long long int Ullong;
typedef unsigned int Uint;

const unsigned long long int seed = 500000;

//taken from Numerical recipes
struct Random {
	Ullong u,v,w;
	Random(Ullong j) : v(4101842887655102017LL), w(1) {
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}
	inline Ullong int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	inline double doub() { return 5.42101086242752217E-20 * int64(); }
	inline Uint int32() { return (Uint)int64(); }
};

extern Random rnd;

//get 1,2,3 in random order
void random_1_2_3 (int rand[3]);
//init array with seeds for device
void init_random_seeds(Random random, hmc_ocl_ran * hmc_rndarray, usetimer * inittime);
// Construct new SU2 matrix using improved alg by Kennedy Pendleton
void SU2Update(hmc_float dst [su2_entries], const hmc_float alpha);

#endif
