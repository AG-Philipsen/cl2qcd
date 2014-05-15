/** @file
 * SpinorStaggeredTester class
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _HARDWARE_CODE_SPINOR_STAGGERED_TESTER_
#define _HARDWARE_CODE_SPINOR_STAGGERED_TESTER_

#include "../hardware/code/kernelTester.hpp"

#include "../../meta/util.hpp"
#include "../../host_functionality/host_random.h"
#include "../../physics/prng.hpp"
#include "spinors_staggered.hpp"
#include "complex.hpp"


class SpinorStaggeredTester : public KernelTester {
public:
	SpinorStaggeredTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1);
	~SpinorStaggeredTester();
	
protected:
	
	void fill_with_zero(su3vec * sf_in, int size);
	hmc_float count_sf(su3vec * in, int size);
	hmc_float calc_var(hmc_float in, hmc_float mean);
	hmc_float calc_var_sf(su3vec * in, int size, hmc_float sum);
	hmc_float count_sf_eo(su3vec * sf_in, int size, bool eo);
	
	su3vec * createSpinorfield(size_t numberOfElements, int seed = 123456);
	su3vec * createSpinorfieldWithOnesAndZerosDependingOnSiteParity();
	void calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Plain<su3vec> * in);
	void calcSquarenormEvenOddAndStoreAsKernelResult(const hardware::buffers::SU3vec * in);
	std::string getSpecificInputfile(std::string inputfileIn);
	
	const hardware::code::Spinors_staggered * code;
	physics::PRNG * prng;
	hardware::buffers::Plain<double> * doubleBuffer;
	size_t spinorfieldElements;
	size_t spinorfieldEvenOddElements;
	bool useRandom;
	bool evenOrOdd;
	bool calcVariance;
	hmc_complex alpha_host;
	hmc_complex beta_host;
	int iterations;
	
	//Maybe the following should be static functions in the .cpp file
	void fill_with_random(su3vec * in, int size, int seed);
	void fill_with_one(su3vec * in, int size);
	void fill_with_one_eo(su3vec * in, int size, bool eo);
};



#endif // _HARDWARE_CODE_SPINOR_STAGGERED_TESTER_