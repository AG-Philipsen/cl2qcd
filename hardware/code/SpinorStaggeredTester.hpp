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

#include "kernelTester.hpp"

#include "../../meta/util.hpp"
#include "../../host_functionality/host_random.h"
#include "../../physics/prng.hpp"
#include "../../physics/lattices/staggeredfield_eo.hpp"
#include "spinors_staggered.hpp"
#include "complex.hpp"


class SpinorStaggeredTester : public KernelTester {
public:
	SpinorStaggeredTester(std::string kernelName, std::string inputfileIn,
			       int numberOfValues = 1, int typeOfComparision = 1);
	//The following constructor is used only for force tests where one inherits from MolecularDynamicsTester
	//and from this class at the same time (to avoid that system, device and parameters are created twice)
	SpinorStaggeredTester(meta::Inputparameters * parameters, const hardware::System * system,
			      hardware::Device * device);
	virtual ~SpinorStaggeredTester();
	
protected:
	//Methods (protected for inheritance resons)
	void fill_with_zero(su3vec * sf_in, int size);
	hmc_float count_sf(su3vec * in, int size);
	hmc_float calc_var(hmc_float in, hmc_float mean);
	hmc_float calc_var_sf(su3vec * in, int size, hmc_float sum);
	hmc_float count_sf_eo(su3vec * sf_in, int size, bool eo);
	su3vec * createSpinorfield(size_t numberOfElements, int seed = 123456);
	su3vec * createSpinorfieldWithOnesAndZerosDependingOnSiteParity();
	su3vec * createSpinorfieldEvenOddWithOnesAndZerosDependingOnSiteParity();
	void calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Plain<su3vec> * in);
	void calcSquarenormEvenOddAndStoreAsKernelResult(const hardware::buffers::SU3vec * in);
	std::string getSpecificInputfile(std::string inputfileIn);
	
	//Members (protected for inheritance resons)
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
	bool allocatedObjects; //Take trace if system, device, inputparameters are allocated
	
	//These methods are used to produce files for the Reference Code (D'Elia et al)
	void print_staggeredfield_to_textfile(std::string outputfile, su3vec * sf);
	void print_staggeredfield_eo_to_textfile(std::string outputfile, su3vec * sf);
	//Utilities methods
	std::vector<hmc_float> reals_from_su3vec(su3vec v);

private:
	su3vec * inputfield;
	void setMembers();
};



#endif // _HARDWARE_CODE_SPINOR_STAGGERED_TESTER_
