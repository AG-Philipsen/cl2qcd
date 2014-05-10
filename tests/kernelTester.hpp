/*
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

#ifndef KERNELTESTER_H_
#define KERNELTESTER_H_

/*
#include <fstream>
#include "../common_header_files/types.h"
#include "../host_functionality/host_use_timer.h"
#include "../executables/exceptions.h"
*/

#include "../physics/prng.hpp"
#include "../hardware/system.hpp"
#include "../physics/lattices/gaugefield.hpp"
#include "../host_functionality/logger.hpp"
#include "../physics/gaugeObservables.h"
//#include "../meta/util.hpp"

class KernelTester
{

public:
  //	virtual ~KernelTester();

	/**
	* Initialize meta::Inputparametes and hardware::System objects
	 */
  KernelTester(std::string kernelNameIn, std::string inputfileIn); 

  std::string kernelName;
  std::string inputfile;
protected:

  meta::Inputparameters * parameters;
  hardware::System * system;
  physics::PRNG * prng;
  physics::lattices::Gaugefield * gaugefield;
};

#endif /* KERNELTESTER_H_ */
