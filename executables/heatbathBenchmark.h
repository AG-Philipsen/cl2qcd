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

#ifndef HEATBATHBENCHMARK_H_
#define HEATBATHBENCHMARK_H_

#include "benchmarkExecutable.h"
#include "../physics/lattices/gaugefield.hpp"
#include "../physics/algorithms/heatbath.hpp"

class heatbathBenchmark : public benchmarkExecutable
{
public:
  heatbathBenchmark(int argc, const char* argv[]);

protected:
	/*
	 * Calls the heatbath and overrelax kernels.
	 * Per iteration, the kernel is called with EVEN and ODD parameters.
	 */
  void performBenchmarkForSpecificKernels();
};

#endif /* HEATBATHBENCHMARK_H_ */

