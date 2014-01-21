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

#ifndef SU3HEATBATHBENCHMARK_H_
#define SU3HEATBATHBENCHMARK_H_

#include "benchmarkExecutable.h"
#include "../physics/lattices/gaugefield.hpp"
#include "../physics/algorithms/su3heatbath.hpp"

class su3heatbathBenchmark : public benchmarkExecutable
{
public:
  su3heatbathBenchmark(int argc, const char* argv[]);

protected:
	/*
	 * Calls the heatbath and overrelax kernels.
	 * Per iteration, the kernel is called with EVEN and ODD parameters.
	 */
  void performBenchmarkForSpecificKernels() override;
};

#endif /* SU3HEATBATHBENCHMARK_H_ */

