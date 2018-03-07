/*
 * Copyright (c) 2014,2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DKSBENCHMARK_H_
#define DKSBENCHMARK_H_

#include "../physics/lattices/staggeredfield_eo.hpp"
#include "../hardware/code/fermions_staggered.hpp"
#include "../hardware/code/spinors_staggered.hpp"
#include "../physics/fermionmatrix/fermionmatrix_stagg.hpp"
#include "benchmarkExecutable.hpp"

class dksBenchmark : public benchmarkExecutable {
public:
  dksBenchmark(int argc, const char* argv[]);

protected:
	const physics::lattices::Staggeredfield_eo * staggeredfield1;
	const physics::lattices::Staggeredfield_eo * staggeredfield2;

	/*
	 * Calls the dks_eo kernel.
	 * Per iteration, the kernel is called with EVEN and ODD parameters.
	 */
	void performBenchmarkForSpecificKernels() override;
	/*
	 * Calls dks_eo on all devices in the system.
	 * Per iteration, the kernel is called with EVEN and ODD parameters.
	 */
	void enqueueSpecificKernelForBenchmarkingMultipleDevices()  override;

	void printProfilingDataToScreen() override;
};

#endif /* DSLASHBENCHMARK_H_ */
