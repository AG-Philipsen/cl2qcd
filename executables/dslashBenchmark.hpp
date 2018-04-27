/*
 * Copyright (c) 2014,2018 Alessandro Sciarra
 * Copyright (c) 2014 Christopher Pinke
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

#ifndef DSLASHBENCHMARK_H_
#define DSLASHBENCHMARK_H_

#include "../hardware/code/fermions.hpp"
#include "../hardware/code/spinors.hpp"
#include "../physics/fermionmatrix/fermionmatrix.hpp"
#include "../physics/lattices/spinorfield_eo.hpp"
#include "benchmarkExecutable.hpp"

class dslashBenchmark : public benchmarkExecutable {
  public:
    dslashBenchmark(int argc, const char* argv[]);

  protected:
    const physics::lattices::Spinorfield_eo* spinorfield1;
    const physics::lattices::Spinorfield_eo* spinorfield2;

    /*
     * Calls the dslash_eo kernel.
     * Per iteration, the kernel is called with EVEN and ODD parameters.
     */
    void performBenchmarkForSpecificKernels() override;
    /*
     * Calls dslash_eo on all devices in the system.
     * Per iteration, the kernel is called with EVEN and ODD parameters.
     */
    void enqueueSpecificKernelForBenchmarkingMultipleDevices() override;

    void printProfilingDataToScreen() override;
};

#endif /* DSLASHBENCHMARK_H_ */
