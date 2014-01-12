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

/*
 * @file
 * Declaration of the benchmarkExecutable class.
 * This class provides features for the benchmarking of kernels.
 */

#ifndef BENCHMARKEXECUTABLE_H_
#define BENCHMARKEXECUTABLE_H_

#include "generalExecutable.h"

class benchmarkExecutable : public generalExecutable
{
public:
	benchmarkExecutable(int argc, const char* argv[]);

	/**
	 * performs to-be-specified kernels a number of times.
	 * This should be used with profiling enabled.
	 */
	void benchmark();

	/**
	 * Calls a kernel on possibly multiple devices.
	 * The total execution time is measured after a warm-up run.
	 * Therefore, it should not be used with profiling enabled.
	 */
	virtual	void benchmarkMultipleDevices() {};

protected:
  hardware::Device * device;
	int benchmarkSteps;

	virtual void performBenchmarkForSpecificKernels() {};
};

#endif /* BENCHMARKEXECUTABLE_H_ */
