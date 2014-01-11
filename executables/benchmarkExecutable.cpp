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


#include "benchmarkExecutable.h"


benchmarkExecutable::benchmarkExecutable(int argc, const char* argv[]) : generalExecutable (argc, argv)
{
	initializationTimer.reset();
	benchmarkSteps = parameters.get_benchmarksteps();;
	initializationTimer.add();
}

void benchmarkExecutable::benchmark()
{
	performanceTimer.reset();
	logger.trace() << "Perform benchmarks..";
	for (int iteration = 0; iteration < benchmarkSteps; iteration += 1)
	{
		performBenchmarkForSpecificKernels();
	}
	logger.trace() << "Benchmarks done";
	performanceTimer.add();
}
