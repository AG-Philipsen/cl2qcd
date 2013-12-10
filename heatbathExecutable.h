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

#ifndef HEATBATHEXECUTABLE_H_
#define HEATBATHEXECUTABLE_H_

#include "./executables/generationExecutable.h"
#include "physics/algorithms/heatbath.hpp"

class heatbathExecutable: public generationExecutable
{
public:
	heatbathExecutable(int argc, const char* argv[]);

private:
	int overrelaxSteps;

	/*
	 * Thermalize the system using the heatbath algorithm.
	 */
	void thermalizeAccordingToSpecificAlgorithm();

	/*
	 * Generate configurations using the heatbath algorithm.
	 */
	void generateAccordingToSpecificAlgorithm();

	void writeHeatbathLogfile();

	void printParametersToScreenAndFile();

	void setIterationParameters();
};

#endif /* HEATBATHEXECUTABLE_H_ */
