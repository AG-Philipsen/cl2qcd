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

#include "su3heatbathExecutable.h"

su3heatbathExecutable::su3heatbathExecutable(int argc, const char* argv[]) :
  generationExecutable(argc, argv)
{
	initializationTimer.reset();
	setIterationParameters();
	printParametersToScreenAndFile();
	initializationTimer.add();
}

void su3heatbathExecutable::printParametersToScreenAndFile()
{
  meta::print_info_observables_gauge_io(parameters);
  meta::print_info_heatbath(parameters);
  writeSu3heatbathLogfile();
}

inline void su3heatbathExecutable::writeSu3heatbathLogfile()
{
	outputToFile.open(filenameForLogfile,
			std::ios::out | std::ios::app);
	if (outputToFile.is_open()) {
		meta::print_info_observables_gauge_io(&outputToFile, parameters);
		meta::print_info_heatbath(&outputToFile, parameters);
		outputToFile.close();
	} else {
		throw File_Exception(filenameForLogfile);
	}
}

void su3heatbathExecutable::setIterationParameters()
{
	generationExecutable::setIterationParameters();
	generationSteps += parameters.get_heatbathsteps();
	overrelaxSteps = parameters.get_overrelaxsteps();
}

void su3heatbathExecutable::thermalizeAccordingToSpecificAlgorithm() {
	physics::algorithms::su3heatbath(*gaugefield, *prng);
}

void su3heatbathExecutable::generateAccordingToSpecificAlgorithm() {
	physics::algorithms::su3heatbath(*gaugefield, *prng, overrelaxSteps);
}

