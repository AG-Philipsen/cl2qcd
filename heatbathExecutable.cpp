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

#include "heatbathExecutable.h"

inline heatbathExecutable::heatbathExecutable(int argc, const char* argv[]) :
		generationExecutable(argc, argv)
{
	initializationTimer.reset();
	meta::print_info_heatbath(ownName, parameters);
	writeHeatbathLogfile();
	setIterationParameters();
	initializationTimer.add();
}

inline void heatbathExecutable::writeHeatbathLogfile()
{
	outputToFile.open(filenameForHeatbathLogfile,
			std::ios::out | std::ios::app);
	if (outputToFile.is_open()) {
		meta::print_info_heatbath(ownName, &outputToFile, parameters);
		outputToFile.close();
	} else {
		throw File_Exception(filenameForHeatbathLogfile);
	}
}

void heatbathExecutable::setIterationParameters()
{
	generationExecutable::setIterationParameters();
	generationSteps += parameters.get_heatbathsteps();
	overrelaxSteps = parameters.get_overrelaxsteps();
}

void heatbathExecutable::thermalizeAccordingToSpecificAlgorithm() {
	physics::algorithms::heatbath(*gaugefield, *prng);
}

void heatbathExecutable::generateAccordingToSpecificAlgorithm() {
	physics::algorithms::heatbath(*gaugefield, *prng, overrelaxSteps);
}

int main(int argc, const char* argv[])
{
	try {
		heatbathExecutable heatbathInstance(argc, argv);
		heatbathInstance.generateConfigurations();
	} //try
	//exceptions from Opencl classes
	catch (Opencl_Error& e) {
		logger.fatal() << e.what();
		exit(1);
	} catch (File_Exception& fe) {
		logger.fatal() << "Could not open file: " << fe.get_filename();
		logger.fatal() << "Aborting.";
		exit(1);
	} catch (Print_Error_Message& em) {
		logger.fatal() << em.what();
		exit(1);
	} catch (Invalid_Parameters& es) {
		logger.fatal() << es.what();
		exit(1);
	}

	return 0;

}
