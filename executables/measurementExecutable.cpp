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


#include "measurementExecutable.h"
#include "../physics/utilities.hpp"

void measurementExecutable::checkStartconditions()
{
  if(parameters.get_read_multiple_configs() ){
    logger.info() << "To work on multiple configurations, this executable requires the following parameter value(s) to work properly:";
    logger.info() << "startcondition:\tcontinue";
    if(parameters.get_startcondition() != common::start_from_source ) {
      logger.fatal() << "Found wrong startcondition! Aborting..";
      throw Invalid_Parameters("Found wrong startcondition!", "continue", parameters.get_startcondition());
    }
  }
}

measurementExecutable::measurementExecutable(int argc, const char* argv[], std::string parameterSet) : generalExecutable (argc, argv, parameterSet)
{
	initializationTimer.reset();
	setIterationVariables();
	checkStartconditions();
	initializationTimer.add();
}

void measurementExecutable::setIterationVariables()
{
	iterationStart =		(parameters.get_read_multiple_configs()) ? parameters.get_config_read_start() : 0;
	iterationEnd = 			(parameters.get_read_multiple_configs()) ? parameters.get_config_read_end() + 1 : 1;
	iterationIncrement =	(parameters.get_read_multiple_configs()) ? parameters.get_config_read_incr() : 1;
	iteration = iterationStart;
}

void measurementExecutable::initializeGaugefieldAccordingToIterationVariable()
{
	currentConfigurationName = physics::buildCheckpointName(parameters.get_config_prefix(), parameters.get_config_postfix(), parameters.get_config_number_digits(), iteration);
	gaugefield = new physics::lattices::Gaugefield(*system, &(interfacesHandler->getInterface<physics::lattices::Gaugefield>()), *prng, currentConfigurationName);
}

void measurementExecutable::initializeGaugefieldAccordingToConfigurationGivenInSourcefileParameter()
{
	currentConfigurationName = parameters.get_sourcefile();
	gaugefield = new physics::lattices::Gaugefield(*system, &(interfacesHandler->getInterface<physics::lattices::Gaugefield>()), *prng);
}

void measurementExecutable::initializeGaugefield()
{
	initializationTimer.reset();
	if (parameters.get_read_multiple_configs()) {
		initializeGaugefieldAccordingToIterationVariable();
	} else {
		initializeGaugefieldAccordingToConfigurationGivenInSourcefileParameter();
	}
	initializationTimer.add();
}

void measurementExecutable::performMeasurements()
{
	logger.trace() << "Perform inversion(s) on device..";
	
	for (;iteration < iterationEnd; iteration += iterationIncrement)
	{
	  performMeasurementsForSpecificIteration();
	}
	logger.trace() << "Inversion(s) done";
}

void measurementExecutable::performMeasurementsForSpecificIteration()
{
	initializeGaugefield();
	performanceTimer.reset();
	performApplicationSpecificMeasurements();
	prng->saveToSpecificFile(iteration);
	delete gaugefield;
	performanceTimer.add();
}





