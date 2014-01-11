#include "measurementExecutable.h"

void measurementExecutable::checkStartconditions()
{
  if(parameters.get_read_multiple_configs() ){
    logger.info() << "To work on multiple configurations, this executable requires the following parameter value(s) to work properly:";
    logger.info() << "startcondition:\tcontinue";
    if(parameters.get_startcondition() != meta::Inputparameters::start_from_source ) {
      logger.fatal() << "Found wrong startcondition! Aborting..";
      throw Invalid_Parameters("Found wrong startcondition!", "continue", parameters.get_startcondition());
    }
  }
}

measurementExecutable::measurementExecutable(int argc, const char* argv[], bool enableProfilingIn) : generalExecutable (argc, argv, enableProfilingIn)
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
	currentConfigurationName = meta::create_configuration_name(parameters, iteration);
	gaugefield = new physics::lattices::Gaugefield(*system, *prng, currentConfigurationName);
}

void measurementExecutable::initializeGaugefieldAccordingToConfigurationGivenInSourcefileParameter()
{
	currentConfigurationName = parameters.get_sourcefile();
	gaugefield = new physics::lattices::Gaugefield(*system, *prng);
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
	performanceTimer.reset();
	logger.trace() << "Perform inversion(s) on device..";

	for (;iteration < iterationEnd; iteration += iterationIncrement)
	{
		performMeasurementsForSpecificIteration();
	}
	logger.trace() << "Inversion(s) done";
	performanceTimer.add();
}

void measurementExecutable::performMeasurementsForSpecificIteration()
{
	initializeGaugefield();
	performApplicationSpecificMeasurements();
	prng->save(iteration);
	delete gaugefield;
}





