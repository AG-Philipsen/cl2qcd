#include "measurementExecutable.h"

measurementExecutable::measurementExecutable(int argc, const char* argv[]) : generalExecutable (argc, argv)
{
	initializationTimer.reset();
	setIterationVariables();
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

	for (iteration; iteration < iterationEnd; iteration += iterationIncrement)
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





