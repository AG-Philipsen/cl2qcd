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
}

void measurementExecutable::initializeGaugefieldAccordingToIterationVariable(int iteration)
{
	currentConfigurationName = meta::create_configuration_name(parameters, iteration);
	gaugefield = new physics::lattices::Gaugefield(*system, *prng, currentConfigurationName);
}

void measurementExecutable::initializeGaugefieldAccordingToConfigurationGivenInSourcefileParameter()
{
	currentConfigurationName = parameters.get_sourcefile();
	gaugefield = new physics::lattices::Gaugefield(*system, *prng);
}

void measurementExecutable::initializeGaugefield(int iteration)
{
	initializationTimer.reset();
	if (parameters.get_read_multiple_configs()) {
		initializeGaugefieldAccordingToIterationVariable(iteration);
	} else {
		initializeGaugefieldAccordingToConfigurationGivenInSourcefileParameter();
	}
	initializationTimer.add();
}

void measurementExecutable::performMeasurements()
{
	performanceTimer.reset();
	logger.trace() << "Perform inversion(s) on device..";

	int iteration = 0;
	for (iteration = iterationStart; iteration < iterationEnd; iteration += iterationIncrement)
	{
		performMeasurementsForSpecificIteration(iteration);
	}
	logger.trace() << "Inversion(s) done";
	performanceTimer.add();
}

void measurementExecutable::performMeasurementsForSpecificIteration(int iteration)
{
	initializeGaugefield(iteration);
	performApplicationSpecificMeasurements();
	prng->save(iteration);
	delete gaugefield;
}





