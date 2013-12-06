/*
 * generalExecutable.cpp
 *
 *  Created on: Dec 1, 2013
 *      Author: christopher
 */

#include "generalExecutable.h"

generalExecutable::generalExecutable(int argc, const char* argv[]) : parameters(argc, argv)
{
	ownName = argv[0];
	totalRuntimeOfExecutable.reset();
	initializationTimer.reset();
	switchLogLevel(parameters.get_log_level());
	system = new hardware::System(parameters);
	prng = new physics::PRNG(*system);
	initializationTimer.add();
}
generalExecutable::~generalExecutable()
{
	totalRuntimeOfExecutable.add();
	printRuntimeInformationToScreenAndFile();
}

void generalExecutable::printRuntimeInformationToScreenAndFile()
{
	printGeneralTimesToScreen();
	printGeneralTimesToFile();
	return;
}

void generalExecutable::printGeneralTimesToScreen()
{
	using namespace std;
	logger.info() << "## *******************************************************************";
	logger.info() << "## General Times [mus]:";
	logger.info() << "## *******************************************************************";
	logger.info() << "## Program Parts:\t" << setfill(' ') << setw(5) << "total" << '\t' << setw(5) << "perc";
	logger.info() << "## Total:\t" << setfill(' ') << setw(12) << totalRuntimeOfExecutable.getTime();
	logger.info() << "## Init.:\t" << setfill(' ') << setw(12) << initializationTimer.getTime() << '\t' << fixed << setw(5) << setprecision(1) << percent(initializationTimer.getTime(), totalRuntimeOfExecutable.getTime()) ;
	logger.info() << "## Perf.:\t" << setfill(' ') << setw(12) << performanceTimer.getTime() << '\t' << fixed << setw(5) << setprecision(1) << percent(performanceTimer.getTime(), totalRuntimeOfExecutable.getTime()) ;
	logger.info() << "## *******************************************************************";
	return;
}
void generalExecutable::printGeneralTimesToFile()
{
	using namespace std;
	logger.info() << "## writing general times to file: \"" << generalTimeOutputFilename << "\"";
	outputToFile.open(generalTimeOutputFilename);
	if(outputToFile.is_open()) {
		outputToFile  << "## *******************************************************************" << endl;
		outputToFile  << "## General Times [mus]:" << endl;
		outputToFile  << "## Total\tInit\tPerformance" << endl;
		outputToFile  << totalRuntimeOfExecutable.getTime() << "\t" << initializationTimer.getTime() << '\t' << performanceTimer.getTime() << endl;
		outputToFile.close();
	} else {
		logger.warn() << "Could not open output file for general time output.";
	}
	return;
}

void generalExecutable::saveCurrentPrngStateToFile()
{
	logger.info() << "saving current prng state to \"" << filenameForCurrentPrngState << "\"";
	prng->store(filenameForCurrentPrngState);
}

void multipleConfigurationExecutable::setIterationVariables()
{
	iterationStart =		(parameters.get_read_multiple_configs()) ? parameters.get_config_read_start() : 0;
	iterationEnd = 			(parameters.get_read_multiple_configs()) ? parameters.get_config_read_end() + 1 : 1;
	iterationIncrement =	(parameters.get_read_multiple_configs()) ? parameters.get_config_read_incr() : 1;
}

void multipleConfigurationExecutable::initializeGaugefieldAccordingToIterationVariable(int iteration)
{
	currentConfigurationName = meta::create_configuration_name(parameters, iteration);
	gaugefield = new physics::lattices::Gaugefield(*system, *prng, currentConfigurationName);
}

void multipleConfigurationExecutable::initializeGaugefieldAccordingToConfigurationGivenInSourcefileParameter()
{
	currentConfigurationName = parameters.get_sourcefile();
	gaugefield = new physics::lattices::Gaugefield(*system, *prng);
}

void multipleConfigurationExecutable::initializeGaugefield(int iteration)
{
	initializationTimer.reset();
	if (parameters.get_read_multiple_configs()) {
		initializeGaugefieldAccordingToIterationVariable(iteration);
	} else {
		initializeGaugefieldAccordingToConfigurationGivenInSourcefileParameter();
	}
	initializationTimer.add();
}

void multipleConfigurationExecutable::performMeasurements()
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

void multipleConfigurationExecutable::performMeasurementsForSpecificIteration(int iteration)
{
	initializeGaugefield(iteration);
	performApplicationSpecificMeasurements();
	saveCurrentPrngStateToFile();
	delete gaugefield;
}


