#include "generalExecutable.h"

void generalExecutable::printParametersToScreenAndFile()
{
	logger.info() << "## Starting executable: " << ownName;
	print_info_global(parameters);
        print_info_configs_io(parameters);
        print_info_prng_io(parameters);

	outputToFile.open(filenameForLogfile);
	if (outputToFile.is_open()) {
		meta::print_info_global(&outputToFile, parameters);
		meta::print_info_configs_io(&outputToFile, parameters);
		meta::print_info_prng_io(&outputToFile, parameters);
		outputToFile.close();
	} else {
		throw File_Exception(filenameForLogfile);
	}
}

generalExecutable::generalExecutable(int argc, const char* argv[]) : parameters(argc, argv)
{
	totalRuntimeOfExecutable.reset();
	initializationTimer.reset();
	ownName = argv[0];
	filenameForLogfile = meta::createLogfileName(ownName);
	switchLogLevel(parameters.get_log_level());
	printParametersToScreenAndFile();
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

