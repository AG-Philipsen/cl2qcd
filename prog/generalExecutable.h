/*
 * generalExecutable.h
 *
 *  Created on: Dec 1, 2013
 *      Author: christopher
 */

#ifndef GENERALEXECUTABLE_H_
#define GENERALEXECUTABLE_H_

class generalExecutable
{

public:
	generalExecutable(int argc, const char* argv[]) : parameters(argc, argv)
	{
		ownName = argv[0];
		totalRuntimeOfExecutable.reset();
		initializationTimer.reset();
		switchLogLevel(parameters.get_log_level());
		system = new hardware::System(parameters);
		initializationTimer.add();
	}
	~generalExecutable()
	{
		totalRuntimeOfExecutable.add();
		printRuntimeInformationToScreenAndFile();
	}

protected:
	const char* ownName;
	usetimer totalRuntimeOfExecutable;
	usetimer initializationTimer;
	usetimer performanceTimer;
	meta::Inputparameters parameters;
	hardware::System * system;
	ofstream outputToFile;
	const char* generalTimeOutputFilename = "general_time_output";

	void printRuntimeInformationToScreenAndFile()
	{
		printGeneralTimesToScreen();
		printGeneralTimesToFile();
		return;
	}

	void printGeneralTimesToScreen()
	{
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
	void printGeneralTimesToFile()
	{
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
};

#endif /* GENERALEXECUTABLE_H_ */
