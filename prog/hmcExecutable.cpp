#include "hmcExecutable.h"

hmcExecutable::hmcExecutable(int argc, const char* argv[]) :
	generationExecutable(argc, argv)
{
	initializationTimer.reset();
	meta::print_info_hmc(ownName, parameters);
	writeHmcLogfile();
	setIterationParameters();
	initializationTimer.add();
}

inline void hmcExecutable::writeHmcLogfile()
{
	outputToFile.open(filenameForHmcLogfile,
			std::ios::out | std::ios::app);
	if (outputToFile.is_open()) {
		meta::print_info_heatbath(ownName, &outputToFile, parameters);
		outputToFile.close();
	} else {
		throw File_Exception(filenameForHmcLogfile);
	}
}

void hmcExecutable::setIterationParameters()
{
	generationExecutable::setIterationParameters();
	generationSteps = parameters.get_hmcsteps();
}

void hmcExecutable::thermalizeAccordingToSpecificAlgorithm()
{

}

void hmcExecutable::generateAccordingToSpecificAlgorithm()
{

}

void hmcExecutable::performOnlineMeasurements(int iteration)
{

}

int main(int argc, const char* argv[])
{
	try {
		hmcExecutable hmcInstance(argc, argv);
		hmcInstance.generateConfigurations();
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
