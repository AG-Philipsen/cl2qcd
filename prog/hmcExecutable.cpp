#include "hmcExecutable.h"

hmcExecutable::hmcExecutable(int argc, const char* argv[]) :
	generationExecutable(argc, argv)
{

}

void hmcExecutable::setIterationParameters()
{
	generationExecutable::setIterationParameters();
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
